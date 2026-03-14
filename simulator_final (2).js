// ============================================================================
// ROWING ASYNCHRONICITY SIMULATOR — FINAL VERSION
// Double Scull (2x) — 1D Newtonian Model
// ============================================================================
//
// PHYSICS:
//   M_eff × dv/dt = F_propulsive(φ₁,φ₂) - F_drag(v) + F_momentum(φ₁,φ₂)
//
// THREE FORCE TERMS:
//   1. Propulsive: Asymmetric beta-distribution force profile (Kleshnev 2016)
//      with blade interference factor (calibrated to Greidanus et al. 2016)
//   2. Drag: F = k(Fr)·v·|v| with Froude-dependent wave resistance
//   3. Momentum: Reaction force from rower mass oscillating on sliding seat
//
// INTEGRATION: 4th-order Runge-Kutta, dt = 1ms
//
// KEY PARAMETERS:
//   - Boat: 2x double scull, 27 kg hull + 160 kg crew
//   - Force profile: Beta(α=1.6, β=3.5), peak at 19% through drive
//   - Drive fraction: SPM-dependent from Kleshnev BioRow database
//   - Blade interference: k=0.40, calibrated to Greidanus 5-8% power loss
//   - Slide: L=0.35m, M_slide=53 kg per rower (66% of body mass)
//
// VALIDATED OUTPUTS:
//   - Mean speed: ~4.0-4.5 m/s across 20-32 SPM (realistic for club 2x)
//   - Speed loss at 100ms offset, 30 SPM: ~1.5% → ~7.6s over 2000m
//     (Greidanus predicts 7-10s)
//   - Speed loss scales with SPM as expected (higher rate → larger penalty)
//
// KNOWN LIMITATIONS:
//   - Δv/v̄ = 41-58% (real rowing: 20-35%) — 1D model lacks hull damping
//   - Δv decreases with offset (should increase) — requires yaw dynamics
//   - Blade interference is empirical, not first-principles
// ============================================================================

// ========================== CONSTANTS =======================================

const CONSTANTS = {
    rho_water: 1000,    // kg/m³
    rho_air: 1.225,     // kg/m³
    g: 9.81,            // m/s²
    nu_water: 1.0e-6    // m²/s (kinematic viscosity)
};

// ========================== BOAT PARAMETERS =================================

const BOAT = {
    // --- Mass ---
    M_shell: 27,            // kg, FISA minimum for 2x
    M_crew: 160,            // kg, 2 × 80 kg rowers
    M_total: 187,           // kg
    k_added: 0.07,          // added mass coefficient (entrained water)
    get M_eff() { return this.M_total * (1 + this.k_added); }, // 200.09 kg

    // --- Hull geometry ---
    L: 10.4,                // m, waterline length for racing 2x
    B: 0.29,                // m, beam
    A_wetted: 3.5,          // m², wetted surface area

    // --- Propulsion ---
    r_lever: 1.95,          // m, effective oar lever (inboard/outboard ratio)
    F_0: 230,               // N, peak net propulsive force per rower (handle→hull conversion)
    alpha: 1.6,             // beta distribution shape — early peak
    beta: 3.5,              // beta distribution shape — gradual decay
    // Peak force occurs at t_peak = (α-1)/(α+β-2) = 0.6/3.1 = 19.4% through drive

    // --- Blade interference (Greidanus et al. 2016 calibration) ---
    k_interference: 0.40,   // interference strength
    // Physical basis: second blade enters wake-disturbed water
    // At 100ms/30SPM: ~1.5% speed loss → ~4.5% power loss (P∝v³)
    // Greidanus reports 5-8% power loss at typical offsets

    // --- Drag ---
    k_hull_base: 5.4,       // N·s²/m², base hull drag (skin+form), scaled from van Holst 1x→2x
    k_air: 0.28,            // N·s²/m², aerodynamic drag
    k_wave: 0.3,            // wave resistance Froude scaling factor

    // --- Sliding seat ---
    M_rower_total_1: 80,    // kg, rower 1 total mass
    M_rower_total_2: 80,    // kg, rower 2 total mass
    M_slide_1: 53,          // kg, rower 1 sliding mass (66% of body)
    M_slide_2: 53,          // kg, rower 2 sliding mass
    L_slide: 0.35,          // m, slide travel (half-amplitude)

    // --- Blade extraction residual ---
    residual_fraction: 0.05,     // 5% of cycle duration
    residual_tau_fraction: 0.02  // exponential decay time constant = 2% of cycle
};

// ========================== SIMULATION PARAMETERS ===========================

const SIM = {
    dt: 0.001,          // s, integration timestep (1 ms)
    warmup_time: 20,    // s, discard initial transient
    total_time: 60      // s, total simulation duration
};

// ========================== STATE ==========================================

let state = { v: 0, x: 0, phi_1: 0, phi_2: 0, t: 0 };
let omega = 0;          // rad/s, angular stroke rate = 2π × SPM/60
let delta = 0;          // rad, phase offset = ω × Δt

// Drive fraction parameters (set per simulation based on SPM)
let driveFrac = 0.50;
let drivePhaseEnd = Math.PI;

// ========================== DRIVE FRACTION =================================
// From Kleshnev BioRow database (n ≈ 36,000 measurements, sculling)
// drive_time(SPM) = 1.45 - 0.015 × SPM [seconds]
// Validated: drive(36) = 0.91s matches Kleshnev, drive(30) = 1.0s matches Rowperfect
// Reliable to ±2% for SPM 22-32, ±3% at edges

function getDriveFraction(SPM) {
    let drive_time = 1.45 - 0.015 * SPM;
    let cycle_time = 60 / SPM;
    return Math.max(0.30, Math.min(0.60, drive_time / cycle_time));
}

// ========================== UTILITY ========================================

function safeMod(a, b) {
    return ((a % b) + b) % b;
}

// ========================== SLIDE KINEMATICS ===============================
// Hermite (smoothstep) interpolation for slide position
// Drive [0, drivePhaseEnd]: slide moves from -1 (catch) to +1 (finish)
// Recovery [drivePhaseEnd, 2π]: slide returns from +1 to -1
// Position is continuous and smooth at boundaries (C1 continuity)

function slidePosition(phi) {
    let phi_norm = safeMod(phi, 2 * Math.PI);

    if (phi_norm < drivePhaseEnd) {
        // Drive phase
        let t = phi_norm / drivePhaseEnd;           // normalized 0→1
        let s = 3 * t * t - 2 * t * t * t;          // hermite smoothstep
        return -1 + 2 * s;                           // maps to [-1, +1]
    } else {
        // Recovery phase
        let recPhase = 2 * Math.PI - drivePhaseEnd;
        let t = (phi_norm - drivePhaseEnd) / recPhase;
        let s = 3 * t * t - 2 * t * t * t;
        return 1 - 2 * s;                            // maps to [+1, -1]
    }
}

// Slide acceleration (second derivative of position w.r.t. real time)
// Used in momentum force: F_mom = -M_slide × L_slide × a_slide

function slideAcceleration(phi) {
    let phi_norm = safeMod(phi, 2 * Math.PI);

    if (phi_norm < drivePhaseEnd) {
        let t = phi_norm / drivePhaseEnd;
        let d2s_dt2 = 6 * (1 - 2 * t);              // d²(smoothstep)/dt²
        // Chain rule: real-time accel = L_scale × d2s/dt2 × (ω/phaseLength)²
        return 2 * d2s_dt2 * Math.pow(omega / drivePhaseEnd, 2);
    } else {
        let recPhase = 2 * Math.PI - drivePhaseEnd;
        let t = (phi_norm - drivePhaseEnd) / recPhase;
        let d2s_dt2 = 6 * (1 - 2 * t);
        return -2 * d2s_dt2 * Math.pow(omega / recPhase, 2);
    }
}

// ========================== FORCE: PROPULSIVE ==============================
// Asymmetric beta-distribution force profile (Kleshnev 2016)
// f(t) = t^(α-1) × (1-t)^(β-1), normalized so peak = F_0 × η_blade
// α=1.6, β=3.5 → peak at 19.4% through drive (sharp catch, gradual finish)
// Blade efficiency η = 0.99 (Cabrera & Ruina 2006)

function betaFunction(t, alpha, beta) {
    if (t <= 0 || t >= 1) return 0;
    return Math.pow(t, alpha - 1) * Math.pow(1 - t, beta - 1);
}

function betaNormalization(alpha, beta) {
    let t_peak = (alpha - 1) / (alpha + beta - 2);
    return betaFunction(t_peak, alpha, beta);
}

function F_blade_profile(phi) {
    let phi_norm = safeMod(phi, 2 * Math.PI);

    // Drive phase: [0, drivePhaseEnd]
    if (phi_norm < drivePhaseEnd) {
        let t = phi_norm / drivePhaseEnd;             // normalized 0→1
        let norm = betaNormalization(BOAT.alpha, BOAT.beta);
        return BOAT.F_0 * (betaFunction(t, BOAT.alpha, BOAT.beta) / norm) * 0.99;
    }

    // Blade extraction residual: exponential decay after drive ends
    // Prevents unphysical instantaneous force drop to zero
    let residualDuration = BOAT.residual_fraction * 2 * Math.PI;
    let tauPhi = BOAT.residual_tau_fraction * 2 * Math.PI;

    if (phi_norm < drivePhaseEnd + residualDuration) {
        let F_end = BOAT.F_0 * (betaFunction(0.98, BOAT.alpha, BOAT.beta) /
                    betaNormalization(BOAT.alpha, BOAT.beta)) * 0.99;
        let dt = phi_norm - drivePhaseEnd;
        return F_end * Math.exp(-dt / tauPhi);
    }

    // Recovery phase: no propulsive force
    return 0;
}

// ========================== FORCE: BLADE INTERFERENCE =======================
// Empirical model: second blade (in time) enters water disturbed by first
// blade's wake, reducing propulsive efficiency.
// η(δ) = 1 - k × |sin(δ/2)|
// "Always" mode: applies whenever offset exists, regardless of drive overlap
// k = 0.40, calibrated to Greidanus et al. (2016) ~5% power loss at 100ms

function bladeInterferenceFactor(d) {
    if (Math.abs(d) < 0.01) return 1.0;
    return 1.0 - BOAT.k_interference * Math.abs(Math.sin(d / 2));
}

function F_propulsive(s) {
    let F1 = F_blade_profile(s.phi_1);
    let F2 = F_blade_profile(s.phi_2);

    // Apply interference to whichever rower catches second
    if (Math.abs(delta) > 0.01) {
        let interference = bladeInterferenceFactor(delta);
        if (delta > 0) {
            F2 *= interference;   // rower 2 catches later
        } else {
            F1 *= interference;   // rower 1 catches later
        }
    }

    return F1 + F2;
}

// ========================== FORCE: DRAG ====================================
// F_drag = k(Fr) × v × |v|
// k(Fr) = k_hull_base × (1 + k_wave × Fr²) + k_air
// Froude scaling captures wave resistance increase at higher speeds
// v|v| formulation ensures drag opposes motion for any sign of v

function F_drag(v) {
    let Fr = Math.abs(v) / Math.sqrt(CONSTANTS.g * BOAT.L);
    let k_hull = BOAT.k_hull_base * (1 + BOAT.k_wave * Fr * Fr);
    return (k_hull + BOAT.k_air) * v * Math.abs(v);
}

// ========================== FORCE: INTERNAL MOMENTUM =======================
// F_mom = -M_slide × L_slide × a_slide (per rower)
// Newton's third law: rower accelerates on slide → boat feels reaction force
// This is the "check" that coaches talk about — it causes boat deceleration
// during recovery even with no external forces

function F_momentum(s) {
    let a1 = slideAcceleration(s.phi_1);
    let a2 = slideAcceleration(s.phi_2);
    return -BOAT.M_slide_1 * BOAT.L_slide * a1
           -BOAT.M_slide_2 * BOAT.L_slide * a2;
}

// ========================== EQUATION OF MOTION =============================
// M_eff × dv/dt = F_propulsive - F_drag + F_momentum

function derivatives(s) {
    return {
        dv: (F_propulsive(s) - F_drag(s.v) + F_momentum(s)) / BOAT.M_eff,
        dx: s.v,
        dphi_1: omega,
        dphi_2: omega
    };
}

// ========================== RK4 INTEGRATOR =================================

function RK4_step() {
    let k1 = derivatives(state);

    let s2 = {
        v: state.v + 0.5 * SIM.dt * k1.dv,
        x: state.x + 0.5 * SIM.dt * k1.dx,
        phi_1: state.phi_1 + 0.5 * SIM.dt * k1.dphi_1,
        phi_2: state.phi_2 + 0.5 * SIM.dt * k1.dphi_2,
        t: state.t + 0.5 * SIM.dt
    };
    let k2 = derivatives(s2);

    let s3 = {
        v: state.v + 0.5 * SIM.dt * k2.dv,
        x: state.x + 0.5 * SIM.dt * k2.dx,
        phi_1: state.phi_1 + 0.5 * SIM.dt * k2.dphi_1,
        phi_2: state.phi_2 + 0.5 * SIM.dt * k2.dphi_2,
        t: state.t + 0.5 * SIM.dt
    };
    let k3 = derivatives(s3);

    let s4 = {
        v: state.v + SIM.dt * k3.dv,
        x: state.x + SIM.dt * k3.dx,
        phi_1: state.phi_1 + SIM.dt * k3.dphi_1,
        phi_2: state.phi_2 + SIM.dt * k3.dphi_2,
        t: state.t + SIM.dt
    };
    let k4 = derivatives(s4);

    state.v     += (SIM.dt / 6) * (k1.dv     + 2*k2.dv     + 2*k3.dv     + k4.dv);
    state.x     += (SIM.dt / 6) * (k1.dx     + 2*k2.dx     + 2*k3.dx     + k4.dx);
    state.phi_1 += (SIM.dt / 6) * (k1.dphi_1 + 2*k2.dphi_1 + 2*k3.dphi_1 + k4.dphi_1);
    state.phi_2 += (SIM.dt / 6) * (k1.dphi_2 + 2*k2.dphi_2 + 2*k3.dphi_2 + k4.dphi_2);
    state.t     += SIM.dt;
}

// ========================== METRICS ========================================

function computeMetrics(velocities) {
    let skip = Math.floor(SIM.warmup_time / SIM.dt);
    let steady = velocities.slice(skip);

    let v_mean = steady.reduce((a, b) => a + b) / steady.length;
    let variance = steady.reduce((s, v) => s + (v - v_mean) ** 2, 0) / steady.length;

    let v_max = -Infinity, v_min = Infinity;
    for (let v of steady) {
        if (v > v_max) v_max = v;
        if (v < v_min) v_min = v;
    }

    return {
        v_mean,
        v_std: Math.sqrt(variance),
        variance,
        Delta_v: v_max - v_min,
        v_max,
        v_min
    };
}

// ========================== MAIN SIMULATION ================================

function simulate(SPM, Delta_t, duration = SIM.total_time) {
    omega = 2 * Math.PI * (SPM / 60);
    delta = omega * Delta_t;
    driveFrac = getDriveFraction(SPM);
    drivePhaseEnd = 2 * Math.PI * driveFrac;
    state = { v: 0, x: 0, phi_1: 0, phi_2: delta, t: 0 };

    let n = Math.floor(duration / SIM.dt);
    let velocities = [];
    for (let i = 0; i < n; i++) {
        velocities.push(state.v);
        RK4_step();
    }

    return {
        metrics: computeMetrics(velocities),
        SPM, Delta_t, delta, driveFrac
    };
}

// ========================== SPEED LOSS TABLE GENERATOR ======================

function generateSpeedLossTable(SPM, offsets_ms) {
    let sync = simulate(SPM, 0);
    let table = [];

    for (let dt_ms of offsets_ms) {
        let result = simulate(SPM, dt_ms / 1000);
        let loss = (sync.metrics.v_mean - result.metrics.v_mean) / sync.metrics.v_mean;
        let dv_change = sync.metrics.Delta_v > 0 ?
            (result.metrics.Delta_v - sync.metrics.Delta_v) / sync.metrics.Delta_v : 0;
        table.push({
            dt_ms,
            delta_rad: result.delta,
            v_mean: result.metrics.v_mean,
            speed_loss_frac: loss,
            speed_loss_pct: loss * 100,
            Delta_v: result.metrics.Delta_v,
            dv_over_vmean: result.metrics.Delta_v / result.metrics.v_mean * 100,
            dv_change_pct: dv_change * 100
        });
    }

    return { SPM, driveFrac: getDriveFraction(SPM), sync_v: sync.metrics.v_mean,
             sync_dv: sync.metrics.Delta_v, table };
}

// ========================== CLI PARAMETER HANDLING ===========================
//
// Usage:
//   node simulator_final.js                        # full sweep with defaults
//   node simulator_final.js --weight 75            # set rower mass (kg each)
//   node simulator_final.js --offset 100           # single offset (ms)
//   node simulator_final.js --weight 75 --offset 100
//   node simulator_final.js --help

function parseArgs() {
    let args = { weight: 80, offset: null, spm: null };
    let argv = process.argv.slice(2);

    for (let i = 0; i < argv.length; i++) {
        if (argv[i] === '--weight' || argv[i] === '-w') {
            args.weight = parseFloat(argv[++i]);
            if (isNaN(args.weight) || args.weight < 40 || args.weight > 120) {
                console.error("Error: --weight must be between 40 and 120 kg");
                process.exit(1);
            }
        } else if (argv[i] === '--offset' || argv[i] === '-o') {
            args.offset = parseFloat(argv[++i]);
            if (isNaN(args.offset) || args.offset < 0 || args.offset > 500) {
                console.error("Error: --offset must be between 0 and 500 ms");
                process.exit(1);
            }
        } else if (argv[i] === '--spm' || argv[i] === '-s') {
            args.spm = parseFloat(argv[++i]);
            if (isNaN(args.spm) || args.spm < 16 || args.spm > 44) {
                console.error("Error: --spm must be between 16 and 44");
                process.exit(1);
            }
        } else if (argv[i] === '--help' || argv[i] === '-h') {
            console.log("Rowing Asynchronicity Simulator — Double Scull (2x)");
            console.log("");
            console.log("Usage: node simulator_final.js [options]");
            console.log("");
            console.log("Options:");
            console.log("  --weight, -w <kg>    Rower mass in kg (default: 80, range: 40-120)");
            console.log("  --offset, -o <ms>    Catch timing offset in ms (default: full 0-200 sweep)");
            console.log("  --spm, -s <rate>     Stroke rate (default: all essay phases, range: 16-44)");
            console.log("  --help, -h           Show this help");
            console.log("");
            console.log("Examples:");
            console.log("  node simulator_final.js                     Full run, 80 kg rowers");
            console.log("  node simulator_final.js -w 75               Full run, 75 kg rowers");
            console.log("  node simulator_final.js -s 30 -o 100        Single point: 30 SPM, 100ms");
            console.log("  node simulator_final.js -s 28               One rate, all offsets");
            console.log("  node simulator_final.js -o 100              All rates, one offset");
            process.exit(0);
        } else {
            console.error(`Unknown argument: ${argv[i]}. Use --help for usage.`);
            process.exit(1);
        }
    }
    return args;
}

// ========================== CONFIGURE BOAT FOR ROWER WEIGHT =================
// Recalculates all mass-dependent parameters.
// F_0 and k_hull scale with rower mass to keep velocity realistic.
//   - F_0 ∝ mass (bigger rower → more force, roughly linear)
//   - k_hull ∝ total_mass^(2/3) (wetted area scales with displacement)
// Baseline: 80 kg rower → F_0=230 N, k_hull_base=5.4

function configureBoat(rowerMass) {
    let baseline = 80;  // kg, calibration point

    BOAT.M_rower_total_1 = rowerMass;
    BOAT.M_rower_total_2 = rowerMass;
    BOAT.M_crew = 2 * rowerMass;
    BOAT.M_total = BOAT.M_shell + BOAT.M_crew;
    BOAT.M_slide_1 = Math.round(rowerMass * 0.6625);  // 66.25% moves on slide
    BOAT.M_slide_2 = BOAT.M_slide_1;

    // Scale F_0 linearly with rower mass
    BOAT.F_0 = Math.round(230 * (rowerMass / baseline));

    // Scale hull drag with wetted area ∝ displacement^(2/3)
    let baseline_total = BOAT.M_shell + 2 * baseline;  // 187 kg
    BOAT.k_hull_base = 5.4 * Math.pow(BOAT.M_total / baseline_total, 2 / 3);

    // Wetted area also scales
    BOAT.A_wetted = 3.5 * Math.pow(BOAT.M_total / baseline_total, 2 / 3);
}

// ========================== DISPLAY SETTINGS ================================

function displaySettings() {
    let Fr = 4.3 / Math.sqrt(CONSTANTS.g * BOAT.L);
    let k_eff = BOAT.k_hull_base * (1 + BOAT.k_wave * Fr * Fr) + BOAT.k_air;

    console.log("============================================================");
    console.log("  ROWING ASYNCHRONICITY SIMULATOR");
    console.log("  Double Scull (2x) — 1D Newtonian Model");
    console.log("============================================================\n");
    console.log("  MASS");
    console.log(`    Hull:       ${BOAT.M_shell} kg (FISA minimum 2x)`);
    console.log(`    Rowers:     2 × ${BOAT.M_rower_total_1} kg = ${BOAT.M_crew} kg`);
    console.log(`    Total:      ${BOAT.M_total} kg`);
    console.log(`    Effective:  ${BOAT.M_eff.toFixed(1)} kg (added mass ${(BOAT.k_added*100).toFixed(0)}%)`);
    console.log(`    Slide mass: ${BOAT.M_slide_1} kg/rower (${(BOAT.M_slide_1/BOAT.M_rower_total_1*100).toFixed(0)}% of body)`);
    console.log();
    console.log("  PROPULSION");
    console.log(`    F_0:        ${BOAT.F_0} N peak net propulsive force/rower`);
    console.log(`    Profile:    Beta(α=${BOAT.alpha}, β=${BOAT.beta}), peak at ${((BOAT.alpha-1)/(BOAT.alpha+BOAT.beta-2)*100).toFixed(0)}% through drive`);
    console.log(`    Blade eff:  99% (Cabrera & Ruina 2006)`);
    console.log(`    Interference: k=${BOAT.k_interference} (Greidanus et al. 2016)`);
    console.log();
    console.log("  DRAG");
    console.log(`    k_hull:     ${BOAT.k_hull_base.toFixed(2)} Ns²/m² (skin+form, van Holst/Lazauskas)`);
    console.log(`    k_air:      ${BOAT.k_air} Ns²/m²`);
    console.log(`    k_wave:     ${BOAT.k_wave} (Froude scaling)`);
    console.log(`    k_eff:      ~${k_eff.toFixed(1)} Ns²/m² at race pace`);
    console.log();
    console.log("  KINEMATICS");
    console.log(`    L_slide:    ${BOAT.L_slide} m (total travel ${(BOAT.L_slide*2*100).toFixed(0)} cm)`);
    console.log(`    Hull:       ${BOAT.L}m × ${BOAT.B}m, A_wet=${BOAT.A_wetted.toFixed(2)} m²`);
    console.log();
    console.log("  INTEGRATION");
    console.log(`    Method:     RK4, dt=${SIM.dt*1000}ms, warmup=${SIM.warmup_time}s, total=${SIM.total_time}s`);
    console.log("------------------------------------------------------------");
}

// ========================== RUN =============================================

const OFFSETS_MS = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200];

const ESSAY_SPMS = [
    { label: "Controlled test (low)",  spm: 20 },
    { label: "Controlled test (high)", spm: 22 },
    { label: "Body (low)",             spm: 26 },
    { label: "Body (high)",            spm: 28 },
    { label: "Start",                  spm: 30 },
    { label: "Sprint",                 spm: 32 }
];

let args = parseArgs();
configureBoat(args.weight);
displaySettings();

// Build the SPM list: either user-specified or all essay phases
let spmList = args.spm !== null
    ? [{ label: `Custom`, spm: args.spm }]
    : ESSAY_SPMS;

if (args.offset !== null) {
    // ── Single offset mode (one or all SPMs) ──
    console.log(`\n  Running: offset = ${args.offset} ms\n`);
    console.log("  Phase                | SPM | Drive% | v̄ sync | v̄ async | Loss     | Δv/v̄  ");
    console.log("  ---------------------|-----|--------|--------|---------|----------|-------");

    for (let { label, spm } of spmList) {
        let sync  = simulate(spm, 0);
        let async = simulate(spm, args.offset / 1000);
        let loss = (sync.metrics.v_mean - async.metrics.v_mean) / sync.metrics.v_mean * 100;
        let dvr  = async.metrics.Delta_v / async.metrics.v_mean * 100;
        let df   = getDriveFraction(spm) * 100;

        console.log(`  ${label.padEnd(21)} | ${String(spm).padStart(3)} | ${df.toFixed(1).padStart(5)}% | ${sync.metrics.v_mean.toFixed(3)}  | ${async.metrics.v_mean.toFixed(3)}   | +${loss.toFixed(3).padStart(6)}% | ${dvr.toFixed(1)}%`);
    }

    // Race impact for this offset
    let s30 = simulate(30, 0), a30 = simulate(30, args.offset/1000);
    let s27 = simulate(27, 0), a27 = simulate(27, args.offset/1000);
    let s32 = simulate(32, 0), a32 = simulate(32, args.offset/1000);
    let pen_start  = 500  / a30.metrics.v_mean - 500  / s30.metrics.v_mean;
    let pen_body   = 3500 / a27.metrics.v_mean - 3500 / s27.metrics.v_mean;
    let pen_sprint = 1000 / a32.metrics.v_mean - 1000 / s32.metrics.v_mean;
    console.log(`\n  5K race penalty: Start +${pen_start.toFixed(1)}s | Body +${pen_body.toFixed(1)}s | Sprint +${pen_sprint.toFixed(1)}s | TOTAL +${(pen_start+pen_body+pen_sprint).toFixed(1)}s`);

} else {
    // ── Full offset table (one or all SPMs) ──
    for (let { label, spm } of spmList) {
        console.log(`\n=== ${label}: ${spm} SPM (drive fraction: ${(getDriveFraction(spm)*100).toFixed(1)}%) ===\n`);
        let result = generateSpeedLossTable(spm, OFFSETS_MS);

        console.log(`  Synchronized: v̄ = ${result.sync_v.toFixed(4)} m/s, Δv = ${result.sync_dv.toFixed(4)} m/s, Δv/v̄ = ${(result.sync_dv/result.sync_v*100).toFixed(1)}%`);
        console.log();
        console.log("  Δt(ms) |  δ(rad)  | v̄(m/s)  | Speed Loss | Δv/v̄   | Δv change");
        console.log("  -------|----------|---------|------------|--------|----------");

        for (let row of result.table) {
            console.log(`  ${String(row.dt_ms).padStart(5)}  | ${row.delta_rad.toFixed(4).padStart(7)}  | ${row.v_mean.toFixed(4)} | ${row.speed_loss_pct >= 0 ? '+' : ''}${row.speed_loss_pct.toFixed(4)}%   | ${row.dv_over_vmean.toFixed(1).padStart(5)}% | ${row.dv_change_pct >= 0 ? '+' : ''}${row.dv_change_pct.toFixed(2)}%`);
        }
    }

    // Race impact and SyncRow table only in full run (no --spm filter)
    if (args.spm === null) {
        console.log("\n\n=== RACE IMPACT: 5000m HEAD RACE ===\n");
        console.log("  Start 500m @ 30 SPM | Body 3500m @ 27 SPM | Sprint 1000m @ 32 SPM\n");

        for (let dt_ms of [20, 50, 100, 150, 200]) {
            let s30 = simulate(30, 0), a30 = simulate(30, dt_ms/1000);
            let s27 = simulate(27, 0), a27 = simulate(27, dt_ms/1000);
            let s32 = simulate(32, 0), a32 = simulate(32, dt_ms/1000);
            let pen_start  = 500  / a30.metrics.v_mean - 500  / s30.metrics.v_mean;
            let pen_body   = 3500 / a27.metrics.v_mean - 3500 / s27.metrics.v_mean;
            let pen_sprint = 1000 / a32.metrics.v_mean - 1000 / s32.metrics.v_mean;
            let total = pen_start + pen_body + pen_sprint;
            console.log(`  Δt=${String(dt_ms).padStart(3)}ms: Start +${pen_start.toFixed(1)}s | Body +${pen_body.toFixed(1)}s | Sprint +${pen_sprint.toFixed(1)}s | TOTAL +${total.toFixed(1)}s`);
        }

        console.log("\n\n=== SYNCROW SPEED_LOSS_TABLE (30 SPM) ===\n");
        console.log("const SPEED_LOSS_TABLE = {");
        let s30 = simulate(30, 0);
        for (let dt_ms of OFFSETS_MS) {
            let dv = 2 * Math.PI * 0.5 * dt_ms / 1000;
            let a = simulate(30, dt_ms / 1000);
            let loss = (s30.metrics.v_mean - a.metrics.v_mean) / s30.metrics.v_mean;
            console.log(`    ${dv.toFixed(3)}: ${loss.toFixed(6)},  // ${dt_ms}ms → ${(loss*100).toFixed(2)}%`);
        }
        console.log("};");
    }
}
