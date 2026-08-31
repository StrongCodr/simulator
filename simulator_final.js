// ============================================================================
// ROWING ASYNCHRONICITY SIMULATOR
// Double Scull (2x) — 1D Newtonian model with emergent blade mechanics
// ============================================================================
//
// PURPOSE: quantify the mean-speed cost of inter-rower CATCH-TIMING OFFSETS
// (Δt) in a 2x, as a physics companion to the SyncRow measurement program.
//
// PHYSICS (equation of motion — two-mass formulation, derivation at the
// EQUATION OF MOTION section):
//
//   M_total(1+k_added) · dv/dt = Σᵢ F_blade,x,i − F_drag(v) + F_momentum
//
//   • v is the HULL velocity; drag acts on it, not on the system CM.
//   • F_blade,x,i — blade-slip closure: each rower's handle-force profile is
//     PRESCRIBED (Beta(2.92, 3.08) in drive time, fitted BY US to Kleshnev's
//     published shape statistics — see the profile block); the oar is a
//     massless lever, so blade-side perpendicular force is F_h·L_in/L_out,
//     and the oar angular velocity is SOLVED each instant so Caplan-Gardner
//     blade hydrodynamics (Cl = Al·sin2α, Cd = Ad·sin²α; Big Blade
//     Al = 1.25, Ad = 2.07 via Atkinson's distillation) supply exactly that
//     force. Blade slip, oar angle, swept arc, release angle and blade
//     efficiency are EMERGENT, and blade force depends on the boat's
//     instantaneous speed when each blade is loaded — so ALL inter-rower
//     coupling is emergent through the shared hull (formulation lineage:
//     van Holst; Atkinson's ROWING; Cabrera-Ruina-Kleshnev 2006).
//   • F_drag: ½·ρ·Cf(Re)·S·v|v|·(1+k_wave·Fr²) + k_air·v|v|, with
//     S = 2.30·√(V·Lwl) computed from displacement each run and ITTC-1957
//     Cf(Re). There is NO drag constant — only cited, measurable quantities.
//   • F_momentum — per-mass body model: −Σ mⱼ·s̈ⱼ over per-rower body masses
//     (de Leva fractions) on quintic (C²) slide motion with per-mass travels,
//     seat leading the handle by the Catch Factor (−15 ms); arm mass rides
//     the emergent handle path.
//
// CALIBRATION PROTOCOL: ONE free scale — EFFORT.F_peak_ref, the
// peak total handle force per 80 kg rower — set by running `--calibrate`,
// which bisects F_peak so the synchronized 2x at 30 SPM reaches the club-pace
// anchor 4.30 m/s [ASSUMED anchor]. Current value 352.0 N [DERIVED, 31 Aug
// 2026, after the sculling catch-angle correction].
// Every output regenerates from this one number.
//
// EXPERIMENT DESIGN: offset comparisons hold PER-ROWER AVERAGE HANDLE
// POWER equal to the synchronized baseline (iterative effort renormalization),
// so the reported penalty is a TIMING cost, not a hidden power cut. --rawpower
// compares at fixed force profiles instead. Steady-state phase speeds only;
// standing starts are not modeled.
//
// VERIFICATION: startup self-checks (sign convention through the live force
// path; closure torque balance), a per-cycle energy audit (W_prop + W_mom −
// W_drag = ΔKE within 0.001% of W_drag — the gate was tightened 100× after the
// 30 Aug 2026 code review drove residuals to ~1e-4%; loud failure), and a
// BLINDED output table: the program prints our numbers ONLY. Literature
// reference values are kept in the accompanying fix specification
// (simulator_fix_spec.md), outside this program, and are consulted once, at
// the end, never at runtime. No output of this program states a literature
// value.
//
// VALIDATION STATUS (end-of-work literature comparison, 30-31 Aug 2026; the
// model was frozen before any of these numbers were looked up):
//   • Comparative behaviour MATCHES on-water measurement: async is smoother
//     (swing 31.6% sync → 3.9% antiphase; de Brouwer, de Poel & Hofmijster
//     2013; Cuijpers, Zaal & de Poel 2015; Boucher, Labbé & Clanet 2017) and
//     sync is FASTER (antiphase −0.57% raw / −2.98% equal-power at 30 SPM;
//     Boucher et al.'s robot crews on water). The model takes the on-water
//     side of the de Brouwer (slides ergometer) vs Boucher (water) disagreement,
//     via the body-motion/hull-velocity phasing term that antiphase cancels.
//   • Stroke-rate dependence MATCHES: fixed peak handle force is what Kleshnev
//     (2021 stroke-rate study, n≈36k) measured below 32 spm ("average force
//     nearly constant"); model speed 20→32 spm +16.4% vs his handle-velocity
//     trend implying ~+17%; blade efficiency rises and velocity efficiency
//     falls with rate, as in Hofmijster et al. 2007. Not modelled: his force
//     shape drifting toward rectangular (+4-6% avg/max over 20-44 spm) and the
//     3.5% arc shortening — sub-percent within 20-32 spm.
//   • Force-peak position 39.3% of drive length vs Kleshnev's ≤40% target /
//     38.2% racing mean (2011, Table 9.2) — inside, after the sculling catch-
//     angle correction of 31 Aug 2026.
//   • Absolute misses, ALL in one direction (model too idealized), with the
//     root cause quantified against Kleshnev's published acceleration curves
//     (RBN 2023/01, Fig. 3): the model's catch deceleration is −2.5 m/s² at
//     30 spm and −4.0 at 36 spm, vs measured −6 (JM8+, 36 spm) to −12 (M1x,
//     37 spm); the drive positive peak (+3.4/+4.8) matches his +3 to +5. The
//     CHECK is 2-3× too shallow, by construction: blade force is zero at the
//     catch (no oar inertia — Kleshnev 2024 averaged curves: 13% of max at
//     the catch) and the body reversal is smoothed. A shallow check raises
//     v_min, so: swing 31% vs ~40% implied for singles (peak-to-peak; note
//     Greidanus et al. 2016's "A = 20-25%" is the HALF-amplitude), fluctuation
//     loss 3.1% vs 5-6% (de Brouwer et al. 2013; Kleshnev 1999 velocity
//     efficiency 93.5-95.5%), blade efficiency 0.89 vs 0.78 ± 0.06 (Kleshnev
//     1999, 21 crews; Hofmijster et al. 2007/2010 losses 20-30%), release
//     49° vs 42-44° elite sculling. One missing piece of physics — oar inertia
//     and a real catch — moves all four the right way. The 2x-specific
//     reference for swing and catch deceleration is not pinned (bracketed by
//     1x and 8+); BioRow's 720-coefficient acceleration model would give it.
//
// INTEGRATION: RK4, dt = 1 ms (--dt override); warmup to stroke-mean
// convergence (<0.001 m/s, cap 400 strokes, loud warning); measurement =
// exactly 20 phase-anchored complete cycles, window opened and closed through
// the same sub-timestep boundary interpolation; means computed from distance
// over exact cycle windows. Numerical floor ±0.0001 pp on losses (convergence
// study 30 Aug 2026, printed under every table); headline dt-invariant.
//
// REPORTING: every equal-power loss ships with its raw-power (fixed force
// profile) companion — they differ 6-12× at 100 ms, and the ratio moves with
// stroke rate, so neither number is meaningful alone.
//
// HULL ATTITUDE (the hull-attitude modes, opt-in, default off): --trim = quasi-static
// trim (Formaggia static limit; wetted area first-order blind to trim on a
// symmetric wall-sided hull — an honest null), --pitch = heave+pitch dynamics
// with pitch radiation reported as a diagnostic, not fed back into surge.
//
// PROVENANCE CONVENTION: every parameter carries exactly one label —
// [CITED source] / [DERIVED — inputs cited, arithmetic shown] / [ASSUMED].
// PLACEHOLDER marks values SyncRow measurements will replace.
//
// KNOWN LIMITATIONS (excluded channels, honest mapping — the header rewrite):
//   - Pitch/trim (vertical plane): EXCLUDED. Maps to exactly our variable —
//     inter-rower offsets shift fore-aft mass phasing between rowers seated
//     ~2.5 m apart, exciting different pitch forcing. Formaggia's group puts
//     secondary heave/pitch at up to ~10% of total dissipation — larger than
//     the effects under study. The one scope expansion with a number behind
//     it (the recorded heave/pitch extension decision; cheap static-trim variant: the static-trim mode).
//   - Yaw/sway (horizontal plane): EXCLUDED. Maps to WITHIN-rower bilateral
//     asymmetry and sweep boats, not to our inter-rower offset in a scull.
//   - Puddle/wake re-entry between rowers: NOT modeled — no published model
//     or data exists (the deleted k_interference pretended otherwise; the interference-term deletion).
//   - Quasi-static blade coefficients: transient/unsteady blade forces differ
//     (Grift et al. 2019); Caplan & Gardner 2007a flume values may also run ~6%
//     high (Sliasas & Tullis 2009 shallow-flume re-analysis) — swept in the sensitivity sweep.
//   - OAR INERTIA AND THE CATCH: EXCLUDED — DECISION RECORDED 31 Aug 2026
//     (Vess): documented, not built. The oar is massless, so the quasi-static
//     torque balance jumps the blade to its zero-force angle the instant the
//     pulse opens; a real oar lags, and during the lag the planted blade brakes
//     the hull — that lag IS the catch check. This one exclusion is the
//     identified common cause of all four absolute misses in VALIDATION STATUS
//     (catch deceleration 2-3× too shallow → v_min too high → swing 31% vs
//     ~40%, fluctuation loss 3.1% vs 5-6%, blade efficiency 0.89 vs 0.78,
//     release 49° vs 42-44°). Also the source of Kleshnev's 13%-of-max handle
//     force at zero handle velocity (2024 averaged curves).
//     WHY IT IS NOT BUILT: the code change is modest — replace the algebraic
//     closure with I·θ̈ = F_h·L_in − F_blade,⊥·L_out per rower, sub-stepped
//     because the hydrodynamic damping makes the oar equation stiff (relaxation
//     ≪ 1 ms); it would also remove the catch-instant discontinuity and the
//     arm-impulse artifact. The PHYSICS is the obstacle: the oar shaft's own
//     inertia (~1.5 kg scull, ~2 kg·m² about the pin; Laschowski et al. 2015)
//     is ~50× smaller than the hydrodynamic ADDED mass of the blades
//     (potential-flow plate estimate ~17 kg per 46 × 21.5 cm blade, ~100 kg·m²
//     for the pair at L_out), and that added mass is NOT a constant during the
//     entry acceleration (Grift et al. 2019; Sliasas & Tullis 2010 needed full
//     unsteady CFD). Building it now means either an [ASSUMED] coefficient
//     known to be wrong in the relevant regime, or a coefficient FITTED to
//     Kleshnev's measured catch deceleration — a declared second calibration
//     that would turn swing, fluctuation loss and blade efficiency from
//     independent validations into consequences of a fit. Neither improves
//     credibility. WHEN TO BUILD IT: when SyncRow's oar-angle-vs-time record
//     gives θ̈ through the catch, so the effective inertia is fitted from the
//     boat being studied, not from a plate formula or a literature target.
//     Keep the massless closure as a switchable comparison mode when that day
//     comes.
//   - Oar-shaft flexibility: neglected (small energy store, ~few %).
//   - Shallow water: excluded; amplifies fluctuation losses where depth < ~3 m.
//   - Standing-start dynamics: excluded; all results are steady-state.
//   - Waterline vs overall length: overall 10.0 m used as Lwl (stated, the drag rebuild).
// ============================================================================

// ========================== CONSTANTS ======================================

const CONSTANTS = {
    rho_water: 997,     // kg/m³ [CITED: Sliasas & Tullis 2009 — rowing blade/hull hydrodynamics; changed from 1000 at the drag rebuild]
    rho_air: 1.225,     // kg/m³
    g: 9.81,            // m/s²
    nu_water: 1.0e-6    // m²/s (kinematic viscosity)
};

// ========================== BOAT PARAMETERS ================================

// REFERENCE ROWER (stated once; all geometric baselines below are quoted for
// this rower and anything scaled for other body sizes scales from here):
// a 6 ft (183 cm), 80 kg male sculler. Kleshnev notes that standard textbook
// rigging dimensions suit a sculler of ~190 cm (Kleshnev 2007, BioRow newsletter No. 3), so 183 cm
// sits just inside the population the published averages describe.

const BOAT = {
    // --- Mass: MOVED (the single mass block) ---
    // Every mass, mass fraction, travel distance, and the added-mass
    // coefficient now lives EXCLUSIVELY in the MASS block below — the single
    // source of truth. The provenance essays that sat here (Zatsiorsky/de Leva
    // 71.8% derivation; the 0.52 m route-3 travel decision and the retired
    // routes; the k_added anchors) are preserved in the MASS block and in
    // the accompanying fix specification (simulator_fix_spec.md). Effective inertia (added water) is computed in the
    // EQUATION OF MOTION from MASS.total and MASS.k_added.

    // --- Hull geometry (the drag rebuild — every value cited, all values CONSUMED) ---
    L: 10.0,                // m [DOUBLE-CITED: Filippi Lido F17 double-scull specification = 10.0 m (85-105 kg
                            // athletes; "the ideal" pro 2x hull); full Filippi 2x range
                            // 9.40 m (F13, 65-75 kg) to 10.14 m (F51, 85-100 kg), keyed
                            // to crew weight — 80 kg crew sits nearest the F17. The old
                            // 10.4 m EXCEEDED the largest double Filippi builds and
                            // traced to a general-audience aggregator, not a manufacturer specification.
                            // World Rowing does not regulate class length (only min
                            // weight, 2x = 27 kg, Bye-Law to Rule 32).
                            // ASSUMPTION STATED: this is OVERALL length used as
                            // waterline length; true Lwl is somewhat shorter.]
    B: 0.372,               // m [DOUBLE-CITED: Filippi F17 (Filippi Lido specifications) = 37.2 cm; F13 = 33.0 cm.
                            // The old 0.29 m matched a SINGLE (Filippi F01 1x = 26.6 cm)
                            // — a boat-class error, unnoticed while beam was unused.]
    S_coeff: 2.30,          // wetted-area relation S = S_coeff·√(V·Lwl) [DERIVED from
                            // three measured racing hulls — US Patents 5474008 (8+),
                            // 5188048 (4-), 5279239 (LW4), the Vespoli-lineage hulls
                            // Scragg measured and Atkinson's ROWING is built on:
                            // individual coefficients 2.168 / 2.345 / 2.377, mean 2.30,
                            // spread <10% across 331-914 kg. Cross-check: displacement^
                            // (2/3) scaling of Sliasas & Tullis 2009's 4- (5.0 m² at 426 kg,
                            // Kleshnev pers. comm.) gives 2.89 m² for the 2x vs 3.15
                            // computed — within 7%.]
    Cf_const: 0.00225,      // constant skin-friction coefficient for the --constcf
                            // comparison mode [CITED: Sliasas & Tullis 2009]

    // --- Handle-force profile shape (the closure and the force-profile refit) ---
    // Beta shape in DRIVE TIME for the PRESCRIBED handle force (the prescribed-force closure). The old
    // prescribed-propulsive F_0 = 230 N and the uncited r_lever = 1.95 were
    // DELETED with the blade-slip closure; the effort scale lives in EFFORT.
    //
    // [The force-profile refit — parameters refit 30 Aug 2026, superseding alpha=1.6/beta=3.5
    // whose "Kleshnev 2016" attribution was false (a real book — Kleshnev 2016, see REFERENCES — but it contains no beta
    // parameterization in it — Kleshnev publishes measured curves and targets,
    // not fitted distributions).]
    // Beta parameterization fitted BY US to Kleshnev's published shape
    // statistics: peak force at 40% of drive LENGTH (BioRow, "Averaged
    // biomechanical curves", Kleshnev 2024, BioRow newsletter, n > 50,000 samples; independently 30-44%
    // per Soper & Hume 2004, Sports Med 34(12):825-848; normative targets per
    // Kleshnev 2011, "Biomechanics of Rowing", in Nolte (ed.) Rowing Faster 2nd
    // ed., Table 9.2 — sculling means 35.6% training / 38.2% racing, target
    // ≤40%). Converted to drive TIME via the handle-velocity shape map:
    // length fraction = normalized ∫handle speed; measured drive R = Vav/Vmax
    // ≈ 0.60 sculling (Kleshnev, "The Biomechanics of the Recovery Phase",
    // BioRow, n = 25,658; sculling 3-4% below the 60-70% band). A trapezoidal
    // profile at that R maps 40% of LENGTH → 48-49% of TIME (insensitive to R
    // across 0.60-0.65: peak shifts < 1% of length). Beta peak position
    // (α−1)/(α+β−2) = 0.48 of drive time → 38-39% of drive length: on the
    // sculling racing mean and inside the ≤40% target.
    // THE DISTRIBUTION FORM IS OURS, NOT KLESHNEV'S. Known misses, stated:
    // force at 60% of length 87% vs 74% target; at 90% of length ~34% vs ~0;
    // 13% catch force unreproducible (any beta with α>1 starts at 0 — the
    // catch force comes from oar inertia, un-modeled; do NOT spline around it).
    // A two-parameter beta cannot hit peak position AND the decay anchors
    // (needs α≈5.2/β≈5.5, far narrower than Kleshnev's 28-40% plateau stat);
    // peak position is this item's namesake defect and Kleshnev's most
    // precisely stated statistic, so it is the anchor matched exactly.
    // The TIME↔LENGTH map is PLACEHOLDER — SyncRow records oar angle vs time
    // directly and replaces it outright. After the blade-slip closure the measured statistics
    // above are VALIDATION CHECKS on the emergent force curve (blinded table,
    // the blinded output table), not inputs.
    alpha: 2.92,            // beta shape [DERIVED — fitted BY US as documented above]
    beta: 3.08,            // beta shape [DERIVED — fitted BY US as documented above]

    // --- Drag ---
    k_air: 0.28,            // N·s²/m², air drag on hull+crew [rebuilt at the drag rebuild]
    k_wave: 0.3,            // wave resistance Froude scaling factor

};

// (Blade-extraction note: the blade-extraction residual placeholder and its unsourced
// residual_fraction / residual_tau_fraction parameters were REMOVED with the blade-slip closure
// — extraction is handled by the closure: the prescribed force pulse ends
// smoothly (beta > 1) and the release angle is emergent. See FORCE: PROPULSIVE.)

// ========================== SIMULATION PARAMETERS ==========================

const SIM = {
    dt: 0.001            // s, integration timestep (1 ms)
};

// ---- Measurement window: exactly N complete cycles (decided 29 Aug 2026) ----
// The old fixed 40 s window contained a whole number of stroke cycles only by
// coincidence (exactly 20 at 30 SPM, but 13.33 at 20 SPM, 17.33 at 26, 21.33 at
// 32). Measured impact (26 Aug 2026 diagnostics): reported speed losses off by
// 1.7-2.8% relative at fractional-window rates vs ~0.3-0.5% at 30 SPM — a
// RATE-DEPENDENT error that systematically flattened the loss-vs-SPM trend and
// produced the non-monotonic recovered-k artifact. Fix: after the convergence
// criterion, average over exactly N complete cycles starting at a defined phase
// (a phi_1 = 0 boundary). v_mean, Δv, and all derived tables come from that
// window only.
const MEASURE = {
    nCycles: 20         // [ASSUMED — window length in complete cycles; the recorded fixed-measurement-window decision suggests N = 20.
                        //  Any integer N works; 20 cycles ≈ 40-60 s at essay rates, matching the
                        //  old window length while eliminating the fractional-cycle bias.]
};

// ---- Warmup convergence criterion (decided 29 Aug 2026) ----
// Warmup ends on a CONVERGENCE CRITERION, not a fixed clock: "the boat has
// settled" is a property of the boat, not of the clock. A fixed warmup encodes
// a settling-time assumption for one configuration — settling scales roughly as
// M_eff/(2·k·v), so it shifts with stroke rate, rower weight, and every The physics rebuild
// drag change. Measured context (30 SPM sync diagnostic): transient time
// constant ~4 s; at t=20 s the boat still sits 0.028 m/s (0.64%) below the
// converged limit cycle; the old 20-60 s window carried a -0.077% bias in mean
// speed.
// RULE: after each complete stroke, compare that stroke's MEAN boat speed to the
// previous stroke's mean (stroke-to-stroke means — the ~1 m/s within-stroke
// surge/slump is permanent, sits inside both averages, and cancels; only the
// start-from-rest transient shrinks stroke to stroke). When the difference falls
// below the threshold, the boat is settled and the measurement window opens.
const WARMUP = {
    // TIGHTENED at the 30 Aug 2026 code review: the original 0.01 m/s threshold
    // carried its own trigger condition — "if any future result depends on
    // differences below ~0.1%, tighten to 1e-3" — and The physics rebuild's emergent
    // headline IS ~0.08%, so the condition fired. Convergence study (review):
    // at 0.01 m/s the 20 ms loss read 0.0004% vs 0.0033% at 1e-3 (8× noise);
    // the 100 ms headline moved only 0.0815% → 0.0817%. Absolute sync speed
    // bias at 1e-3 is ~0.16 mm/s. Cost: longer warmup (~2× runtime).
    settleThreshold_mps: 0.001, // m/s [ASSUMED — tightened 30 Aug 2026, see above]
    maxStrokes: 400             // [ASSUMED — safety cap raised with the tighter threshold; exceeded => loud warning]
};

// ========================== OAR & BLADE ====================================

const DEG = Math.PI / 180;  // rad per degree [DERIVED — unit conversion]

const OAR = {
    length: 2.88,           // m, sculling oar overall [CITED: typical sculls 2.84-2.90 m — Concept2 scull specs; "Oar (sport rowing)" overview]
    inboard_nominal: 0.88,  // m [CITED: Kleshnev standard sculling rigging, Kleshnev 2007, BioRow newsletter No. 3]
    L_in: 0.84,             // m, ACTUAL inboard (pin to hand) [CITED: Kleshnev, "Amplitude and power of body segments", Kleshnev 2014, BioRow newsletter — the 110° sculling arc is quoted at actual inboard 0.84 m]
    blade_length: 0.46,     // m [CITED: Concept2 Inc., Smoothie2 Plain Edge scull blade specification — "industry standard at all levels"; see REFERENCES]
    blade_width: 0.215,     // m [CITED: same Concept2 specification]
    blade_fill: 0.89,       // planform area / bounding rectangle [ASSUMED — hatchet planform; swept ±20% in the sensitivity sweep]
    get L_out() {           // m, pin to blade CENTER [DERIVED: length − nominal inboard − blade_length/2 = 2.88 − 0.88 − 0.23 = 1.77. Replaces the uncited r_lever = 1.95, DELETED]
        return this.length - this.inboard_nominal - this.blade_length / 2;
    },
    get A_blade() {         // m² per blade [DERIVED: 0.46 × 0.215 × 0.89 ≈ 0.088]
        return this.blade_length * this.blade_width * this.blade_fill;
    },
    oars_per_rower: 2,      // sculling — two mirrored oars per rower [CITED — definitional].
    // Mirrored port/starboard oars: y-forces cancel; x-forces double. Since every
    // force relation is linear in blade area at fixed kinematics, the pair is
    // modeled as ONE equivalent oar with doubled blade area carrying the rower's
    // TOTAL handle force — exactly equivalent, stated here per the prescribed-force closure.

    // Blade hydrodynamics: Cl = Al·sin(2α), Cd = Ad·sin²(α), no stall at any
    // attack angle. [CITED: Caplan & Gardner 2007a, "A fluid dynamic
    // investigation of the Big Blade and Macon oar blade designs in rowing
    // propulsion", J Sports Sciences 25(6):643-650 — closed forms over the full
    // sweep range; amplitudes per Atkinson's (ROWING) distillation: Big Blade
    // Al = 1.25, Ad = 2.07 (flat plate 1.3/2.0; Macon 1.24/1.90).
    // KNOWN BIAS: Sliasas & Tullis 2009 CFD re-analysis finds the flume likely too
    // shallow, so coefficients may run ~6% high vs deep water — swept in item
    // the sensitivity sweep. Atkinson: flat-vs-curved coefficient difference moves mean shell
    // speed only ~2%, so precision risk lives in the closure, not here.]
    Al: 1.25,
    Ad: 2.07,
    theta_catch: -63 * DEG        // rad. CORRECTED 31 Aug 2026: the previous −55° was a SWEEP
                                  // catch angle applied to a sculling model. Elite SCULLING catch
                                  // angles are 61-66° (finish 42-44°, total 104-110°) [CITED:
                                  // BioRow measurement database as summarized in Rowing
                                  // Victoria's rigging guide (Wilson); BioRow/NK Oar Angle Guide
                                  // targets −65/+40 for Olympic-level scullers at 20-28 spm, noting
                                  // club rowers and higher rates run 3-5° shorter]. −63° is the
                                  // midpoint of the elite average range [midpoint choice ASSUMED;
                                  // the club/high-rate shortening is stated, not applied].
                                  // Release angle is EMERGENT; its literature reference value
                                  // (42-44° sculling) is kept in the accompanying fix specification
                                  // (simulator_fix_spec.md), not in this program (blind analysis).
};

// ========================== EFFORT =========================================
// One calibrated parameter in the whole model: peak TOTAL handle force per
// rower (both hands), shaped in drive time by the beta profile below.
// CALIBRATION PROTOCOL (the blade-slip closure): choose F_peak_ref ONCE so the synchronized
// 2x with two 80 kg rowers at cal_SPM reaches v_cal_target; regenerate every
// output from it. Rerun via --calibrate. Recalibrated after the drag rebuild,
// which changed total resistance, so the effort scale had to be re-fitted against it.

const EFFORT = {
    F_peak_ref: 352.0,      // N [DERIVED — calibration protocol above, rerun 31 Aug 2026 via --calibrate after the sculling catch-angle correction (−55°→−63°): sync 2x, 80 kg rowers, 30 SPM → 4.3000 m/s. Prior: 680 (estimate), 553.1 (pre-drag-rebuild), 341.8, 341.5 (sweep catch angle).]
    v_cal_target: 4.30,     // m/s [ASSUMED — club-level 2x pace anchor at 30 SPM; the single effort anchor, documented]
    cal_SPM: 30             // [ASSUMED — calibration stroke rate]
};

// Per-rower effort (needed for unequal rowers and the prescribed-force closure power renormalization).
// Scales linearly with rower mass [ASSUMED — carried over from the old F_0 ∝ mass rule].
let F_peak = [EFFORT.F_peak_ref, EFFORT.F_peak_ref];

// ========================== CLOSURE SOLVER PARAMS ==========================
const CLOSURE = {
    bracket_lo_frac: 0.7,   // warm-start bracket = cache × [lo_frac, hi_frac] [ASSUMED — numerical]
    bracket_hi_frac: 1.4,   // [ASSUMED — numerical]
    bracket_hi_init: 8,     // rad/s, cold-start upper bracket ≈ 4× typical peak oar angular velocity [ASSUMED — numerical]
    bracket_grow: 2,        // bracket expansion factor [ASSUMED — numerical]
    bracket_grow_max: 40,   // expansion cap [ASSUMED — numerical]
    bisect_resolution: 1e-7,// rad/s [ASSUMED — numerical]
    bisect_iters: 60,       // [ASSUMED — numerical]
    w2_eps: 1e-12           // m²/s², zero-flow guard [ASSUMED — numerical]
};
// ========================== MASS MODEL =====================================
// SINGLE SOURCE OF TRUTH — project architectural decision, 29 Aug 2026: ALL
// mass and ALL travel distances live in THIS block. Nothing outside it defines
// or redefines a mass, a mass fraction, or a travel distance. (The old model
// scattered them across the BOAT literal, configureBoat(), the slide functions
// and the momentum term — which is exactly how the 66%-fraction / 0.70 m-travel
// partial cancellation stayed invisible.)
//
// Segment mass fractions [CITED: de Leva, P. (1996), "Adjustments to
// Zatsiorsky-Seluyanov's segment inertia parameters", J. Biomechanics
// 29(9):1223-1230, Table 4, male — living young fit subjects]:
//   trunk 43.46% (= UPT 15.96 + MPT 16.33 + LPT/pelvis 11.17), thigh 14.16
//   each, shank 4.33 each, foot 1.37 each, head 6.94, upper arm 2.71 each,
//   forearm 1.62 each, hand 0.61 each.
// Translating fraction = trunk + both thighs = 43.46 + 28.32 = 71.78%
// [DERIVED]. Shanks and feet excluded as approximately offsetting (feet
// strapped to the stretcher; shanks pivot about a near-fixed ankle). Head is
// carried with the hull-frame remainder — an approximation absorbed in the
// same offset statement. Arms (9.88%) travel with the HANDLE, not the seat —
// attached to the emergent handle path (from the solved oar angle) below.
//
// Travel distances [project decision recorded in the accompanying fix specification (simulator_fix_spec.md); segment amplitudes from Kleshnev,
// "Amplitude and power of body segments", Kleshnev 2014, BioRow newsletter, n = 5437 — legs
// 33%, trunk 31%, arms 36% of the 1.61 m sculling stroke arc. Using segment
// amplitudes as centre-of-mass travel distances is a stated, recorded approximation]:
//   seat/pelvis  0.52 m total  [DERIVED — route 3: 0.33 × 1.61 = 0.53 ≈ 0.52.
//                               The 0.70 m figure is the TRACK length, an
//                               upper bound, NOT the travel]
//   trunk (rest) 0.50 m total  [DERIVED — 0.31 × 1.61]
//   thighs       0.26 m total  [ASSUMED — pinned at the knee, CM ≈ half the
//                               seat travel — recorded approximation]
//   arms         handle path   (emergent from the blade-slip closure's oar angle; x_handle =
//                               L_in·sinθ, ~1.17 m projected travel)
// PLACEHOLDER — every travel above is replaced by SyncRow/shell measurement
// when available.
//
// THE INTERACTION (why fraction and travel had to change together): the old 66% fraction (too low) and 0.70 m
// travel (~35% too long) were PARTLY CANCELLING; they are corrected TOGETHER
// here, with each mass given its own travel, so the accidental balance is
// replaced by the composed physical values.

const MASS = {
    shell: 27,             // kg [CITED: World Rowing Rules of Racing, Rule 32 Bye-Law — 2x minimum 27 kg; see REFERENCES]
    seat_hw: 2.5,          // kg seat hardware, travels with the seat [ASSUMED — typical seat+wheels; swept via the sensitivity sweep mass terms]
    rower_1: 80,           // kg [default; --weight]
    rower_2: 80,           // kg [default; --weight2 for unequal rowers]
    frac_pelvis: 0.1117,   // [CITED: de Leva Table 4, LPT]
    frac_trunk_rest: 0.3229, // [DERIVED: trunk 43.46% − pelvis 11.17%]
    frac_thighs: 0.2832,   // [CITED: de Leva, 2 × 14.16%]
    frac_arms: 0.0988,     // [CITED: de Leva, 2 × (2.71 + 1.62 + 0.61)%]
    travel_seat: 0.52,     // m total [DERIVED — see block header]  PLACEHOLDER
    travel_trunk: 0.50,    // m total [DERIVED — see block header]  PLACEHOLDER
    travel_thighs: 0.26,   // m total [ASSUMED — see block header]  PLACEHOLDER
    k_added: 0.07,         // entrained-water added mass [ASSUMED — no published value
                           // exists for racing shells. Anchors: Atkinson treats added
                           // mass as small & lumps it; BioRow (Kleshnev) quotes total
                           // inertial sensitivity ~0.15%/kg of displacement while pure
                           // mass-scaling of drag gives ~0.24%/kg — the gap bounds the
                           // added-mass share. Swept 0.03/0.07/0.10 in the sensitivity sweep]
    get frac_translating() { return this.frac_pelvis + this.frac_trunk_rest + this.frac_thighs; }, // 71.78% [DERIVED]
    get crew() { return this.rower_1 + this.rower_2; },
    get total() { return this.shell + this.crew + 2 * this.seat_hw; }
};

// ========================== ENERGY AUDIT PARAMS ============================
const AUDIT = {
    // Gate TIGHTENED 100× (30 Aug 2026 code review): after the exact arm-work
    // accounting and the symmetric window opening, measured worst residuals are
    // ~1e-4% of W_drag across 20-32 SPM, sync and offset (was 0.07% generic,
    // 3.1% at the 32 SPM quadrature breach). Gate set 10× above the observed
    // worst so any future dead term or quadrature break fires immediately —
    // the original 0.1% gate would now hide two orders of regression.
    tol_frac: 0.00001,   // fail threshold: |residual| > 0.001% of the cycle's W_drag [DERIVED — see above; the recorded energy-audit design originally 0.1%]
    wdrag_floor_J: 1.0   // J — denominator floor for near-zero-drag cycles (startup) [ASSUMED — numerical]
};

// ========================== STATE ==========================================
// NOTE on state.x: the 30 Aug 2026 code audit flagged x as integrated-but-never-read
// and slated it for removal. That premise no longer holds: the warmup criterion and the fixed measurement window made x
// LOAD-BEARING — cycle mean speeds are computed exactly as boundary-to-boundary
// distance over cycle time (see makeCycleAccumulator), which is what removes the
// sample-quantization error a sample-based mean would carry. x stays.
// Decision confirmed 30 Aug 2026.

// theta_1/theta_2 (the blade-slip closure): oar angle per rower. Integrated from the
// closure-solved angular velocity during the drive; prescribed analytically
// during the recovery (blade out of water) and hard-reset to theta_catch at
// each cycle boundary in RK4_step.
let state = { v: 0, x: 0, phi_1: 0, phi_2: 0, theta_1: OAR.theta_catch, theta_2: OAR.theta_catch,
              z: 0, w_z: 0, tau: 0, q_tau: 0, t: 0 };   // z/tau: the hull-attitude modes (heave m, pitch rad; bow-down +)
// Per-rower emergent release angle (captured when the force pulse ends) and
// closure warm-start cache. Diagnostics, not dynamics inputs.
let thetaRelease = [OAR.theta_catch, OAR.theta_catch];
let thetaDotCache = [1, 1];   // rad/s warm start [ASSUMED — numerical seed]
// Arm/handle finite-difference state (the single mass block): previous handle x and velocity
// per rower, and the frozen-per-step handle acceleration the momentum term uses.
let armFD = [{ x: null, vx: 0 }, { x: null, vx: 0 }];
let armAccel = [0, 0];
let omega = 0;          // rad/s, angular stroke rate = 2π × SPM/60
let delta = 0;          // rad, phase offset = ω × Δt

// Drive fraction parameters (set per simulation based on SPM)
let driveFrac = 0.50;
let drivePhaseEnd = Math.PI;

// ========================== DRIVE FRACTION =================================
// drive_time(SPM) = 1.45 - 0.015 × SPM  [seconds]
// [CITED: Kleshnev, BioRow, "Effect of stroke rate on rowing technique" —
//  reproduces the published regression: drive time shortens 15 ms per 1 SPM,
//  singles at 0.91 s at 36 SPM. The fitted line is for the 1x; applied here to
//  a 2x — same sculling family, expected error ~1-2%.]
// Cross-checks: drive(36) = 0.91 s matches Kleshnev; drive(30) = 1.0 s matches
// Rowperfect. Clamp bounds 0.30-0.60 on the drive fraction are [ASSUMED —
// sanity limits only, never active at legal SPM 16-44].

function getDriveFraction(SPM) {
    let drive_time = 1.45 - 0.015 * SPM;
    let cycle_time = 60 / SPM;
    return Math.max(0.30, Math.min(0.60, drive_time / cycle_time));
}

// ========================== UTILITY ========================================

function safeMod(a, b) {
    return ((a % b) + b) % b;
}

// ========================== BODY KINEMATICS ================================
// QUINTIC smoothstep S₅(t) = 6t⁵ − 15t⁴ + 10t³ for all prescribed body
// motions: S₅'' = 60t(2t−1)(t−1) is ZERO at both segment ends, so the body's
// acceleration is continuous (C²) across the drive/recovery boundaries. The
// old cubic had S'' = ±6 at the ends — an instantaneous acceleration step
// (~55 N per rower, twice per stroke) that no body can produce, kinking the
// hull-velocity trace and locally degrading RK4's smooth-RHS assumption.
// This is a claim about the REVERSAL taking finite time, not about seat speed:
// rowers reverse sharply (Kleshnev: the best crews check the boat hard but
// briefly at the catch) — sharp stays; instantaneous goes.
//
// CATCH FACTOR (Kleshnev CF): the SEAT reverses BEFORE the HANDLE at the
// catch — the two reversals are offset, which the old model collapsed to one
// simultaneous, zero-time event. Implemented as a uniform phase LEAD of the
// body motion relative to the force/oar phase: phi_body = phi − ω·CF_s (CF is
// negative, so the body runs ahead). Uniform lead shifts the finish reversal
// by the same amount — stated simplification; the published number is
// specifically the catch-side offset.
//   Default CF = −15 ms   [CITED: Kleshnev's sculling target; measured crews
//   −16 and −30 ms (Kleshnev 2025, BioRow newsletter No. 4, two M2- crews); −40 ms flagged
//   excessive; −50 ms and beyond = "slide shooting" with significant leg-power
//   loss. Physical scale (BioRow, "Visualisation of Catch Factor and Blade
//   Slip"): at optimal CF the seat travels ≤1.5 cm before the handle reverses
//   (0.5-1.0 cm in the measured scullers); 2-4 cm at −50 ms; 4-8 cm at −80.]
//   Per-rower parameter — SyncRow's oar data measures the handle side
//   directly, so measured CF can replace the default per rower once the seat
//   side is estimated from the Kleshnev curves (the rower-mass provenance work).
// PLACEHOLDER — to be populated/validated with SyncRow measured data (oar-only
// sensing cannot measure the seat side).

const BODY = {
    CF_s: [-0.015, -0.015]   // s, Catch Factor per rower [CITED — see block above]
};

function quinticSmoothstep(t) {
    return t * t * t * (10 + t * (-15 + 6 * t));            // 6t⁵−15t⁴+10t³
}

function quinticSmoothstepD2(t) {
    return 60 * t * (2 * t - 1) * (t - 1);                   // d²S₅/dt²
}

function bodyPhase(rower, phi) {
    return phi - omega * BODY.CF_s[rower];   // CF < 0 → body leads
}

function slidePosition(phi) {
    let phi_norm = safeMod(phi, 2 * Math.PI);
    if (phi_norm < drivePhaseEnd) {
        let t = phi_norm / drivePhaseEnd;
        return -1 + 2 * quinticSmoothstep(t);               // catch −1 → finish +1
    } else {
        let recPhase = 2 * Math.PI - drivePhaseEnd;
        let t = (phi_norm - drivePhaseEnd) / recPhase;
        return 1 - 2 * quinticSmoothstep(t);                // finish +1 → catch −1
    }
}

// Normalized-position second derivative w.r.t. real time (chain rule through
// the phase map), used by F_momentum with per-mass travels from MASS.
function slideAcceleration(phi) {
    let phi_norm = safeMod(phi, 2 * Math.PI);
    if (phi_norm < drivePhaseEnd) {
        let t = phi_norm / drivePhaseEnd;
        return 2 * quinticSmoothstepD2(t) * Math.pow(omega / drivePhaseEnd, 2);
    } else {
        let recPhase = 2 * Math.PI - drivePhaseEnd;
        let t = (phi_norm - drivePhaseEnd) / recPhase;
        return -2 * quinticSmoothstepD2(t) * Math.pow(omega / recPhase, 2);
    }
}

// ========================== FORCE: PROPULSIVE ==============================
// CLOSURE (the prescribed-force closure — Atkinson's method): the handle-force profile is PRESCRIBED
// (beta shape in drive time, one calibrated peak); oar kinematics are SOLVED.
// Massless-oar torque balance about the pin sets the blade-side force
// perpendicular to the shaft:
//     F_perp_required = F_handle × L_in / L_out
// (handle force assumed perpendicular to the shaft [ASSUMED — small-angle hand
// path approximation]). The hydrodynamic resultant on the blade must supply
// exactly that perpendicular component, which fixes the oar angular velocity
// θ̇ at each instant — solved by bracketed bisection with a warm-started
// bracket (Atkinson's ROWING solves the same balance ~500 instants per drive).
// Slip, θ(t), swept arc, release angle, and blade efficiency EMERGE.
//
// Blade kinematics (2D horizontal, x = travel; pin moves with hull at v):
//     u_x = v − L_out·θ̇·cosθ        blade-center velocity rel. water
//     u_y =   − L_out·θ̇·sinθ
//     w = −u  (flow over the blade),  face normal n̂ = (cosθ, sinθ)
//     sinα = (w·n̂)/|w|
//     F = q·[ Cl(α)·l̂ + Cd(α)·ŵ ],  q = ½·ρ·A_eff·|w|²,  l̂ ⊥ ŵ toward n̂
// Propulsive force on the hull = x-component of F (handle and stretcher forces
// are internal to the hull+rower+oar system).
//
// EXTRACTION (how the blade leaves the water under the prescribed-force closure): the drive ends when
// the prescribed force pulse ends (t_drive = 1) — the beta profile goes
// smoothly to zero because beta > 1, so no residual bolt-on is needed and the
// old inert placeholder block is REMOVED here as the blade-extraction handling anticipated. The
// RELEASE ANGLE is now EMERGENT and prints in the blinded table as the check
// against Kleshnev's release-slip framework (Kleshnev 2018, BioRow newsletter; measured
// magnitudes Kleshnev 2025, BioRow newsletter No. 4). A per-rower release-angle override returns only if
// SyncRow data shows the pulse-end convention misplaces the release.
//
// BLADE EFFICIENCY (replaces a deleted ×0.99 "blade efficiency" fudge factor that carried a false citation): derived per stroke —
//     E_blade = W_propulsive / W_handle
// W_handle = ∫ F_handle × (L_in·θ̇) dt over the drive (handle velocity in the
// hull frame — the rower's frame; stretcher work is not counted here, see the
// system-efficiency line in the blinded table). W_propulsive = ∫ F_prop,x × v dt.
// Efficiency is diagnostic OUTPUT only — never a target, never an input.

function betaFunction(t, alpha, beta) {
    if (t <= 0 || t >= 1) return 0;
    return Math.pow(t, alpha - 1) * Math.pow(1 - t, beta - 1);
}

function betaNormalization(alpha, beta) {
    let t_peak = (alpha - 1) / (alpha + beta - 2);
    return betaFunction(t_peak, alpha, beta);
}

// Prescribed handle force for one rower (N, total both hands) at phase phi.
function handleForce(rower, phi) {
    let phi_norm = safeMod(phi, 2 * Math.PI);
    if (phi_norm >= drivePhaseEnd) return 0;
    let t = phi_norm / drivePhaseEnd;   // drive-TIME fraction
    let norm = betaNormalization(BOAT.alpha, BOAT.beta);
    return F_peak[rower] * betaFunction(t, BOAT.alpha, BOAT.beta) / norm;
}

// Hydrodynamic force on one rower's blade pair at oar angle theta, angular
// velocity thetaDot, hull speed v. Returns x-force on hull and the component
// perpendicular to the shaft (the one the handle torque must balance).
function bladeHydro(theta, thetaDot, v) {
    let A_eff = OAR.oars_per_rower * OAR.A_blade;   // mirrored-pair equivalence, see OAR block
    let ux = v - OAR.L_out * thetaDot * Math.cos(theta);
    let uy = -OAR.L_out * thetaDot * Math.sin(theta);
    let wx = -ux, wy = -uy;
    let W2 = wx * wx + wy * wy;
    if (W2 < CLOSURE.w2_eps) return { Fx: 0, Fperp: 0 };
    let W = Math.sqrt(W2);
    let nhx = Math.cos(theta), nhy = Math.sin(theta);
    let sinA_signed = (wx * nhx + wy * nhy) / W;
    let sinA = Math.min(1, Math.abs(sinA_signed));
    let Cl = OAR.Al * 2 * sinA * Math.sqrt(Math.max(0, 1 - sinA * sinA)); // Al·sin(2α)
    let Cd = OAR.Ad * sinA * sinA;                                        // Ad·sin²(α)
    let q = 0.5 * CONSTANTS.rho_water * A_eff * W2;
    let wxh = wx / W, wyh = wy / W;
    // Lift ⊥ flow, oriented so the plate's normal-force component points to the
    // DOWNSTREAM side of the face: sign(l̂·n̂) = sign(w·n̂). (A plate is pushed
    // along +n̂ when the flow impinges on the −n̂ face, i.e. when w·n̂ > 0, and
    // along −n̂ otherwise — e.g. a freshly planted blade with the boat running
    // is swept toward the finish by the water, it does not resist.)
    let lx = -wyh, ly = wxh;
    if ((lx * nhx + ly * nhy) * sinA_signed < 0) { lx = -lx; ly = -ly; }
    let Fx = q * (Cl * lx + Cd * wxh);
    let Fy = q * (Cl * ly + Cd * wyh);
    return { Fx: Fx, Fperp: Fx * nhx + Fy * nhy };
}

// Solve θ̇ ≥ 0 such that the blade's perpendicular-to-shaft hydrodynamic force
// equals F_perp_req. F_perp is monotone increasing in θ̇ (∝ θ̇² at large θ̇), so
// bracketed bisection is robust; the bracket warm-starts from the previous
// solution since θ, v change little between calls.
function solveThetaDot(rower, theta, v, F_perp_req) {
    // KNOWN MODELING ARTIFACT (30 Aug 2026 code review): at the instant the force
    // pulse opens, this closure jumps the oar from rest to the flow-matching
    // rate θ̇* = v·cosθ/L_out — the handle velocity is DISCONTINUOUS at every
    // catch (~0.7 m/s in one step), which the quintic body-motion / Catch Factor upgrade's own principle ("nothing in
    // a human body changes force instantaneously") forbids for bodies. The
    // missing physics is OAR INERTIA (also the source of Kleshnev's 13%-of-max
    // catch force, the force-profile refit). Verified consequences are bounded: the arm-mass
    // reaction is a dt-invariant impulse (dt 1.0→0.5 ms leaves the headline
    // loss unchanged to 4 decimals) and the per-cycle energy audit integrates it
    // exactly. Listed in the header limitations; fix would be an oar-inertia
    // ODE, not a smoothing hack.
    if (bladeHydro(theta, 0, v).Fperp >= F_perp_req) {
        // θ̇ ≥ 0 clamp. With the corrected force orientation a static blade in
        // moving water is swept TOWARD the finish (Fperp(0) < 0 whenever v > 0),
        // so this branch cannot fire in normal running — it survives only as a
        // guard against a negative-rate solution (blade being backed), which the
        // model clamps to a held oar rather than integrating backing.
        return 0;
    }
    let lo = thetaDotCache[rower] * CLOSURE.bracket_lo_frac;
    let hi = thetaDotCache[rower] * CLOSURE.bracket_hi_frac;
    if (!(hi > 0)) { lo = 0; hi = CLOSURE.bracket_hi_init; }
    let g = 0;
    while (bladeHydro(theta, lo, v).Fperp > F_perp_req && g++ < CLOSURE.bracket_grow_max) lo /= CLOSURE.bracket_grow;
    g = 0;
    while (bladeHydro(theta, hi, v).Fperp < F_perp_req && g++ < CLOSURE.bracket_grow_max) hi *= CLOSURE.bracket_grow;
    for (let i = 0; i < CLOSURE.bisect_iters; i++) {
        let mid = 0.5 * (lo + hi);
        if (bladeHydro(theta, mid, v).Fperp < F_perp_req) lo = mid; else hi = mid;
        if (hi - lo < CLOSURE.bisect_resolution) break;
    }
    let td = 0.5 * (lo + hi);
    thetaDotCache[rower] = td;
    return td;
}

// Prescribed recovery oar path: blade out of the water, oar returns from the
// captured release angle to the catch angle over the recovery phase.
// PLACEHOLDER — to be populated/validated with SyncRow measured data.
function recoveryTheta(rower, phi_norm) {
    let recPhase = 2 * Math.PI - drivePhaseEnd;
    let u = (phi_norm - drivePhaseEnd) / recPhase;
    let s = quinticSmoothstep(u);   // C² return path (the quintic body-motion / Catch Factor upgrade)
    return thetaRelease[rower] + (OAR.theta_catch - thetaRelease[rower]) * s;
}

// One rower's propulsion state at (phi, theta, v): x-force on hull, solved
// angular velocity, and the prescribed handle force (for work integrals).
function rowerBlade(rower, phi, theta, v) {
    let phi_norm = safeMod(phi, 2 * Math.PI);
    if (phi_norm >= drivePhaseEnd) return { Fx: 0, thetaDot: 0, Fh: 0 };  // recovery: blade out
    let Fh = handleForce(rower, phi);
    let F_perp_req = Fh * OAR.L_in / OAR.L_out;   // massless-oar torque balance
    let td = solveThetaDot(rower, theta, v, F_perp_req);
    let h = bladeHydro(theta, td, v);
    return { Fx: h.Fx, thetaDot: td, Fh: Fh };
}
// ========================== FORCE: PROPULSIVE TOTAL ========================
// (The interference-term deletion: k_interference and bladeInterferenceFactor DELETED — the origin is
// de Brouwer, de Poel & Hofmijster 2013's own untested Discussion-section speculation about an
// EIGHT in ANTIPHASE, with no measurement, no coefficient, no functional form;
// the magnitude 0.40 and the sin(δ/2) shape were invented, and the geometry
// fails at our offsets — at 100 ms the second blade enters ~0.45 m AHEAD of
// the first blade's puddle. After the blade-slip closure, all inter-rower coupling is
// EMERGENT from blade hydrodynamics reading the shared hull velocity. If a
// phenomenological wake coupling is ever wanted for sensitivity exploration,
// it must return as an explicitly labelled free parameter citing de Poel as
// the ORIGIN OF AN UNTESTED HYPOTHESIS, never as calibration.)

function F_propulsive(s) {
    let b1 = rowerBlade(0, s.phi_1, s.theta_1, s.v);
    let b2 = rowerBlade(1, s.phi_2, s.theta_2, s.v);
    return { total: b1.Fx + b2.Fx, b1: b1, b2: b2 };
}
// ========================== HULL 2D GEOMETRY ===============================
// Shared machinery for the static-trim (--trim) and heave/pitch (
// --pitch) modes. Both follow Formaggia et al. 2009, "A model for the
// dynamics of rowing boats", Int. J. Numer. Meth. Fluids [CITED] — symmetry-
// plane motion, water pressure on the wetted surface linked to HYDROSTATICS
// (the paper's own first approximation for sculls), friction on the ITTC-1957
// line (their choice AND ours — consistent).
//
// Hull shape model [ASSUMED, stated]: parabolic waterplane b(x) = B·(1−(2x/L)²),
// wall-sided (rectangular) sections, uniform baseline draft d0 = V/A_wp.
// This shape is used ONLY to compute RATIOS and hydrostatic properties:
//   • S(z,τ)/S(0,0) multiplies the drag rebuild's cited wetted area, so baseline drag
//     is EXACTLY the cited-geometry value and only the fractional variation comes
//     from this model.
//   • A_wp (waterplane area), I_wp (waterplane 2nd moment) feed the
//     hydrostatic restoring terms — standard naval hydrostatics, which
//     Formaggia's static limit sanctions.
// KNOWN PROPERTY (stated before running): on a wall-sided hull with a
// symmetric waterplane, wetted area is FIRST-ORDER INSENSITIVE to trim
// (∫ x·b'(d) dx vanishes); only second-order end-emergence effects remain.
// If --trim reports ~nothing, that is the geometry speaking, not a bug.

const HULL2D = {
    nStations: 41,          // odd, Simpson-friendly [ASSUMED — numerical]
    x: [], b: [], w: [],    // stations, local beam, trapezoid weights
    A_wp: 0, I_wp: 0, d0: 0, S0_geom: 0,
    seat_x: [-1.25, 1.25],  // m, rower fore-aft seat positions about the
                            // waterplane centroid: rower 1 = stroke (stern,
                            // −x), rower 2 = bow (+x) — the later-catching
                            // rower for δ>0 sits in the bow [ASSUMED symmetric
                            // placement; ~2.5 m spacing per the recorded heave/pitch extension decision.
                            // WARNING (review): if the lagging seat is really
                            // the STROKE, the sign of every pitch-Δt coupling
                            // in the heave/pitch mode flips. Pin which seat lags with
                            // SyncRow before interpreting pitch results.]
    h_cm: 0.35,             // m, seated body-CM height above waterline for the
                            // dynamic pitch couple [ASSUMED — swept]
    zeta: 0.3,              // heave/pitch damping, fraction of critical
                            // [ASSUMED — Formaggia's damping came from RANSE
                            // fits not liftable from the paper; this value is
                            // a placeholder against resonance artifacts and is
                            // swept via --zeta]
    init() {
        this.x = []; this.b = []; this.w = [];
        let L = BOAT.L, N = this.nStations;
        for (let i = 0; i < N; i++) {
            let x = -L / 2 + L * i / (N - 1);
            this.x.push(x);
            this.b.push(BOAT.B * (1 - Math.pow(2 * x / L, 2)));
            this.w.push((i === 0 || i === N - 1) ? L / (N - 1) / 2 : L / (N - 1));
        }
        this.A_wp = 0; this.I_wp = 0;
        for (let i = 0; i < N; i++) {
            this.A_wp += this.w[i] * this.b[i];
            this.I_wp += this.w[i] * this.b[i] * this.x[i] * this.x[i];
        }
        this.d0 = (MASS.total / CONSTANTS.rho_water) / this.A_wp;
        this.S0_geom = this.S_geom(0, 0);
    },
    // Prism-model wetted area at heave z (down +) and trim tau (bow-down +)
    S_geom(z, tau) {
        let S = 0;
        for (let i = 0; i < this.nStations; i++) {
            let d = this.d0 + z + tau * this.x[i];
            if (d > 0) S += this.w[i] * (this.b[i] + 2 * d);
        }
        return S;
    },
    // Instantaneous crew pitch moment about the waterplane centroid (N·m):
    // static weight-shift part (g·Σ m·x) — what drives --trim — plus, in --pitch
    // mode, the dynamic couple from horizontal body accelerations at height
    // h_cm. Uses the same prescribed motions as F_momentum (single source).
    crewMoment(s, includeDynamic) {
        let M = 0;
        let phis = [s.phi_1, s.phi_2];
        let rowers = [MASS.rower_1, MASS.rower_2];
        for (let r = 0; r < 2; r++) {
            let pos = slidePosition(bodyPhase(r, phis[r]));   // [-1,1] shape (CF-shifted, the quintic body-motion / Catch Factor upgrade)
            let posAcc = slideAcceleration(bodyPhase(r, phis[r]));
            let mArm = MASS.frac_arms * rowers[r];
            let mSeat  = MASS.frac_pelvis * rowers[r] + MASS.seat_hw;
            let mTrunk = MASS.frac_trunk_rest * rowers[r];
            let mThigh = MASS.frac_thighs * rowers[r];
            let sSeat  = (MASS.travel_seat  / 2) * pos;
            let sTrunk = (MASS.travel_trunk / 2) * pos;
            let sThigh = (MASS.travel_thighs / 2) * pos;
            let xh = OAR.L_in * Math.sin(r === 0 ? s.theta_1 : s.theta_2);
            // static: weight × fore-aft position (seat base + displacement)
            M += CONSTANTS.g * (
                mSeat  * (this.seat_x[r] + sSeat) +
                mTrunk * (this.seat_x[r] + sTrunk) +
                mThigh * (this.seat_x[r] + sThigh) +
                mArm   * (this.seat_x[r] + xh) +
                (1 - MASS.frac_translating - MASS.frac_arms) * rowers[r] * this.seat_x[r]
            );
            if (includeDynamic) {
                // horizontal accel at height h above waterline → bow-down couple
                // when masses accelerate sternward [sign: F_react = −m·s̈ on hull
                // at height ~0; couple = +m·s̈·h_cm about y]
                let acc = (mSeat * (MASS.travel_seat / 2) +
                           mTrunk * (MASS.travel_trunk / 2) +
                           mThigh * (MASS.travel_thighs / 2)) * posAcc
                          + mArm * armAccel[r];
                M += acc * this.h_cm;
            }
        }
        // hull + oars assumed centered on the waterplane centroid [ASSUMED]
        return M;
    }
};

let trimMode = false;    // --trim  (the static-trim mode: static/quasi-static trim → S ratio)
let pitchMode = false;   // --pitch (23: heave+pitch dynamics → S ratio + diagnostics)

// ========================== FORCE: DRAG ====================================
// k_hull IS DELETED — NOT A PARAMETER (project decision, 29 Aug 2026). Drag is COMPUTED
// from measurable, cited quantities:
//
//   F_drag = ½·ρ·Cf(Re)·S · v|v| · (1 + k_wave·Fr²)  +  k_air·v|v|
//
//   S  = 2.30·√(V·Lwl), computed from displacement EACH RUN (V = M_total/ρ) —
//        so --weight scales drag by derivation, retiring the old hand-coded
//        displacement^(2/3) hack. At 192 kg default: S ≈ 3.16 m² (3.11-3.21
//        across the 9.8-10.4 m length question — robust to it).
//   Cf(Re) = 0.075 / (log₁₀Re − 2)²   [ITTC-1957 ship correlation line — the
//        same line Formaggia et al. use for rowing hulls; rowing-precedented],
//        Re = v·Lwl/ν. Speed-dependent: the hull surges/slows ~25% within each
//        stroke, so Cf drifts through the cycle. --constcf swaps in the
//        constant c = 0.00225 [Sliasas & Tullis 2009] for comparison.
//   Wave + air: the existing Froude term and k_air are RETAINED as the
//        ~15-20% non-friction remainder [skin friction is ~80-85% of racing-
//        shell total drag — Labbé et al. 2019, citing Lazauskas/Tuck; standard
//        in the builder literature]. Multiplicative wave-term structure
//        carried over [ASSUMED structural form].
//
// CHECKS, not inputs (the published constants the deleted k_hull collapsed):
//   this computation's equivalent at 4.3 m/s ≈ ½·997·Cf·S ≈ 3.5-3.6 N·s²/m²;
//   published anchors: 3.53 derived / 3.47 Sliasas-Tullis-scaled / 5.40 old
//   asserted / 5.97 van Holst-scaled — a 1.72× published disagreement. The
//   --dragscale sweep (equivalents 3.5 / 4.7 / 6.0) quantifies how much the
//   HEADLINE RATIO depends on it; mean speed certainly moves, the sync-async
//   ratio must be DEMONSTRATED robust, not assumed.
//
// Within-run variation: hull shape is fixed while rowing — only Cf(Re), the
// v|v| factor, and the Froude term vary with speed. Between-run variation:
// crew weight changes displacement → S recomputed. Wetted area also varies in
// reality as the hull pitches; that needs the vertical DOF the 1D model lacks
// (see LIMITATIONS in the header; the --trim mode is the cheap static approximation).

let constCfMode = false;   // --constcf
let dragScale = 1;         // --dragscale multiplier on the friction term (the drag rebuild sweep)

function wettedArea() {
    let V = MASS.total / CONSTANTS.rho_water;              // m³ displaced
    return BOAT.S_coeff * Math.sqrt(V * BOAT.L);           // m²
}

// Trim/pitch S-ratio (the hull-attitude modes): multiplies the drag rebuild's cited S. In --trim
// mode (z, τ) are the QUASI-STATIC hydrostatic-balance values; in --pitch mode
// they are the integrated heave/pitch states. Baseline (both off): ratio = 1.
let S_ratio_now = 1;

// Pure attitude→S-ratio map (the hull-attitude modes). In --trim mode the quasi-static
// hydrostatic balance is evaluated here (Formaggia's static limit:
// τ = M_crew/(ρ·g·I_wp); symmetric waterplane → z = 0) WITHOUT writing into
// the state. Every consumer of F_drag must refresh S_ratio_now via this
// function for its own state (derivatives and the work-integral sampler do).
function attitudeRatio(st) {
    if (trimMode) {
        let tau = HULL2D.crewMoment(st, false) / (CONSTANTS.rho_water * CONSTANTS.g * HULL2D.I_wp);
        return HULL2D.S_geom(0, tau) / HULL2D.S0_geom;
    }
    if (pitchMode) return HULL2D.S_geom(st.z, st.tau) / HULL2D.S0_geom;
    return 1;
}

function F_drag(v) {
    let av = Math.abs(v);
    let Fr = av / Math.sqrt(CONSTANTS.g * BOAT.L);
    let Cf;
    if (constCfMode || av < DRAG_NUM.v_eps) {
        Cf = BOAT.Cf_const;   // constant-c mode, and the v→0 guard (Re→0 blows up the ITTC line; at v≈0 drag ≈ 0 anyway)
    } else {
        let Re = av * BOAT.L / CONSTANTS.nu_water;
        let d = Math.log10(Re) - 2;
        Cf = 0.075 / (d * d);
    }
    let F_friction = 0.5 * CONSTANTS.rho_water * Cf * wettedArea() * S_ratio_now * v * av * dragScale;
    return F_friction * (1 + BOAT.k_wave * Fr * Fr) + BOAT.k_air * v * av;
}

const DRAG_NUM = {
    v_eps: 0.01   // m/s — below this, use Cf_const (ITTC line singular at Re→0) [ASSUMED — numerical guard]
};

// ========================== FORCE: INTERNAL MOMENTUM =======================
// Reaction on the hull from body masses accelerating relative to it:
//     F_mom = − Σ_rowers Σ_masses  m_j × s̈_j
// where s_j is each mass's displacement relative to the hull. Masses and
// travels come EXCLUSIVELY from the MASS block. Seat-borne (pelvis + seat
// hardware), trunk, and thighs follow the prescribed slide shape scaled to
// their own travels; ARMS follow the emergent handle path
// (x_handle = L_in·sinθ), whose acceleration is finite-differenced once per
// step in RK4_step (frozen within a step — explicit O(dt) treatment, dt = 1 ms;
// the per-cycle energy audit bounds the consequence).
// This is the "check" coaches feel — hull deceleration from crew motion alone.

// Arm reaction force from the frozen-per-step FD handle acceleration. This is
// PIECEWISE CONSTANT over each integration step by construction (armAccel
// updates post-step), so RK4 applies it as an exact impulse — and the per-cycle energy audit
// audit must integrate its work EXACTLY as F·Δx per step, not by trapezoid.
// (The 30 Aug 2026 code review: the trapezoid credited only HALF of each ~5 kN catch-instant
// spike's work — ~11 J/catch — which surfaced as a 1.7-3.1% audit breach at
// 32 SPM, where the catch phases coherently with the step grid. The DYNAMICS
// were always impulse-exact: dt-halving left the loss unchanged.)
function armForceNow() {
    return -(MASS.frac_arms * MASS.rower_1 * armAccel[0]
           + MASS.frac_arms * MASS.rower_2 * armAccel[1]);
}

function F_momentum_body(s) {
    let F = 0;
    let phis = [s.phi_1, s.phi_2];
    let rowers = [MASS.rower_1, MASS.rower_2];
    for (let r = 0; r < 2; r++) {
        let posAcc = slideAcceleration(bodyPhase(r, phis[r]));   // normalized shape accel, CF-led (the quintic body-motion / Catch Factor upgrade)
        let m_seat  = MASS.frac_pelvis * rowers[r] + MASS.seat_hw;
        let m_trunk = MASS.frac_trunk_rest * rowers[r];
        let m_thigh = MASS.frac_thighs * rowers[r];
        F -= (m_seat  * (MASS.travel_seat  / 2)
            + m_trunk * (MASS.travel_trunk / 2)
            + m_thigh * (MASS.travel_thighs / 2)) * posAcc;
    }
    return F;
}

function F_momentum(s) {
    return F_momentum_body(s) + armForceNow();
}
// ========================== EQUATION OF MOTION =============================
// TWO-MASS FORMULATION (van Holst-style multibody), derivation:
// Bodies: the HULL FRAME (shell + riggers + non-translating body mass, moving
// at v) and, per rower, point masses m_j at hull position + s_j(t), where each
// s_j is the PRESCRIBED motion relative to the hull (seat-borne, trunk, thigh
// masses on the slide shape; arm mass on the emergent handle path — all from
// the single mass block). Newton for the whole system, external x-forces only
// (handle, stretcher, pin and seat forces are internal):
//
//   d/dt [ M_total·v + Σ_j m_j·ṡ_j ] = Σ F_blade,x − F_drag(v) − F_air(v)
//
// Rearranged with the prescribed-motion terms on the right:
//
//   M_total·(1 + k_added)·dv/dt = F_propulsive − F_drag(v) − Σ_j m_j·s̈_j
//                                                            └─ = +F_momentum ─┘
//
// which is exactly what derivatives() integrates. Properties this buys:
//   • MOMENTUM-CONSERVING BY CONSTRUCTION — the internal-motion term is a
//     total derivative of Σ m_j ṡ_j, so over a closed cycle it transfers no
//     net impulse; no bolted-on reaction constant can drift the budget.
//   • DRAG ACTS ON THE HULL VELOCITY v, not the system centre of mass — the
//     whole point of the two-mass split. Internal mass motion reshapes the
//     hull-velocity waveform, and because drag is nonlinear (v|v|, plus
//     Cf(Re) after the drag rebuild), the PHASING of the two rowers' recoveries
//     changes mean dissipation. This is Boucher, Labbé & Clanet's (2017) proposed mechanism for the
//     async robot crews being slower, and — with blade slip — one of the two
//     legitimate emergent mean-speed channels in a 1D hull.
//   • Prior-art motivation: Caplan & Gardner's (2007b) own single-mass simulator
//     validated on MEAN velocity but reproduced the INSTANTANEOUS velocity
//     trace poorly — the documented failure this formulation fixes.
// Added mass k_added (entrained water reacting to hull acceleration) is
// applied to the total inertia [ASSUMED structural form — no published racing-
// shell value exists; swept in the sensitivity sweep].

function derivatives(s) {
    // The hull-attitude modes: refresh the wetted-area ratio for THIS state before the
    // surge equation (the 30 Aug 2026 code review: attitudeRatio is a pure function of the
    // state — no stage mutation, no stale global left for later readers).
    S_ratio_now = attitudeRatio(s);
    let prop = F_propulsive(s);
    return {
        // Effective inertia: total system mass + entrained water (the single mass block).
        // Structural form carried over [ASSUMED]; formalized in the two-mass formulation's
        // two-mass derivation below.
        dv: (prop.total - F_drag(s.v) + F_momentum(s)) / (MASS.total * (1 + MASS.k_added)),
        dx: s.v,
        dphi_1: omega,
        dphi_2: omega,
        // Drive: theta advances at the closure-solved rate. Recovery: theta is
        // prescribed analytically and assigned post-step in RK4_step, so its
        // derivative here is 0 (blade out of water, no dynamic role).
        dtheta_1: prop.b1.thetaDot,
        dtheta_2: prop.b2.thetaDot,
        // The heave/pitch mode heave/pitch dynamics (states inert unless --pitch):
        //   (M+A33)·ẇ = −ρg·A_wp·z − B33·w              (heave; forcing ≈ 0 —
        //     vertical body motion un-modeled, stated PLACEHOLDER)
        //   (Iyy+A55)·q̇ = M_crew(t) − ρg·I_wp·τ − B55·q  (pitch)
        // Restoring: standard hydrostatics (Formaggia's static-pressure
        // linkage). Added mass A33/A55: strip theory on the stated section
        // model. Damping B: [ASSUMED ζ of critical — see HULL2D.zeta].
        dz: pitchMode ? s.w_z : 0,
        dw_z: pitchMode ? pitch23.heaveAccel(s) : 0,
        dtau: pitchMode ? s.q_tau : 0,
        dq_tau: pitchMode ? pitch23.pitchAccel(s) : 0
    };
}

// The heave/pitch mode dynamic coefficients — derived once per simulate() from HULL2D and
// the MASS layout. See derivations at each term.
const pitch23 = {
    A33: 0, A55: 0, Iyy: 0, B33: 0, B55: 0, Wpitch: 0,
    init() {
        // Strip-theory heave added mass per unit length for a wall-sided
        // section of beam b: a33(x) = ρ·π·b(x)²/8 (flat-plate/half-cylinder
        // result) [DERIVED under the stated section model]
        let A33 = 0, A55 = 0;
        for (let i = 0; i < HULL2D.nStations; i++) {
            let a = CONSTANTS.rho_water * Math.PI * HULL2D.b[i] * HULL2D.b[i] / 8;
            A33 += HULL2D.w[i] * a;
            A55 += HULL2D.w[i] * a * HULL2D.x[i] * HULL2D.x[i];
        }
        this.A33 = A33; this.A55 = A55;
        // Pitch inertia about the waterplane centroid [DERIVED from the mass
        // layout: shell distributed over L (uniform-rod L²/12), crew lumped at
        // seat positions, seat hw with crew]
        this.Iyy = MASS.shell * BOAT.L * BOAT.L / 12
                 + (MASS.rower_1 + MASS.seat_hw) * HULL2D.seat_x[0] * HULL2D.seat_x[0]
                 + (MASS.rower_2 + MASS.seat_hw) * HULL2D.seat_x[1] * HULL2D.seat_x[1];
        // Damping at ζ of critical for each mode [ASSUMED — see HULL2D.zeta]
        let k_h = CONSTANTS.rho_water * CONSTANTS.g * HULL2D.A_wp;
        let k_p = CONSTANTS.rho_water * CONSTANTS.g * HULL2D.I_wp;
        this.B33 = 2 * HULL2D.zeta * Math.sqrt(k_h * (MASS.total + this.A33));
        this.B55 = 2 * HULL2D.zeta * Math.sqrt(k_p * (this.Iyy + this.A55));
        this.Wpitch = 0;
    },
    heaveAccel(s) {
        return (-CONSTANTS.rho_water * CONSTANTS.g * HULL2D.A_wp * s.z
                - this.B33 * s.w_z) / (MASS.total + this.A33);
    },
    pitchAccel(s) {
        return (HULL2D.crewMoment(s, true)
                - CONSTANTS.rho_water * CONSTANTS.g * HULL2D.I_wp * s.tau
                - this.B55 * s.q_tau) / (this.Iyy + this.A55);
    }
};

// ========================== RK4 INTEGRATOR =================================
// theta handling (the blade-slip closure): during the drive, theta integrates at the
// closure-solved rate carried in the stage derivatives. On the drive→recovery
// transition the emergent RELEASE ANGLE is captured; during recovery theta is
// assigned from the prescribed return path; on the recovery→drive wrap theta
// is hard-reset to theta_catch. Stage-level inconsistency at the transitions
// is bounded by one dt (1 ms) and is negligible.

function stageState(k, f) {
    return {
        v: state.v + f * k.dv,
        x: state.x + f * k.dx,
        phi_1: state.phi_1 + f * k.dphi_1,
        phi_2: state.phi_2 + f * k.dphi_2,
        theta_1: state.theta_1 + f * k.dtheta_1,
        theta_2: state.theta_2 + f * k.dtheta_2,
        z: state.z + f * k.dz,
        w_z: state.w_z + f * k.dw_z,
        tau: state.tau + f * k.dtau,
        q_tau: state.q_tau + f * k.dq_tau,
        t: state.t + f
    };
}

function RK4_step() {
    let prevPhi = [safeMod(state.phi_1, 2 * Math.PI), safeMod(state.phi_2, 2 * Math.PI)];

    let k1 = derivatives(state);
    let k2 = derivatives(stageState(k1, 0.5 * SIM.dt));
    let k3 = derivatives(stageState(k2, 0.5 * SIM.dt));
    let k4 = derivatives(stageState(k3, SIM.dt));

    state.v     += (SIM.dt / 6) * (k1.dv     + 2*k2.dv     + 2*k3.dv     + k4.dv);
    state.x     += (SIM.dt / 6) * (k1.dx     + 2*k2.dx     + 2*k3.dx     + k4.dx);
    state.phi_1 += (SIM.dt / 6) * (k1.dphi_1 + 2*k2.dphi_1 + 2*k3.dphi_1 + k4.dphi_1);
    state.phi_2 += (SIM.dt / 6) * (k1.dphi_2 + 2*k2.dphi_2 + 2*k3.dphi_2 + k4.dphi_2);
    state.theta_1 += (SIM.dt / 6) * (k1.dtheta_1 + 2*k2.dtheta_1 + 2*k3.dtheta_1 + k4.dtheta_1);
    state.theta_2 += (SIM.dt / 6) * (k1.dtheta_2 + 2*k2.dtheta_2 + 2*k3.dtheta_2 + k4.dtheta_2);
    if (pitchMode) {
        state.z     += (SIM.dt / 6) * (k1.dz     + 2*k2.dz     + 2*k3.dz     + k4.dz);
        state.w_z   += (SIM.dt / 6) * (k1.dw_z   + 2*k2.dw_z   + 2*k3.dw_z   + k4.dw_z);
        state.tau   += (SIM.dt / 6) * (k1.dtau   + 2*k2.dtau   + 2*k3.dtau   + k4.dtau);
        state.q_tau += (SIM.dt / 6) * (k1.dq_tau + 2*k2.dq_tau + 2*k3.dq_tau + k4.dq_tau);
        // Radiated/damped pitch+heave power (diagnostic only — NOT fed back
        // into surge; see the loose-coupling limitation in the header):
        pitch23.Wpitch += (pitch23.B55 * state.q_tau * state.q_tau
                         + pitch23.B33 * state.w_z * state.w_z) * SIM.dt;
    }
    state.t     += SIM.dt;

    let newPhi = [safeMod(state.phi_1, 2 * Math.PI), safeMod(state.phi_2, 2 * Math.PI)];
    let thetas = [state.theta_1, state.theta_2];
    for (let r = 0; r < 2; r++) {
        if (prevPhi[r] < drivePhaseEnd && newPhi[r] >= drivePhaseEnd) {
            thetaRelease[r] = thetas[r];          // emergent release angle (blinded diagnostic)
        }
        if (newPhi[r] >= drivePhaseEnd) {
            thetas[r] = recoveryTheta(r, newPhi[r]);   // prescribed recovery path
        }
        if (newPhi[r] < prevPhi[r]) {
            thetas[r] = OAR.theta_catch;          // new cycle: blade back at the catch
        }
    }
    state.theta_1 = thetas[0];
    state.theta_2 = thetas[1];

    // Arm mass on the handle path: update the arm/handle finite-difference chain (position →
    // velocity → acceleration), frozen until the next step.
    for (let r = 0; r < 2; r++) {
        let xh = OAR.L_in * Math.sin(thetas[r]);
        if (armFD[r].x === null) { armFD[r].x = xh; armFD[r].vx = 0; armAccel[r] = 0; continue; }
        let vx = (xh - armFD[r].x) / SIM.dt;
        armAccel[r] = (vx - armFD[r].vx) / SIM.dt;
        armFD[r].x = xh; armFD[r].vx = vx;
    }
}
// ========================== CYCLE ACCUMULATOR ==============================
// ONE phase-anchored cycle accumulator, shared by the warmup criterion (the warmup convergence criterion)
// and the measurement window (the fixed 20-cycle measurement window). Cycles are anchored at phi_1 = 0
// boundaries. Since dphi_1/dt = omega is constant and phi_1 starts at 0, the
// k-th boundary falls at exactly t = k · T_cycle — but that instant rarely
// lands on an integer multiple of dt, so boundary state (x, v) is interpolated
// sub-timestep between the surrounding integration steps (linear interpolation;
// error is O(dt²) per boundary, negligible at dt = 1 ms).
// A cycle's MEAN speed is computed from distance: (x_boundary_end −
// x_boundary_start) / T_cycle — exact for the interpolated boundaries, with no
// sample-quantization error.

function makeCycleAccumulator(T_cycle) {
    return {
        T_cycle: T_cycle,
        nextBoundaryT: T_cycle,
        prevBoundary: { t: 0, x: 0, v: 0 },
        lastMeanV: null,
        // Feed one integration step (state before, state after). Returns the
        // closed cycle's stats if a boundary was crossed, else null.
        step: function (before, after) {
            if (after.t < this.nextBoundaryT) return null;
            let frac = (this.nextBoundaryT - before.t) / SIM.dt;
            let xB = before.x + frac * (after.x - before.x);
            let vB = before.v + frac * (after.v - before.v);
            let meanV = (xB - this.prevBoundary.x) / this.T_cycle;
            let out = {
                meanV: meanV,
                prevMeanV: this.lastMeanV,
                frac: frac,   // sub-step fraction of the crossing (the per-cycle energy audit splits work integrals here)
                boundary: { t: this.nextBoundaryT, x: xB, v: vB }
            };
            this.lastMeanV = meanV;
            this.prevBoundary = out.boundary;
            this.nextBoundaryT += this.T_cycle;
            return out;
        }
    };
}

// ========================== METRICS ========================================

function computeMetrics(steady) {
    let v_mean = steady.reduce((a, b) => a + b) / steady.length;
    let variance = steady.reduce((sum, v) => sum + (v - v_mean) ** 2, 0) / steady.length;

    let v_max = -Infinity, v_min = Infinity;
    for (let v of steady) {
        if (v > v_max) v_max = v;
        if (v < v_min) v_min = v;
    }

    // Fluctuation power loss (the blinded output table): with F_drag ∝ v², dissipated power
    // ∝ v³, so at fixed mean speed a fluctuating hull dissipates ⟨v³⟩/⟨v⟩³ − 1
    // more than a steady one (this is the effect Greidanus et al. 2016 MEASURE —
    // their correct role, as an output check).
    let v3_mean = steady.reduce((sum, v) => sum + v * v * v, 0) / steady.length;
    let fluctuationLoss = v3_mean / (v_mean * v_mean * v_mean) - 1;

    return {
        v_mean,
        v_std: Math.sqrt(variance),
        variance,
        Delta_v: v_max - v_min,
        v_max,
        v_min,
        fluctuationLoss
    };
}

// ========================== MAIN SIMULATION ================================

function simulate(SPM, Delta_t) {
    omega = 2 * Math.PI * (SPM / 60);
    delta = omega * Delta_t;
    driveFrac = getDriveFraction(SPM);
    drivePhaseEnd = 2 * Math.PI * driveFrac;
    // SIGN CONVENTION: delta > 0 means rower 2 catches delta/omega seconds AFTER rower 1.
    // A larger phase means further along the cycle, so the LAGGING rower must start at a
    // SMALLER phase: phi_2 = -delta. (Previous code used phi_2 = +delta, which put rower 2
    // AHEAD while F_propulsive penalized rower 2 as the late one — invisible today because
    // both rowers are identical, wrong the moment they become distinguishable.)
    // Convention verified at startup by runStartupSelfCheck().
    state = { v: 0, x: 0, phi_1: 0, phi_2: -delta,
              theta_1: OAR.theta_catch, theta_2: OAR.theta_catch,
              z: 0, w_z: 0, tau: 0, q_tau: 0, t: 0 };
    if (trimMode || pitchMode) HULL2D.init();
    if (pitchMode) pitch23.init();
    S_ratio_now = 1;
    thetaRelease = [OAR.theta_catch, OAR.theta_catch];
    thetaDotCache = [1, 1];
    armFD = [{ x: null, vx: 0 }, { x: null, vx: 0 }];
    armAccel = [0, 0];

    let T_cycle = 2 * Math.PI / omega;
    let acc = makeCycleAccumulator(T_cycle);

    // ---- Warmup: run until stroke-to-stroke mean speed converges (the warmup convergence criterion) ----
    let settled = false;
    let strokes = 0;
    while (!settled && strokes < WARMUP.maxStrokes) {
        let before = { t: state.t, x: state.x, v: state.v };
        RK4_step();
        let cyc = acc.step(before, state);
        if (cyc !== null) {
            strokes++;
            if (cyc.prevMeanV !== null &&
                Math.abs(cyc.meanV - cyc.prevMeanV) < WARMUP.settleThreshold_mps) {
                settled = true;
            }
        }
    }
    if (!settled) {
        console.error(`WARNING: warmup did NOT converge within ${WARMUP.maxStrokes} strokes ` +
                      `(SPM=${SPM}, offset=${Delta_t}s, threshold=${WARMUP.settleThreshold_mps} m/s). ` +
                      `Proceeding with unconverged state — results suspect.`);
    }

    // ---- Measurement window: exactly MEASURE.nCycles complete cycles (the fixed 20-cycle measurement window) ----
    // The window opens at the phi_1 = 0 boundary where convergence was declared
    // (acc.prevBoundary, whose x/v were interpolated sub-timestep) and closes
    // exactly nCycles boundaries later.
    // v_mean comes from distance over the exact window: (x_end - x_start) / (N·T)
    // — no sample-quantization error.
    // v_max/v_min include every in-window sample plus both interpolated boundary
    // values. variance is computed over the in-window samples; its window edges
    // are quantized to the timestep (edge error ≤ dt / (N·T) ≈ 0.003%, accepted).
    //
    // WORK INTEGRALS (the blade-efficiency derivation — efficiency ships with the closure; the per-cycle energy audit
    // energy audit reuses them): accumulated once per step over this window
    // using left-endpoint sums (O(dt) quadrature — same order as the sampling
    // quantization already accepted above).
    //   W_handle,i = ∫ F_h,i · (L_in·θ̇_i) dt   (hull-frame handle work, drive)
    //   W_prop     = ∫ F_prop,x,total · v dt
    //   W_drag     = ∫ F_drag · v dt
    //   W_mom      = ∫ F_momentum · v dt
    // The 30 Aug 2026 code review (audit, first-cycle): do NOT open the window at the warmup-exit
    // boundary — that crossing's post-boundary tail (up to one full step,
    // depending on where float drift puts the boundary within the step) was
    // integrated inside the warmup loop and never accumulated, leaving the
    // first measured cycle short by ~0.5-1 J (0.07-0.17% — surfaced by the
    // the code review's exact arm accounting driving interior residuals to ~1e-4 J).
    // Instead the measurement loop runs to the NEXT boundary and opens the
    // window there, through exactly the same interpolated-split code path as
    // every interior boundary. Cost: ≤1 extra cycle per simulation.
    let windowOpen = false;
    let windowStart = null;
    let samples = [];
    let cyclesDone = 0;
    let windowEnd = null;

    // Work-rate integrands at a state (the per-cycle energy audit; trapezoidal accumulation —
    // O(dt²), so the audit gate below tests physics, not quadrature).
    let integr = (st) => {
        S_ratio_now = attitudeRatio(st);   // the 30 Aug 2026 code review: no stale attitude in the audit integrals
        let prop = F_propulsive(st);
        return {
            ph1: prop.b1.Fh * OAR.L_in * prop.b1.thetaDot,
            ph2: prop.b2.Fh * OAR.L_in * prop.b2.thetaDot,
            pp: prop.total * st.v,
            pd: F_drag(st.v) * st.v,
            pm: F_momentum_body(st) * st.v   // arm work integrated exactly per step below
        };
    };
    let zero = () => ({ ph1: 0, ph2: 0, pp: 0, pd: 0, pm: 0 });
    // Integrands linearly interpolated to the boundary at fraction f of the step
    let lerpI = (a, b, f) => ({ ph1: a.ph1 + f * (b.ph1 - a.ph1), ph2: a.ph2 + f * (b.ph2 - a.ph2),
                                pp: a.pp + f * (b.pp - a.pp), pd: a.pd + f * (b.pd - a.pd), pm: a.pm + f * (b.pm - a.pm) });
    let win = zero();            // window totals
    let cyc19 = zero();          // current-cycle totals (per-cycle energy audit)
    let cycStartV = 0;           // set when the window opens
    let M_lhs = MASS.total * (1 + MASS.k_added);
    let residualWorst = 0;       // worst |residual|/W_drag over the window (blinded table row)
    let residualFail = false;
    // The blinded output table: capture rower 1's handle force vs oar angle over the LAST
    // measured cycle — the emergent force-curve statistics come from this.
    let forceTrace = [];
    let prevI = integr(state);

    while (cyclesDone < MEASURE.nCycles) {
        let before = { t: state.t, x: state.x, v: state.v };
        let Farm = armForceNow();   // the value RK4 will apply over this step
        if (windowOpen && cyclesDone === MEASURE.nCycles - 1) {
            let pn = safeMod(state.phi_1, 2 * Math.PI);
            if (pn < drivePhaseEnd) forceTrace.push({ theta: state.theta_1, Fh: handleForce(0, state.phi_1) });
        }
        RK4_step();
        let currI = integr(state);
        let cyc = acc.step(before, state);
        let addTo = (acc_, a, b, w) => {
            acc_.ph1 += 0.5 * (a.ph1 + b.ph1) * w; acc_.ph2 += 0.5 * (a.ph2 + b.ph2) * w;
            acc_.pp  += 0.5 * (a.pp  + b.pp ) * w; acc_.pd  += 0.5 * (a.pd  + b.pd ) * w;
            acc_.pm  += 0.5 * (a.pm  + b.pm ) * w;
        };
        if (cyc === null) {
            if (windowOpen) {
                addTo(win, prevI, currI, SIM.dt);
                addTo(cyc19, prevI, currI, SIM.dt);
                let armW = Farm * (state.x - before.x);   // exact for piecewise-constant force
                win.pm += armW; cyc19.pm += armW;
                samples.push(state.v);
            }
        } else if (!windowOpen) {
            // ---- WINDOW OPEN at this boundary, via the standard split ----
            windowOpen = true;
            windowStart = cyc.boundary;
            cycStartV = cyc.boundary.v;
            samples.push(cycStartV);
            if (pitchMode) pitch23.Wpitch = 0;
            let f = cyc.frac;
            let bI = lerpI(prevI, currI, f);
            addTo(win, bI, currI, (1 - f) * SIM.dt);
            addTo(cyc19, bI, currI, (1 - f) * SIM.dt);
            let armW_open = Farm * (state.x - cyc.boundary.x);
            win.pm += armW_open; cyc19.pm += armW_open;
            samples.push(state.v);
        } else {
            // Split this step's contribution at the interpolated boundary.
            let f = cyc.frac;
            let bI = lerpI(prevI, currI, f);
            addTo(win, prevI, bI, f * SIM.dt);
            addTo(cyc19, prevI, bI, f * SIM.dt);
            let armW_pre = Farm * (cyc.boundary.x - before.x);   // exact, pre-boundary share
            win.pm += armW_pre; cyc19.pm += armW_pre;
            // ---- The per-cycle energy audit: per-cycle energy audit at the closed boundary ----
            // The EOM guarantees M_lhs·v̇·v = (F_prop − F_drag + F_mom)·v, so per
            // cycle: W_prop + W_mom − W_drag = ΔKE_hull-frame (≈ 0 at the limit
            // cycle). Any dead term, double count, or broken quadrature breaks
            // this identity — the check that would have caught the inert
            // residual and the double-counted Greidanus et al. 2016 fluctuation loss on its own.
            let dKE = 0.5 * M_lhs * (cyc.boundary.v * cyc.boundary.v - cycStartV * cycStartV);
            let residual = cyc19.pp + cyc19.pm - cyc19.pd - dKE;
            let rel = Math.abs(residual) / Math.max(cyc19.pd, AUDIT.wdrag_floor_J);
            if (rel > residualWorst) residualWorst = rel;
            if (rel > AUDIT.tol_frac) {
                residualFail = true;
                console.error(`ENERGY AUDIT FAILURE: cycle residual ${(rel*100).toFixed(4)}% of W_drag ` +
                              `(residual=${residual.toExponential(3)} J, W_drag=${cyc19.pd.toFixed(1)} J, ` +
                              `SPM=${SPM}, offset=${Delta_t}s, cycle ending t=${cyc.boundary.t.toFixed(2)}s)`);
            }
            cyc19 = zero();
            cycStartV = cyc.boundary.v;
            cyclesDone++;
            if (cyclesDone === MEASURE.nCycles) {
                // The 30 Aug 2026 code review: stop the window EXACTLY at the boundary — the old
                // code added the post-boundary sliver of this step to win.
                windowEnd = cyc.boundary;
                samples.push(windowEnd.v);
                prevI = currI;
                break;
            }
            addTo(win, bI, currI, (1 - f) * SIM.dt);
            addTo(cyc19, bI, currI, (1 - f) * SIM.dt);
            let armW_post = Farm * (state.x - cyc.boundary.x);
            win.pm += armW_post; cyc19.pm += armW_post;
            samples.push(state.v);
        }
        prevI = currI;
    }

    let metrics = computeMetrics(samples);
    // Exact distance-based mean over the integer-cycle window supersedes the
    // sample mean (which is edge-quantized):
    metrics.v_mean = (windowEnd.x - windowStart.x) / (MEASURE.nCycles * T_cycle);

    let windowT = MEASURE.nCycles * T_cycle;
    let Wh = [win.ph1, win.ph2];
    let energetics = {
        W_handle_1: Wh[0], W_handle_2: Wh[1],
        W_handle: Wh[0] + Wh[1],
        W_prop: win.pp, W_drag: win.pd, W_mom: win.pm,
        P_rower_1: Wh[0] / windowT,           // mean handle power per rower (the prescribed-force closure renorm)
        P_rower_2: Wh[1] / windowT,
        E_blade: win.pp / (Wh[0] + Wh[1]),    // blinded output — see the blinded output table
        thetaRelease_1: thetaRelease[0],      // emergent release angles (blinded)
        thetaRelease_2: thetaRelease[1],
        W_pitch_rad: pitchMode ? pitch23.Wpitch : 0,   // 23 diagnostic (measurement window)
        residualWorst: residualWorst,         // worst per-cycle |audit residual| / W_drag
        residualFail: residualFail,
        forceTrace: forceTrace,               // rower 1, last measured cycle (blinded stats input)
        window_T: windowT
    };

    return {
        metrics, energetics,
        SPM, Delta_t, delta, driveFrac
    };
}

// ========================== POWER RENORMALIZATION ==========================
// An offset rower slips differently, so sync-vs-async at fixed handle force
// conflates the timing penalty with a power difference. For offset runs, each
// rower's effort scale is iterated until their mean handle power matches the
// synchronized baseline's per-rower power. Equal-power results are the
// headline; --rawpower disables (raw = fixed handle-force profiles).

// Numerical noise floor statement printed under every loss table [DERIVED —
// convergence study 30 Aug 2026: at these WARMUP/RENORM defaults the 20 ms /
// 30 SPM loss moved 0.00330% → 0.00335% under a further 10× tightening of
// both tolerances, and the 100 ms headline was dt-invariant under dt 1.0 →
// 0.5 ms. Reproduce by rerunning that study after any solver change.]
const NOISE_FLOOR_NOTE =
    "  Numerical floor: losses converged to ~±0.0001 pp at current tolerances " +
    "(convergence study 30 Aug 2026); treat values below ~0.0005% as at-floor.";

const RENORM = {
    // TIGHTENED 30 Aug 2026 (the 30 Aug 2026 code review): the old tol 0.005 left per-rower
    // power residuals of ±0.17% — through the cube law that is ±0.05% speed each,
    // the same order as the 0.08% headline. The headline survived only because
    // the residuals are ANTI-SYMMETRIC (renorm over/undershoots the two rowers
    // oppositely) and cancel in net crew power. That cancellation is now
    // ASSERTED, not assumed: see the net-power gate in simulatePair.
    iters: 8,        // outer iteration cap [ASSUMED — numerical; breaks early on tol]
    tol_frac: 0.0005,// per-rower relative power tolerance [ASSUMED — numerical]
    exp_damp: 0.7,   // update exponent: scale *= (P_ref/P)^exp_damp [ASSUMED — numerical damping]
    net_gate: 0.0005 // |net crew power residual| gate — breach prints loudly
                     // [ASSUMED — sized so the implied speed bias (~1/3 of it)
                     // stays an order below the headline]
};

let equalPowerMode = true;   // set from CLI (--rawpower turns off)

// syncResult may be passed in to avoid re-simulating the baseline (review perf
// fix — the old version re-ran sync on every call, ~40% of the default run).
// Renorm targets each rower's OWN sync power (the 30 Aug 2026 code review: the old single
// P_rower_1 reference silently mis-normalized unequal --weight2 crews).
function simulatePair(SPM, Delta_t, syncResult, equalPower) {
    if (equalPower === undefined) equalPower = equalPowerMode;   // CLI default; callers may force raw
    let basePeak = F_peak.slice();
    let sync = syncResult !== undefined ? syncResult : simulate(SPM, 0);
    let off;
    if (Delta_t === 0) {
        off = sync;
    } else if (!equalPower) {
        off = simulate(SPM, Delta_t);
    } else {
        let P_ref = [sync.energetics.P_rower_1, sync.energetics.P_rower_2];
        for (let i = 0; i < RENORM.iters; i++) {
            off = simulate(SPM, Delta_t);
            let r1 = P_ref[0] / off.energetics.P_rower_1;
            let r2 = P_ref[1] / off.energetics.P_rower_2;
            if (Math.abs(r1 - 1) < RENORM.tol_frac && Math.abs(r2 - 1) < RENORM.tol_frac) break;
            F_peak[0] *= Math.pow(r1, RENORM.exp_damp);
            F_peak[1] *= Math.pow(r2, RENORM.exp_damp);
        }
        off = simulate(SPM, Delta_t);
        // The 30 Aug 2026 code review gate: the per-rower residuals are anti-symmetric and
        // cancel in net crew power — assert it instead of assuming it.
        let netResid = (off.energetics.P_rower_1 + off.energetics.P_rower_2)
                     / (P_ref[0] + P_ref[1]) - 1;
        if (Math.abs(netResid) > RENORM.net_gate) {
            console.error(`RENORM NET-POWER GATE BREACH: net crew power residual ${(netResid*100).toFixed(4)}% ` +
                          `exceeds ${(RENORM.net_gate*100).toFixed(2)}% (SPM=${SPM}, offset=${Delta_t}s). ` +
                          `The equal-power loss at this point carries a power bias of the same order — do not trust it.`);
        }
        off.energetics.net_power_residual = netResid;
    }
    F_peak = basePeak.slice();   // restore configured efforts
    return { sync: sync, off: off };
}
// ========================== BLINDED OUTPUT TABLE ===========================
// BLIND ANALYSIS (project decision, 29 Aug 2026): the program prints the QUANTITY and OUR
// NUMBER — never literature values, never expected ranges, never PASS/FAIL.
// Reference values live in the accompanying fix specification (simulator_fix_spec.md) (blade efficiency, speed swing, force-curve shape, drag anchors, with their sources) and
// are consulted ONCE, when a result is final. INPUTS stay cited and visible in
// the settings block; OUTPUTS below are what we would tune toward, so they are
// blinded. The only visible internal row is the energy residual (the per-cycle energy audit —
// internal consistency, not a literature value).

// Emergent force-curve statistics from the captured drive trace: drive-LENGTH
// fraction is swept oar angle (θ − θ_catch)/(θ_release − θ_catch); the handle
// force is plotted against it, which is exactly Kleshnev's axis.
function forceCurveStats(trace, thetaRel) {
    if (!trace || trace.length < 3) return null;
    let arc = thetaRel - OAR.theta_catch;
    if (!(arc > 0)) return null;
    let Fmax = 0, peakLen = 0;
    for (let p of trace) if (p.Fh > Fmax) { Fmax = p.Fh; peakLen = (p.theta - OAR.theta_catch) / arc; }
    if (!(Fmax > 0)) return null;
    let atFrac = (frac) => {
        let best = null, bestD = Infinity;
        for (let p of trace) {
            let lf = (p.theta - OAR.theta_catch) / arc;
            let d = Math.abs(lf - frac);
            if (d < bestD) { bestD = d; best = p; }
        }
        return best.Fh / Fmax;
    };
    return { peakLen: peakLen, catchFrac: atFrac(0), at60Frac: atFrac(0.60) };
}

function printBlindedTable(sync, off) {
    let e = sync.energetics, m = sync.metrics;
    // System efficiency: useful dissipation (hull drag work) over total rower
    // work = handle work + body-motion work (the crew's work moving its own
    // mass against hull acceleration, ∫Σmⱼṡⱼ·v̇ dt — accumulated as W_mom's
    // hull-side twin; here taken as W_mom since per cycle they balance).
    // Definition stated here once; a diagnostic OUTPUT, not a target.
    let sysEff = e.W_drag / (e.W_handle + Math.max(0, e.W_mom));
    let fc = forceCurveStats(e.forceTrace, e.thetaRelease_1);
    let rel1 = e.thetaRelease_1 / DEG, arc1 = (e.thetaRelease_1 - OAR.theta_catch) / DEG;

    console.log("\n  ============== BLINDED OUTPUT TABLE ==============");
    console.log("  Quantity and OUR number only — literature reference values are kept outside this program (simulator_fix_spec.md).");
    console.log(`    Blade efficiency               ${e.E_blade.toFixed(3)}`);
    console.log(`    System efficiency              ${sysEff.toFixed(3)}   (def: W_drag/(W_handle+W_hullmom); W_hullmom = ∫F_mom·v dt, hull-side — NOT crew metabolic work)`);
    if (fc) {
        console.log(`    Force peak position            ${(fc.peakLen*100).toFixed(1)}% of drive length`);
        console.log(`    Force at catch                 ${(fc.catchFrac*100).toFixed(1)}% of Fmax`);
        console.log(`    Force at 60% of drive length   ${(fc.at60Frac*100).toFixed(1)}% of Fmax`);
    }
    console.log(`    Release angle / arc (emergent) ${rel1.toFixed(1)}° / ${arc1.toFixed(1)}°`);
    console.log(`    Speed swing Δv/v̄               ${(m.Delta_v/m.v_mean*100).toFixed(1)}%`);
    console.log(`    Fluctuation power loss         ${(m.fluctuationLoss*100).toFixed(2)}%   (⟨v³⟩/⟨v⟩³−1)`);
    console.log(`    Mean speed (sync)              ${m.v_mean.toFixed(4)} m/s @ ${sync.SPM} SPM`);
    if (off) {
        let faster = off.metrics.v_mean < sync.metrics.v_mean ? "sync faster" :
                     off.metrics.v_mean > sync.metrics.v_mean ? "async faster" : "equal";
        console.log(`    Sync vs async (${(off.Delta_t*1000).toFixed(0)} ms)         ${faster}: ${sync.metrics.v_mean.toFixed(4)} vs ${off.metrics.v_mean.toFixed(4)} m/s`);
        console.log(`    Async speed swing Δv/v̄        ${(off.metrics.Delta_v/off.metrics.v_mean*100).toFixed(1)}%`);
    }
    console.log(`    Energy residual (worst cycle)  ${(e.residualWorst*100).toExponential(2)}% of W_drag ${e.residualFail ? "*** EXCEEDED GATE ***" : "(gate " + (AUDIT.tol_frac*100) + "%)"}`);
    if (pitchMode) console.log(`    Pitch+heave radiated work      ${(e.W_pitch_rad / e.W_drag * 100).toFixed(3)}% of W_drag   (diagnostic only — NOT fed back into surge)`);
    console.log("  ================================================================");
}

// ========================== SPEED LOSS TABLE GENERATOR =====================
// Offset rows use simulatePair: equal-power renormalization ON by default
// (the prescribed-force closure — the headline experiment holds per-rower average power fixed);
// --rawpower compares at fixed handle-force profiles instead.

function generateSpeedLossTable(SPM, offsets_ms) {
    let sync = simulate(SPM, 0);
    let table = [];

    for (let dt_ms of offsets_ms) {
        let result = dt_ms === 0 ? sync : simulatePair(SPM, dt_ms / 1000, sync).off;
        let loss = (sync.metrics.v_mean - result.metrics.v_mean) / sync.metrics.v_mean;
        // Raw-power companion (the 30 Aug 2026 code review found that equal-power vs raw differ
        // ~7× at 100 ms — a design-dependent number must ship with its twin):
        let rawLoss = 0;
        if (dt_ms !== 0) {
            let rawOff = simulatePair(SPM, dt_ms / 1000, sync, false).off;
            rawLoss = (sync.metrics.v_mean - rawOff.metrics.v_mean) / sync.metrics.v_mean;
        }
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
            dv_change_pct: dv_change * 100,
            raw_loss_pct: rawLoss * 100
        });
    }

    return { SPM, driveFrac: getDriveFraction(SPM), sync_v: sync.metrics.v_mean,
             sync_dv: sync.metrics.Delta_v, table };
}
// ========================== STARTUP SELF-CHECK =============================
// There is no test harness; unit assertions run here at startup and exit loudly.
// CHECK 1 (the offset sign-convention fix, reworked at the interference-term deletion — interference is gone): the sign
// convention itself is asserted through the live code path. With delta > 0 and
// phases initialized as simulate() does (phi_1 = 0, phi_2 = -delta), at a time
// after rower 1's catch but before rower 2's, rower 1 must be pulling and
// rower 2 must still be in recovery (zero handle force, blade out).
// CHECK 2 (the blade-slip closure): the closure honors the torque balance — the solved
// blade force perpendicular to the shaft equals F_handle × L_in / L_out.

const SELF_CHECK = {
    spm: 30,             // [ASSUMED — arbitrary representative rate; any legal rate works]
    offset_s: 0.100,     // [ASSUMED — arbitrary nonzero offset]
    catchGapFrac: 0.5,   // evaluate midway between the two catches [ASSUMED — any value in (0,1) works]
    drivePointFrac: 0.4, // [ASSUMED — arbitrary mid-drive evaluation point for CHECK 2]
    theta_test: -10 * DEG, // [ASSUMED — arbitrary mid-drive oar angle]
    v_test: 4.0,         // m/s [ASSUMED — arbitrary hull speed in the working range]
    F_peak_test: 600,    // N [ASSUMED — arbitrary nonzero effort for the check]
    tolerance_rel: 1e-6  // [ASSUMED — relative tolerance for the torque-balance check]
};

function runStartupSelfCheck() {
    // Save and later restore the module-level configuration this check mutates
    let saved = { omega, delta, driveFrac, drivePhaseEnd, F_peak: F_peak.slice() };

    omega = 2 * Math.PI * (SELF_CHECK.spm / 60);
    delta = omega * SELF_CHECK.offset_s;
    driveFrac = getDriveFraction(SELF_CHECK.spm);
    drivePhaseEnd = 2 * Math.PI * driveFrac;
    F_peak = [SELF_CHECK.F_peak_test, SELF_CHECK.F_peak_test];

    // CHECK 1: between the catches, only rower 1 is pulling
    let phi_1 = delta * SELF_CHECK.catchGapFrac;   // after rower 1's catch…
    let phi_2 = phi_1 - delta;                     // …before rower 2's
    let Fh1 = handleForce(0, phi_1);
    let Fh2 = handleForce(1, phi_2);
    let ok1 = Fh1 > 0 && Fh2 === 0;

    // CHECK 2: torque balance at the solved point
    let phiMid = drivePhaseEnd * SELF_CHECK.drivePointFrac;
    let b = rowerBlade(0, phiMid, SELF_CHECK.theta_test, SELF_CHECK.v_test);
    let Freq = b.Fh * OAR.L_in / OAR.L_out;
    let Fperp = bladeHydro(SELF_CHECK.theta_test, b.thetaDot, SELF_CHECK.v_test).Fperp;
    let ok2 = b.thetaDot === 0
        ? Fperp >= Freq   // clamp branch (see solveThetaDot)
        : Math.abs(Fperp - Freq) / Freq < SELF_CHECK.tolerance_rel;

    omega = saved.omega; delta = saved.delta;
    driveFrac = saved.driveFrac; drivePhaseEnd = saved.drivePhaseEnd;
    F_peak = saved.F_peak;

    if (!ok1 || !ok2) {
        console.error("SELF-CHECK FAILED:");
        if (!ok1) console.error(`  (sign convention) with delta > 0 rower 2 must catch LATER: between the catches expected Fh1 > 0 and Fh2 = 0, got Fh1=${Fh1} Fh2=${Fh2}`);
        if (!ok2) console.error(`  (item the prescribed-force closure) torque balance violated: Fperp=${Fperp} required=${Freq} thetaDot=${b.thetaDot}`);
        process.exit(1);
    }
}
// ========================== CLI PARAMETER HANDLING =========================
//
// Usage:
//   node simulator_final.js                        # full sweep with defaults
//   node simulator_final.js --weight 75            # set rower mass (kg each)
//   node simulator_final.js --offset 100           # single offset (ms)
//   node simulator_final.js --weight 75 --offset 100
//   node simulator_final.js --help

function parseArgs() {
    let args = { weight: 80, weight2: null, offset: null, spm: null, calibrate: false, sensitivity: false };
    let argv = process.argv.slice(2);

    for (let i = 0; i < argv.length; i++) {
        if (argv[i] === '--weight' || argv[i] === '-w') {
            args.weight = parseFloat(argv[++i]);
            if (isNaN(args.weight) || args.weight < 40 || args.weight > 120) {
                console.error("Error: --weight must be between 40 and 120 kg");
                process.exit(1);
            }
        } else if (argv[i] === '--weight2') {
            // the single mass block: unequal rowers (rower 2's mass; defaults to --weight)
            args.weight2 = parseFloat(argv[++i]);
            if (isNaN(args.weight2) || args.weight2 < 40 || args.weight2 > 120) {
                console.error("Error: --weight2 must be between 40 and 120 kg");
                process.exit(1);
            }
        } else if (argv[i] === '--offset' || argv[i] === '-o') {
            args.offset = parseFloat(argv[++i]);
            if (isNaN(args.offset) || args.offset < 0 || args.offset > 2000) {
                // Cap raised from 500 ms (31 Aug 2026) so the ANTIPHASE case — half a
                // cycle, e.g. 1000 ms at 30 SPM — can be run as a validation point
                // against Boucher, Labbé & Clanet 2017 (on-water robots) and
                // de Brouwer, de Poel & Hofmijster 2013 (ergometers on slides).
                console.error("Error: --offset must be between 0 and 2000 ms");
                process.exit(1);
            }
        } else if (argv[i] === '--spm' || argv[i] === '-s') {
            args.spm = parseFloat(argv[++i]);
            if (isNaN(args.spm) || args.spm < 16 || args.spm > 44) {
                console.error("Error: --spm must be between 16 and 44");
                process.exit(1);
            }
        } else if (argv[i] === '--constcf') {
            // The drag rebuild: constant skin-friction mode (c = BOAT.Cf_const) for
            // diffing against the ITTC-1957 speed-dependent line.
            constCfMode = true;
        } else if (argv[i] === '--dragscale') {
            // The drag rebuild sweep: multiplies the computed friction term so the
            // published anchor disagreement (equivalents ~3.5/4.7/6.0) can be
            // swept; report how the sync-async RATIO moves.
            dragScale = parseFloat(argv[++i]);
            if (isNaN(dragScale) || dragScale <= 0 || dragScale > 3) {
                console.error("Error: --dragscale must be in (0, 3]");
                process.exit(1);
            }
        } else if (argv[i] === '--sensitivity') {
            // The sensitivity sweep: ±20% parameter sensitivity table on the headline loss.
            args.sensitivity = true;
        } else if (argv[i] === '--dt') {
            // Review: integration timestep override for convergence studies
            SIM.dt = parseFloat(argv[++i]);
            if (isNaN(SIM.dt) || SIM.dt < 0.0001 || SIM.dt > 0.005) {
                console.error("Error: --dt must be in [0.0001, 0.005] s");
                process.exit(1);
            }
        } else if (argv[i] === '--trim') {
            // The static-trim mode: quasi-static trim → wetted-area ratio (Formaggia static limit)
            trimMode = true;
        } else if (argv[i] === '--pitch') {
            // The heave/pitch mode: heave+pitch dynamics (symmetry-plane) → S ratio + diagnostics
            pitchMode = true;
        } else if (argv[i] === '--zeta') {
            // The heave/pitch mode: heave/pitch damping fraction of critical [ASSUMED default 0.3]
            HULL2D.zeta = parseFloat(argv[++i]);
            if (isNaN(HULL2D.zeta) || HULL2D.zeta <= 0 || HULL2D.zeta > 1) {
                console.error("Error: --zeta must be in (0, 1]");
                process.exit(1);
            }
        } else if (argv[i] === '--rawpower') {
            // the prescribed-force closure: disable equal-power renormalization — compare at fixed
            // handle-force profiles instead of fixed per-rower average power.
            equalPowerMode = false;
        } else if (argv[i] === '--calibrate') {
            // The blade-slip closure calibration protocol: reproduce EFFORT.F_peak_ref.
            args.calibrate = true;
        } else if (argv[i] === '--help' || argv[i] === '-h') {
            console.log("Rowing Asynchronicity Simulator — Double Scull (2x)");
            console.log("");
            console.log("Usage: node simulator_final.js [options]");
            console.log("");
            console.log("Options:");
            console.log("  --weight, -w <kg>    Rower mass in kg (default: 80, range: 40-120)");
            console.log("  --weight2 <kg>       Rower 2 mass if different (unequal crews)");
            console.log("  --offset, -o <ms>    Catch timing offset in ms (default: full 0-200 sweep)");
            console.log("  --spm, -s <rate>     Stroke rate (default: all essay phases, range: 16-44)");
            console.log("  --rawpower           Disable equal-power renormalization: compare at fixed handle-force profiles");
            console.log("  --calibrate          Run the effort calibration protocol; prints F_peak_ref to hardcode");
            console.log("  --constcf            Constant skin friction c=0.00225 (Sliasas & Tullis 2009) instead of the ITTC-1957 line");
            console.log("  --dragscale <f>      Multiply the friction term (published-drag-anchor sweep: 0.93/1.245/1.59 ~ k 3.5/4.7/6.0 N s2/m2)");
            console.log("  --sensitivity        +/-20% one-at-a-time parameter sensitivity sweep on the 100ms/30SPM headline");
            console.log("  --trim               Quasi-static trim drag modulation (Formaggia et al. 2009, static limit)");
            console.log("  --pitch              Heave+pitch dynamics (symmetry-plane, Formaggia et al. 2009); exclusive with --trim");
            console.log("  --zeta <f>           Heave/pitch damping, fraction of critical (default 0.3 [ASSUMED], range (0,1])");
            console.log("  --dt <s>             Integration timestep override for convergence studies (default 0.001)");
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

// ========================== CONFIGURE BOAT FOR ROWER WEIGHT ================
// Recalculates all mass-dependent parameters.
// Effort and hull drag scale with rower mass to keep velocity realistic.
//   - F_peak ∝ mass (bigger rower → more force, roughly linear) [ASSUMED]
//   - Drag scales with crew weight by DERIVATION: S = 2.30·√(V·L) recomputed
//     from displacement at evaluation time (the drag rebuild)
// Baseline: 80 kg rower → EFFORT.F_peak_ref

function configureBoat(rowerMass1, rowerMass2) {
    let baseline = 80;  // kg, calibration point (see RIG_REF / EFFORT)
    if (rowerMass2 === undefined) rowerMass2 = rowerMass1;

    // the single mass block: the MASS block is the only place masses live; configureBoat
    // only sets the per-rower inputs. (the 30 Aug 2026 code audit: no rounding.)
    MASS.rower_1 = rowerMass1;
    MASS.rower_2 = rowerMass2;

    // Effort scale (the prescribed-force closure): F_peak ∝ rower mass [ASSUMED — carried over from the
    // old F_0 ∝ mass rule]; per-rower so unequal rowers and power renorm work.
    F_peak[0] = EFFORT.F_peak_ref * (rowerMass1 / baseline);
    F_peak[1] = EFFORT.F_peak_ref * (rowerMass2 / baseline);

    // (The drag rebuild: drag needs no per-weight scaling here — wettedArea() computes
    // S from displacement at evaluation time, so --weight scales drag by
    // derivation. The old ^(2/3) k_hull/A_wetted hack is deleted.)
}
// ========================== DISPLAY SETTINGS ===============================

function displaySettings() {
    // The drag rebuild: report the COMPUTED drag pipeline (k_hull no longer exists).
    let v_show = 4.3;   // m/s, representative pace for the equivalent-k display only
    let Re_show = v_show * BOAT.L / CONSTANTS.nu_water;
    let Cf_show = 0.075 / Math.pow(Math.log10(Re_show) - 2, 2);
    let k_fric_eq = 0.5 * CONSTANTS.rho_water * Cf_show * wettedArea();

    console.log("============================================================");
    console.log("  ROWING ASYNCHRONICITY SIMULATOR");
    console.log("  Double Scull (2x) — 1D Newtonian Model");
    console.log("============================================================\n");
    console.log("  MASS");
    console.log(`    Hull:       ${MASS.shell} kg shell (FISA minimum 2x) + 2 × ${MASS.seat_hw} kg seat hardware`);
    console.log(`    Rowers:     ${MASS.rower_1} + ${MASS.rower_2} = ${MASS.crew} kg`);
    console.log(`    Total:      ${MASS.total} kg`);
    console.log(`    Effective:  ${(MASS.total*(1+MASS.k_added)).toFixed(1)} kg (added mass ${(MASS.k_added*100).toFixed(0)}%)`);
    console.log(`    Translating: ${(MASS.frac_translating*100).toFixed(1)}% of body (de Leva) + seat hw; arms ${(MASS.frac_arms*100).toFixed(1)}% on handle path`);
    console.log();
    console.log("  PROPULSION");
    console.log(`    F_peak:     ${F_peak[0].toFixed(0)} N peak handle force/rower (both hands; the one calibrated scale)`);
    console.log(`    Gearing:    L_in=${OAR.L_in} m, L_out=${OAR.L_out.toFixed(2)} m (oar ${OAR.length} m)`);
    console.log(`    Blade:      2 × ${(OAR.A_blade*1e4).toFixed(0)} cm²/rower, Al=${OAR.Al}, Ad=${OAR.Ad} (Caplan & Gardner 2007a via Atkinson)`);
    console.log(`    Arc input:  catch ${(OAR.theta_catch/DEG).toFixed(0)}°; release/arc EMERGENT`);
    console.log(`    Profile:    Beta(α=${BOAT.alpha}, β=${BOAT.beta}), peak at ${((BOAT.alpha-1)/(BOAT.alpha+BOAT.beta-2)*100).toFixed(0)}% through drive`);
    console.log();
    console.log("  DRAG");
    console.log(`    Friction:   ½ρ·Cf(Re)·S·v|v| — ITTC-1957 line${constCfMode ? ' [OVERRIDDEN: --constcf c=' + BOAT.Cf_const + ']' : ''}, ρ=${CONSTANTS.rho_water}`);
    console.log(`    S (wetted): ${wettedArea().toFixed(2)} m² = ${BOAT.S_coeff}·√(V·L), from displacement ${MASS.total} kg`);
    console.log(`    Cf @4.3m/s: ${Cf_show.toFixed(5)} (Re=${(Re_show/1e6).toFixed(1)}e6) → friction-k equivalent ${k_fric_eq.toFixed(2)} Ns²/m²${dragScale !== 1 ? ' × dragScale ' + dragScale.toFixed(3) : ''}`);
    console.log(`    Remainder:  k_wave=${BOAT.k_wave} (Froude, multiplicative), k_air=${BOAT.k_air} Ns²/m² — the ~15-20% non-friction share`);
    console.log();
    console.log("  KINEMATICS");
    console.log(`    Travels:    seat ${MASS.travel_seat} m, trunk ${MASS.travel_trunk} m, thighs ${MASS.travel_thighs} m (totals; Kleshnev Kleshnev 2014, BioRow newsletter route)`);
    console.log(`    Hull:       ${BOAT.L} m × ${BOAT.B} m (Filippi F17 (Filippi Lido specifications)); S computed per run`);
    if (trimMode || pitchMode) console.log(`    Attitude:   ${trimMode ? 'QUASI-STATIC TRIM' : 'HEAVE+PITCH DYNAMICS (ζ=' + HULL2D.zeta + ' [ASSUMED])'} — S-ratio drag modulation, Formaggia et al. 2009 formulation`);
    console.log();
    console.log("  INTEGRATION");
    console.log(`    Method:     RK4, dt=${SIM.dt*1000}ms, warmup: stroke-mean convergence < ${WARMUP.settleThreshold_mps} m/s (cap ${WARMUP.maxStrokes} strokes), window=${MEASURE.nCycles} complete cycles`);
    console.log("------------------------------------------------------------");
}

// ========================== RUN ============================================

const OFFSETS_MS = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200];

const ESSAY_SPMS = [
    { label: "Controlled test (low)",  spm: 20 },
    { label: "Controlled test (high)", spm: 22 },
    { label: "Body (low)",             spm: 26 },
    { label: "Body (high)",            spm: 28 },
    { label: "Start",                  spm: 30 },
    { label: "Sprint",                 spm: 32 }
];

runStartupSelfCheck();

let args = parseArgs();
if (trimMode && pitchMode) {
    console.error("Error: --trim and --pitch are exclusive (--pitch subsumes the static channel)");
    process.exit(1);
}
configureBoat(args.weight, args.weight2 === null ? args.weight : args.weight2);
displaySettings();

if (args.sensitivity) {
    // ================= PARAMETER SENSITIVITY TABLE (the sensitivity sweep) =================
    // Each uncertain parameter varied ±20% (k_added additionally at its recorded
    // sweep anchors 0.03/0.10); the change in the HEADLINE LOSS (equal-power
    // speed loss at 100 ms / 30 SPM) is printed. Converts every remaining
    // "where does that constant come from" into "demonstrated (in)sensitive".
    // Our numbers only — no thresholds, no verdicts (the blinded output table blinding).
    const SENS = { spm: 30, offset_s: 0.100, frac: 0.20 };   // [the recorded sensitivity-sweep design]
    let headline = () => {
        let p = simulatePair(SENS.spm, SENS.offset_s);
        return (p.sync.metrics.v_mean - p.off.metrics.v_mean) / p.sync.metrics.v_mean * 100;
    };
    let base = headline();
    console.log(`\n=== SENSITIVITY: headline loss @ ${SENS.offset_s*1000} ms / ${SENS.spm} SPM ===`);
    console.log(`  Baseline loss: ${base.toFixed(4)}%\n`);
    console.log("  Parameter                −20% → loss      +20% → loss");
    let rows = [
        ["MASS.k_added",       (f) => { let o = MASS.k_added;       MASS.k_added = o * f;       return () => MASS.k_added = o; }],
        ["translating masses", (f) => { let o = [MASS.frac_pelvis, MASS.frac_trunk_rest, MASS.frac_thighs];
                                        MASS.frac_pelvis = o[0]*f; MASS.frac_trunk_rest = o[1]*f; MASS.frac_thighs = o[2]*f;
                                        return () => { MASS.frac_pelvis = o[0]; MASS.frac_trunk_rest = o[1]; MASS.frac_thighs = o[2]; }; }],
        ["travel distances",   (f) => { let o = [MASS.travel_seat, MASS.travel_trunk, MASS.travel_thighs];
                                        MASS.travel_seat = o[0]*f; MASS.travel_trunk = o[1]*f; MASS.travel_thighs = o[2]*f;
                                        return () => { MASS.travel_seat = o[0]; MASS.travel_trunk = o[1]; MASS.travel_thighs = o[2]; }; }],
        ["OAR.Al (lift amp)",  (f) => { let o = OAR.Al;             OAR.Al = o * f;             return () => OAR.Al = o; }],
        ["OAR.Ad (drag amp)",  (f) => { let o = OAR.Ad;             OAR.Ad = o * f;             return () => OAR.Ad = o; }],
        ["OAR.blade_fill",     (f) => { let o = OAR.blade_fill;     OAR.blade_fill = o * f;     return () => OAR.blade_fill = o; }],
        ["BOAT.k_wave",        (f) => { let o = BOAT.k_wave;        BOAT.k_wave = o * f;        return () => BOAT.k_wave = o; }]
    ];
    for (let [name, apply] of rows) {
        let restore = apply(1 - SENS.frac); let lo = headline(); restore();
        restore = apply(1 + SENS.frac);     let hi = headline(); restore();
        console.log(`  ${name.padEnd(22)} ${lo.toFixed(4)}% (Δ${(lo-base).toFixed(4)})   ${hi.toFixed(4)}% (Δ${(hi-base).toFixed(4)})`);
    }
    // k_added recorded anchor values (see MASS block)
    {
        let o = MASS.k_added;
        MASS.k_added = 0.03; let a = headline();
        MASS.k_added = 0.10; let b = headline();
        MASS.k_added = o;
        console.log(`  k_added anchors        0.03: ${a.toFixed(4)}%        0.10: ${b.toFixed(4)}%`);
    }
    process.exit(0);
}

if (args.calibrate) {
    // CALIBRATION PROTOCOL (the blade-slip closure): secant iteration on F_peak_ref until the
    // synchronized 2x (two 80 kg rowers) at EFFORT.cal_SPM hits EFFORT.v_cal_target.
    // Prints the value to hardcode into EFFORT.F_peak_ref. Deterministic.
    let lo = EFFORT.F_peak_ref * 0.5, hi = EFFORT.F_peak_ref * 2;   // bracket [ASSUMED — numerical]
    let f = (Fp) => { F_peak = [Fp, Fp]; return simulate(EFFORT.cal_SPM, 0).metrics.v_mean - EFFORT.v_cal_target; };
    for (let i = 0; i < 40; i++) {   // bisection iterations [ASSUMED — numerical]
        let mid = 0.5 * (lo + hi);
        if (f(mid) < 0) lo = mid; else hi = mid;
        if (hi - lo < 0.01) break;   // N resolution [ASSUMED — numerical]
    }
    let Fcal = 0.5 * (lo + hi);
    F_peak = [Fcal, Fcal];
    let vchk = simulate(EFFORT.cal_SPM, 0).metrics.v_mean;
    console.log(`\nCALIBRATION: F_peak_ref = ${Fcal.toFixed(1)} N gives v_mean = ${vchk.toFixed(4)} m/s at ${EFFORT.cal_SPM} SPM (target ${EFFORT.v_cal_target}).`);
    console.log("Hardcode this into EFFORT.F_peak_ref with a [DERIVED] label and the run date.");
    process.exit(0);
}


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
        let pair = simulatePair(spm, args.offset / 1000);
        let syncRun = pair.sync, offsetRun = pair.off;
        let loss = (syncRun.metrics.v_mean - offsetRun.metrics.v_mean) / syncRun.metrics.v_mean * 100;
        let dvr  = offsetRun.metrics.Delta_v / offsetRun.metrics.v_mean * 100;
        let df   = getDriveFraction(spm) * 100;

        console.log(`  ${label.padEnd(21)} | ${String(spm).padStart(3)} | ${df.toFixed(1).padStart(5)}% | ${syncRun.metrics.v_mean.toFixed(3)}  | ${offsetRun.metrics.v_mean.toFixed(3)}   | +${loss.toFixed(3).padStart(6)}% | ${dvr.toFixed(1)}%`);
    }

    // Race impact for this offset
    // ASSUMPTION: penalties use steady-state phase speeds; standing-start
    // acceleration is not modeled.
    let s30r = simulate(30, 0), s27r = simulate(27, 0), s32r = simulate(32, 0);
    let a30 = simulatePair(30, args.offset/1000, s30r).off;
    let a27 = simulatePair(27, args.offset/1000, s27r).off;
    let a32 = simulatePair(32, args.offset/1000, s32r).off;
    let s30 = s30r, s27 = s27r, s32 = s32r;
    let pen_start  = 500  / a30.metrics.v_mean - 500  / s30.metrics.v_mean;
    let pen_body   = 3500 / a27.metrics.v_mean - 3500 / s27.metrics.v_mean;
    let pen_sprint = 1000 / a32.metrics.v_mean - 1000 / s32.metrics.v_mean;
    console.log(`\n  5K race penalty: Start +${pen_start.toFixed(1)}s | Body +${pen_body.toFixed(1)}s | Sprint +${pen_sprint.toFixed(1)}s | TOTAL +${(pen_start+pen_body+pen_sprint).toFixed(1)}s`);

    // The blinded output table: blinded table for the requested (or 30 SPM) condition
    printBlindedTable(s30, a30);

} else {
    // ── Full offset table (one or all SPMs) ──
    for (let { label, spm } of spmList) {
        console.log(`\n=== ${label}: ${spm} SPM (drive fraction: ${(getDriveFraction(spm)*100).toFixed(1)}%) ===\n`);
        let result = generateSpeedLossTable(spm, OFFSETS_MS);

        console.log(`  Synchronized: v̄ = ${result.sync_v.toFixed(4)} m/s, Δv = ${result.sync_dv.toFixed(4)} m/s, Δv/v̄ = ${(result.sync_dv/result.sync_v*100).toFixed(1)}%`);
        console.log();
        console.log("  Δt(ms) |  δ(rad)  | v̄(m/s)  | EqPwr Loss | Raw Loss  | Δv/v̄   | Δv change");
        console.log("  -------|----------|---------|------------|-----------|--------|----------");

        for (let row of result.table) {
            console.log(`  ${String(row.dt_ms).padStart(5)}  | ${row.delta_rad.toFixed(4).padStart(7)}  | ${row.v_mean.toFixed(4)} | ${row.speed_loss_pct >= 0 ? '+' : ''}${row.speed_loss_pct.toFixed(4)}%   | ${row.raw_loss_pct >= 0 ? '+' : ''}${row.raw_loss_pct.toFixed(4)}% | ${row.dv_over_vmean.toFixed(1).padStart(5)}% | ${row.dv_change_pct >= 0 ? '+' : ''}${row.dv_change_pct.toFixed(2)}%`);
        }
        console.log(NOISE_FLOOR_NOTE);
    }

    // Race impact and SyncRow table only in full run (no --spm filter)
    if (args.spm === null) {
        // ASSUMPTION: penalties use steady-state phase speeds; standing-start
        // acceleration is not modeled.
        console.log("\n\n=== RACE IMPACT: 5000m HEAD RACE ===\n");
        console.log("  Start 500m @ 30 SPM | Body 3500m @ 27 SPM | Sprint 1000m @ 32 SPM\n");

        let s30r = simulate(30, 0), s27r = simulate(27, 0), s32r = simulate(32, 0);
        for (let dt_ms of [20, 50, 100, 150, 200]) {
            let s30 = s30r, a30 = simulatePair(30, dt_ms/1000, s30r).off;
            let s27 = s27r, a27 = simulatePair(27, dt_ms/1000, s27r).off;
            let s32 = s32r, a32 = simulatePair(32, dt_ms/1000, s32r).off;
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
            let a = dt_ms === 0 ? s30 : simulatePair(30, dt_ms / 1000, s30).off;
            // the 30 Aug 2026 code audit: delta comes from the simulation result (derived once from omega
            // in simulate()), not re-derived here. The old line hardcoded 0.5 Hz
            // (2π·0.5·Δt) — correct at 30 SPM only, and a second independent
            // derivation of delta.
            let dv = a.delta;
            let loss = (s30.metrics.v_mean - a.metrics.v_mean) / s30.metrics.v_mean;
            // Raw-power companion for the comment (review: design-dependence must ship with the number)
            let rawTxt = "";
            if (dt_ms !== 0) {
                let ra = simulatePair(30, dt_ms / 1000, s30, false).off;
                rawTxt = `, raw ${((s30.metrics.v_mean - ra.metrics.v_mean) / s30.metrics.v_mean * 100).toFixed(3)}%`;
            }
            console.log(`    ${dv.toFixed(3)}: ${loss.toFixed(6)},  // ${dt_ms}ms → ${(loss*100).toFixed(3)}% equal-power${rawTxt}`);
        }
        console.log("};");
        console.log("// Values are EQUAL-POWER losses (per-rower mean handle power held equal to the synchronized baseline); raw-power companions in comments.");
        console.log("//" + NOISE_FLOOR_NOTE);
    }
}

// ============================================================================
// REFERENCES
// Every source used anywhere in this file, in full. Inline text refers to
// these by author and year. Web sources give the organization and the page
// name; retrieval date is the last verification (30 Aug 2026). Where a source
// is a self-published practitioner document (BioRow newsletters, Atkinson's
// ROWING software notes, van Holst's site) it is cited as such — these are
// the primary sources for their numbers, not peer-reviewed literature, and the
// text says so where it matters.
//
// Atkinson W. ROWING — Rowing Computer Research (simulation software and
//   documentation, incl. tabulated wetted-surface data for measured racing
//   hulls and the blade-force closure method). http://www.atkinsopht.com/row/
//   (self-published; verified 30 Aug 2026).
// Boucher J-P, Labbé R, Clanet C. Row bots. Physics Today 2017;70(6):82-83.
//   (Robotic crews: synchronous rowing faster, asynchronous smoother.)
// Cabrera D, Ruina A, Kleshnev V. A simple 1+ dimensional model of rowing
//   mimics observed forces and motions. Human Movement Science
//   2006;25(2):192-220.
// Caplan N, Gardner TN. A fluid dynamic investigation of the Big Blade and
//   Macon oar blade designs in rowing propulsion. Journal of Sports Sciences
//   2007;25(6):643-650.  [cited as Caplan & Gardner 2007a — blade coefficients]
// Caplan N, Gardner TN. A mathematical model of the oar blade–water
//   interaction in rowing. Journal of Sports Sciences 2007;25(9):1025-1034.
//   [cited as Caplan & Gardner 2007b — single-mass boat model]
// Concept2 Inc. Sculling oars: specifications and blade dimensions
//   (Smoothie2 Plain Edge scull blade 46 × 21.5 cm; scull lengths 284-290 cm).
//   https://www.concept2.com/oars/sculls (verified 30 Aug 2026).
// de Brouwer AJ, de Poel HJ, Hofmijster MJ. Don't rock the boat: how antiphase
//   crew coordination affects rowing. PLoS ONE 2013;8(1):e54996. (Ergometers
//   on slides: antiphase cuts fluctuation losses and raises useful power.)
// Cuijpers LS, Zaal FTJM, de Poel HJ. Rowing crew coordination dynamics at
//   increasing stroke rates. PLoS ONE 2015;10(7):e0133527.
// de Leva P. Adjustments to Zatsiorsky-Seluyanov's segment inertia
//   parameters. Journal of Biomechanics 1996;29(9):1223-1230.
// Filippi Lido S.r.l. Double scull model range specifications (F13, F17,
//   F51: length, beam, crew weight class). https://www.filippiboats.com
//   (verified 30 Aug 2026).
// Formaggia L, Miglio E, Mola A, Montano A. A model for the dynamics of rowing
//   boats. International Journal for Numerical Methods in Fluids
//   2009;61(2):119-143. doi:10.1002/fld.1940 (online 2008).
// Grift EJ, Vijayaragavan NB, Tummers MJ, Westerweel J. Drag force on an
//   accelerating submerged plate. Journal of Fluid Mechanics 2019;866:369-398.
//   doi:10.1017/jfm.2019.102. (Rowing-blade-inspired: drag during acceleration
//   is not captured by a single added-mass coefficient — why this file's
//   quasi-static blade coefficients are a stated limitation.)
// Greidanus AJ, Delfos R, Westerweel J. Drag and power-loss in rowing due to
//   velocity fluctuations. Procedia Engineering 2016;147:317-323.
// ITTC. 1957 model-ship correlation line, Cf = 0.075/(log10 Re − 2)².
//   Proceedings of the 8th International Towing Tank Conference, Madrid, 1957.
// Kleshnev V. Rowing Biomechanics Newsletter (BioRow), self-published,
//   www.biorow.com. Issues used: "Rigging" (2007, No. 3); "Effect of stroke
//   rate on rowing technique"; "Amplitude and power of body segments" (2014,
//   n = 5,437 samples); "The biomechanics of the recovery phase" (n = 25,658);
//   release-slip framework (2018); "Averaged biomechanical curves" (2024,
//   n > 50,000); Catch Factor measurements (2025, No. 4).
// Kleshnev V. The Biomechanics of Rowing. Ramsbury: Crowood Press; 2016.
//   (Cited only as the book to which an earlier version of this code falsely
//   attributed a beta-distribution force profile; it contains none.)
// Kleshnev V. Biomechanics of rowing. In: Nolte V, ed. Rowing Faster. 2nd ed.
//   Champaign, IL: Human Kinetics; 2011: chapter 9 (Table 9.2).
// Laschowski B, Nolte V, Adamovsky M, Alexander R. The effects of oar-shaft
//   stiffness and length on rowing biomechanics. Proceedings of the IMechE,
//   Part P: Journal of Sports Engineering and Technology 2015;229(4):239-247.
//   (Oar mass/stiffness data — the source for the shaft-inertia term if the
//   oar-inertia extension is ever built.)
// Labbé R, Boucher J-P, Clanet C, Benzaquen M. Physics of rowing oars. New
//   Journal of Physics 2019;21:093050. (Imposed-FORCE oar framework "closer
//   to human constraints" — independent precedent for this file's closure.)
// Lazauskas L. A performance prediction model for rowing races. Technical
//   report, Department of Applied Mathematics, University of Adelaide; 1997.
// Sliasas A, Tullis S. Numerical modelling of rowing blade hydrodynamics.
//   Sports Engineering 2009;12(1):31-40.
// Sliasas A, Tullis S. The dynamic flow behaviour of an oar blade in motion
//   using a hydrodynamics-based shell-velocity-coupled model of a rowing
//   stroke. Proceedings of the IMechE, Part P: Journal of Sports Engineering
//   and Technology 2010;224(1):9-24.
// Soper C, Hume PA. Towards an ideal rowing technique for performance: the
//   contributions from biomechanics. Sports Medicine 2004;34(12):825-848.
// van Holst M. On rowing (simulation notes and drag estimates).
//   http://home.hccnet.nl/m.holst/ (self-published; verified 30 Aug 2026).
// Vespoli M, Nelson B, Scragg C. Eight man rowing shell. US Patent 5,474,008;
//   1995. (Wetted-surface data for this and two related Vespoli-lineage hulls,
//   US Patents 5,188,048 and 5,279,239, as tabulated by Atkinson.)
// World Rowing (FISA). Rules of Racing, Rule 32 and Bye-Law: minimum boat
//   weights (2x = 27 kg). https://worldrowing.com (verified 30 Aug 2026).
// ============================================================================
