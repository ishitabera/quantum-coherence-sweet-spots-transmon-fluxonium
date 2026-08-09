# Lab Notebook

## June 30, 2026 — Day 1: Pipeline validation (PASSED)
- Transmon (EJ=15, EC=0.3): freq 5.683 GHz, anharm -0.345 GHz — matches sqrt(8*EJ*EC)-EC theory.
- Fluxonium (EJ=8.9, EC=2.5, EL=0.5): energy-vs-flux flat at flux=0.5 — sweet spot confirmed.
- Transmon coherence: T1 28.0 us, T2 56.0 us  (T2 = 2*T1).
- Fluxonium coherence: T1 771.1 us, T2 1542.2 us = 1.54 ms  (T2 = 2*T1).
- Note: T2 = 2*T1 in both — transmon charge-insensitivity and fluxonium sweet spot make
  dephasing negligible, so relaxation sets the ceiling. Core physics of the project.
- scqubits returns times in NANOSECONDS (convert to us with *1e-3).
- Transmon has no flux-noise channel; fluxonium does (checked supported_noise_channels()).

## July 5, 2026 — Week 2: Transmon EJ/EC sweep (DONE)
- 800-pt sweep at ng=0 saved (data/transmon_full_sweep.csv). Relaxation-limited EVERYWHERE
  -> at the charge sweet spot (ng=0), dielectric loss is the ceiling across the whole plane. (Q2 finding)
- Caveat: best-T2 corner (EJ/EC=0.5) is charge-qubit regime, not a real transmon target.
- Off sweet spot (ng=0.25): charge noise dominates at LOW EJ/EC, relaxation at HIGH EJ/EC
  -> the dominant-channel CROSSOVER. This is the Q2 story. (data/transmon_sweep_ng025.csv)
- T2 heatmap: figures/transmon_T2_map.png
- - Ran TWO ng values on purpose: ng=0.0 is the charge sweet spot (charge noise suppressed
  -> relaxation-limited everywhere); ng=0.25 is off the sweet spot, which wakes up charge
  noise at low EJ/EC and reveals the dominant-channel crossover at EJ/EC ~ 25. (Q2)
- Figures: transmon_T2_map.png (Q1), transmon_dominant_map.png (Q2).

## July 10, 2026 — Week 3: Fluxonium flux-bias sweep (DONE)
- Extended pipeline to fluxonium: T1 = dielectric + quasiparticle; Tphi = flux + crit-current.
- Sanity at flux=0.5: T1 771.1 us (dielectric-limited; quasiparticle 3e30 us = negligible),
  T2 1136.9 us. Note T2 < 2*T1 here because critical-current dephasing (4325 us) also bites.
- Flux sweep 0.30->0.50 (101 pts) saved. T2 peaks sharply at flux=0.5 sweet spot.
- SENSITIVITY (Q1/Q3): T2 = 1136.9 us at sweet spot, falls 10x by flux=0.4980 — only
  0.0020 flux-quanta away. The sweet spot is a razor-thin peak.
- Dominant channel (Q2): flux noise dominates 100/101 points; relaxation dominates ONLY at
  exactly flux=0.5. Crossover is essentially instantaneous — the moment you detune, flux noise wins.
- Design map (Q1): best T2 = 5578.1 us at EJ=15.0, EL=0.2 (high EJ, low EL). = 5.6 ms.
- Figures: fluxonium_flux_falloff.png, fluxonium_channels_vs_flux.png, fluxonium_design_map.png

## July 20, 2026 — Week 4: Cross-architecture analysis & robustness (DONE)
FINDINGS:
1. Peak coherence. Fluxonium reaches ~1137 µs at its sweet spot (up to 5578 µs / 5.6 ms for
   the best EJ-EL design) vs transmon ~16.7 µs. Physics: the flux sweet spot nulls 1/f flux
   dephasing to first order, while the transmon stays dielectric-loss limited.
   Implication: fluxonium wins outright when maximum coherence is the goal.

2. Robustness (the trade-off). Transmon T2 varies only 1.01x across its ENTIRE charge-bias
   range (16.56–16.66 µs); fluxonium varies 528x across just flux 0.45–0.50. Physics: transmon
   charge dispersion is exponentially suppressed at EJ/EC=50; fluxonium flux sensitivity grows
   steeply off Φ/Φ0=0.5. Implication: transmon is dramatically more forgiving of bias drift and
   fabrication spread — why IBM/Google build on it.

3. Robustness to MODEL uncertainty. Across 60 Monte Carlo draws of noise amplitudes over
   literature ranges, fluxonium beats transmon in 100% of samples, ratio median 216x,
   90% CI [51.8x, 371.1x]. The peak-coherence conclusion is robust to input uncertainty.

4. Sweet-spot scaling. Fluxonium T2 ~ (flux detuning)^-0.68 over 0.001–0.05; 10x drop within
   0.002 flux quanta (Week 3). Physics: dephasing rate ~ (dω/dΦ)^2, which climbs sharply off the
   extremum. (Note: pure-flux-limited theory predicts ~-1; the shallower -0.68 reflects relaxation
   and crit-current capping T2 near the spot — worth stating honestly.)

5. Dominant channel (Q2). Transmon: charge noise at low EJ/EC off sweet spot → relaxation in the
   transmon regime (crossover ~EJ/EC=25). Fluxonium: relaxation only exactly at Φ=0.5 → flux noise
   everywhere else.