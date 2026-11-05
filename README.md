# Coherence Sweet Spots Across Superconducting Qubit Architectures

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23988521.svg)](https://doi.org/10.5281/zenodo.23988521)

The version used for the preprint is release **v1.0**, archived at Zenodo: [doi:10.5281/zenodo.23988521](https://doi.org/10.5281/zenodo.23988521) (all versions: [10.5281/zenodo.23988520](https://doi.org/10.5281/zenodo.23988520)).

A parameter-sensitivity study of the **transmon** and **fluxonium** under realistic noise
channels, using [scqubits](https://scqubits.readthedocs.io). All results are model-based
estimates computed on a laptop; the value is in relative trends across parameters and
architectures, not in absolute coherence times.

## Research questions

| | Question |
|---|---|
| **Q1** | Across the design parameter space (EJ, EC for the transmon; EJ, EC, EL, flux for the fluxonium), how do T1, T2 and Tφ vary, and where are the coherence sweet spots? |
| **Q2** | Which noise channel dominates where, and how does the dominant channel shift across the space? |
| **Q3** | Does the fluxonium's flux sweet spot buy more coherence robustness than the transmon's charge insensitivity, and at what design cost? |

## Headline results

All numbers below are traceable to a file and a filter; the ledger in `05_summary.ipynb` §2
recomputes each one from its named source and asserts it.

> **Convention change, 2026-11-04.** Critical-current noise is now read as **fractional**
> (A_cc = δI_c/I_c, passed to scqubits as A_cc·E_J). scqubits applies A_cc to the dimensionless
> ∂H/∂E_J, i.e. as an absolute E_J fluctuation in GHz, which made the channel 8.9× (fluxonium) to
> 13× (transmon) too weak. Every number below is for the fractional reading; the earlier numbers
> (e.g. fluxonium T2 = 1130 µs, raw ratio 18.1×, matched-frequency median 0.86×) are reproducible
> with `cc_convention="absolute"` and are listed in the paper's §3.1.

| Metric | Transmon | Fluxonium |
|---|---|---|
| Operating point | EJ=13.0, EC=0.26 (EJ/EC=50), f01=4.92 GHz | EJ=8.9, EC=2.5, EL=0.5, f01=0.363 GHz |
| T2 at the operating point | 44.1 µs | 369 µs |
| Best T2, frequency-constrained | 147 µs (EJ/EC ≥ 20; a grid-edge point at 1.45 GHz) | 719 µs (f01 in the published 0.17–0.83 GHz band; on the EL grid edge) |
| Anharmonicity \|α\| | 299 MHz | 8660 MHz (24× its own f01 — see the gate-count caveat) |
| Thermal occupation n̄_th at 15 mK | 1.4×10⁻⁷ | 0.45 |
| Bias window for a 2× T2 loss | 100% of the charge period (T2 varies 1.004×) | 2.2×10⁻⁴ of the flux period |
| Dominant noise channel | dielectric relaxation (critical current ~⅓ of the rate) | critical-current dephasing (76% of 1/T2 at Φ=0.5; 358 of 400 design points); first-order flux noise dominates beyond 8.7×10⁻⁵ Φ₀ |
| Sweet-spot falloff exponent | n/a (no sharp knob) | not a constant: −0.22 to −0.93, reaching −0.93 where 1/f theory predicts −1 |
| Practical character | simple, forgiving; no bias to hold | delicate; must hold Φ=0.5 to ~10⁻⁴ Φ₀, and is thermally populated |

**Q3 has no single answer.** The fluxonium/transmon ratio depends on the metric by a factor of
**~390×** across defensible metrics at the default noise amplitudes (242× against 0.62×) — and
the gate-count metric hides a further **393×** inside itself, in the choice of what bounds the
minimum pulse length.

Monte Carlo over literature noise ranges, N=400 per correlation mode (shared draws shown), 90% CI
on the median by bootstrap:

| Metric | Deterministic | MC median [5–95%] | Fluxonium wins |
|---|---|---|---|
| Raw T2, own operating points | 8.4× | 5.1× [3.6, 13.1] | 100% |
| **T2 at matched f01 (MC at 1 GHz)** | **0.62–0.82×** | **0.37× [0.25, 0.99]** | **5%** |
| Gate count, matched f01, 1/\|α\| bound | 364–647× | 262× [179, 698] | 100% |
| Gate count, own op. points, 1/\|α\| bound | 242× | 148× [105, 379] | 100% |
| Gate count, own op. points, 1/min(\|α\|,f01) | 10.2× | 6.2× [4.4, 15.9] | 100% |
| **Gate count, own op. points, 1/f01 bound** | **0.62×** | **0.38× [0.27, 0.97]** | **4%** |
| **Gate count, pulse simulation** | **22×** [17–27] | — | — |
| **Gate count, matched f01, 1/f01 bound** | 0.62–0.81× | **0.37× [0.26, 0.99]** | **5%** |

The raw coherence advantage is **a frequency effect**: under T2 ∝ 1/f01 the 13.6× frequency
ratio would more than account for the raw 8.4×, and at matched f01 the transmon is ahead at every
frequency tested and in 95% of noise draws. Matching frequency requires transmons at 0.5–2 GHz
with EJ/EC = 325–5050, outside any design envelope, so this is a control for the confound rather
than a device comparison. What survives for the fluxonium is converting coherence into operations
— but only to the extent that anharmonicity, rather than the qubit's own transition frequency,
sets how short a pulse can be. Its \|α\| = 8.7 GHz is 24× its 0.363 GHz transition frequency, so
the 1/\|α\| bound assumes a gate 24× shorter than one period of the transition being driven.

**`04_analysis.ipynb` §3.1 estimates it by simulation.** Integrating the lab-frame Schrödinger
equation with no RWA and numerically optimizing amplitude, DRAG and Stark detuning at each
duration gives a shortest tested error-<10⁻³ gate of 8 ns (transmon, bracketed 6–8) and 3 ns
(fluxonium, bracketed 2.5–3), so **N = T₂/T_gate gives 22×, bracketed 17–27×** — between the two
naive estimates *within this model*. The mechanism the α-only argument misses: the fluxonium's
drive couples 10.2× more strongly to its leakage transition than to the qubit transition
(transmon: 1.37×). Durations are upper bounds — the calibration is coarse and the optimizer is not
converged — so this is a model-dependent estimate, not an uncertainty interval. The price, on
every metric, is operating-point precision.

> **Superseded.** Earlier versions of this README reported "best T2 = 5578 µs", "bias sensitivity
> 528×", "T2 ∝ detuning^−0.68" and "Monte Carlo median 216×". All four were artefacts of method,
> not physics, and are corrected above; see `04_analysis.ipynb` Appendix A and
> `05_summary.ipynb` §4. The numbers of the 2026-10 version (62.4 µs, 1130 µs, 18.1×, 0.86×, 48×)
> are superseded by the critical-current convention change above.

## Repository layout

```
.
├── notebooks/     analysis notebooks, run in numbered order
│   ├── 01_validation.ipynb        reproduce published values (run this first)
│   ├── 02_transmon_sweep.ipynb    800-point EJ/EC sweep, dominant-channel map
│   ├── 03_fluxonium_sweep.ipynb   flux-bias sweep + EJ/EL design map
│   ├── 04_analysis.ipynb          Q3 comparison: three metrics, Monte Carlo robustness
│   └── 05_summary.ipynb           master table, provenance ledger, summary figure
├── build_nb0N.py   source of truth for each notebook — edit these, not the .ipynb
├── src/            reusable pipeline (qubits, coherence, plotting)
├── tests/          pytest suite: conventions + regression anchors
├── data/           generated CSVs + run metadata (JSON)
├── figures/        publication figures, PNG at dpi=300 + .pdf siblings
├── docs/           lab notebook and findings worksheet (git-ignored, local only)
└── paper/          LaTeX source: preprint.tex (arXiv) · main.tex (old IEEE draft, frozen)
```

**Notebooks are generated, not hand-edited.** Each `notebooks/0N_*.ipynb` is produced by the
corresponding `build_nb0N.py`; editing the `.ipynb` directly is undone by the next build. To
change a notebook:

```bash
python build_nb03.py && jupyter nbconvert --execute --inplace notebooks/03_fluxonium_sweep.ipynb
```

## Reproducing every figure

```bash
conda env create -f environment.yml     # one-time
conda activate qubits
jupyter lab
```

Check the pipeline first — `pytest tests/ -q` (35 tests, ~5 s) verifies the units, rate-addition
and channel conventions, pins the validated anchor values, and gates the published-number anchors
(Nguyen's second-order flux bound, Manucharyan's three parameter sets, the gate-count pulse bound). Then run the notebooks in order
`01 → 05`. Each notebook begins with a bootstrap cell that
sets the working directory to the project root, so `data/` and `figures/` resolve correctly
whether you launch from the repo root or from `notebooks/`.

Every sweep writes its raw results to `data/*.csv` **before** plotting, and saves a
`*_metadata.json` recording the parameters, noise amplitudes, and units for that run.

## Units and conventions

- **scqubits returns coherence times in nanoseconds.** All results here are converted to
  microseconds (`× 1e-3`). This is the single easiest mistake to make.
- **Rates add, not times:** `1/T2 = 1/(2·T1) + Σ 1/Tφ_k`. The shortest timescale dominates.
- A plain transmon has no flux loop, so `tphi_1_over_f_flux` is not a supported channel for
  it; the fluxonium supports it. Channels are chosen per architecture accordingly.
- **Critical-current noise is fractional** (A_cc = δI_c/I_c): `src/coherence.py` passes
  A_cc·E_J to scqubits, which otherwise applies A_cc as an absolute E_J fluctuation in GHz.

## Noise channels modelled

| Channel | scqubits method | Affects |
|---|---|---|
| Dielectric loss | `t1_capacitive` | T1 |
| Quasiparticle tunneling | `t1_quasiparticle_tunneling` | T1 (fluxonium) |
| 1/f flux noise, first order | `tphi_1_over_f_flux` | Tφ (fluxonium) |
| 1/f flux noise, second order | `tphi_flux_second_order` — **this project, not scqubits** | Tφ (fluxonium) |
| Charge (offset) noise | `tphi_1_over_f_ng` | Tφ (transmon) |
| Critical-current noise | `tphi_1_over_f_cc` | Tφ (both) |

scqubits implements 1/f dephasing to **first order only** and returns a divergent time wherever
∂ω₀₁/∂λ = 0 — documented behaviour, not a bug (Groszkowski & Koch 2021, §5.1 fn. 4). The
second-order channel is implemented in `src/coherence.py` from Nguyen et al. (2019) Eq. (A3).
It carries 0.2% of the rate budget at the nominal sweet spot, but 46–47% at the interior design
optimum of the extended (EJ, EL) scan (5–30 MHz) — an optimum that is invisible without it.

Because that channel is the one piece of physics written here rather than taken from scqubits, it
is validated against a published number: run on Nguyen's nine devices at *their* measured
A_Φ = 1.8 µΦ₀, it gives Tφ⁽²⁾ = 15–804 ms, clearing their stated bound of "> 10 ms for all qubits"
(Appendix A 2) on every device. See `01_validation.ipynb` §6.1.

## Limitations

These are model-based estimates using literature noise amplitudes, not measurements.

**Array modes.** The superinductance is ideal, so junction-array modes are absent. Explicit
diagonalization against Manucharyan's two modes (8.18 and 10.78 GHz) bounds the cost: f01 shifts
−16% at Φ_ext = 0 but only **−0.43%** at the flux sweet spot, where every coherence number here is
computed — smaller than the 9.5% spread across their own published parameter sets.

**Spectrum.** Against Manucharyan et al. (2009) the model reproduces f01 and f02 at Φ_ext = 0 to
better than 0.7% under *every* one of the paper's three published parameter sets (main text
9.0/2.5/0.52; Table 1 junction-model fit 8.93/2.39/0.52; Table 1 extended-model fit
8.97/2.47/0.52), and a single four-parameter fit reproduces all six digitized lines to an RMS of
0.0102 GHz. f03 at Φ_ext = 0 is excluded for a stated physical reason: it is strongly hybridized
with the junction-array mode at 10.78 GHz (g/|Δ| = 1.5–3.7), which is why Manucharyan needed their
extended model there — see `01_validation.ipynb` §9.2(iii). Note that
those three sets span **9.5%** in f01 at the flux sweet spot, so the frequently quoted "0.012%
agreement" is agreement with the *rounded main-text values*, not with the paper's fits — see
`01_validation.ipynb` §9.2.

**Calibration.** Run on the nine measured fluxonium devices of Nguyen et al. (2019), Table I
(published version), the model is optimistic by a median factor of **2.0× on T1 and 1.5× on T2**
(worst case 3.2×, device F; 2.6× and 5.7× under the old absolute critical-current reading). The
gap is not uniform: **device I**, the device nearest this project's nominal design in frequency
(0.395 vs 0.363 GHz), is the one on which the model is *pessimistic* on both, T1 (402 µs against
500 µs measured) and T2 (291 µs against 510 µs).

The gap is not uniform across devices, so it is not a single missing constant — the largest
identified contributors are the flat Q_cap = 10⁶ against Nguyen's
measured tan δ_C = (2.0–3.6)×10⁻⁶ at 6 GHz, the absence of Purcell decay, photon shot noise and
inductive loss, the fact that Nguyen's T2 is a Hahn echo while this model computes free-induction
decay, and a temperature mismatch (this model runs at scqubits' 15 mK; Nguyen state 20 mK, worth
~1.5× in thermal occupation at these frequencies).

**On the flat Q_cap, which is the project's largest single modelling choice.** Passing `Q_cap`
explicitly installs a frequency-*independent* quality factor in place of scqubits' default
Q_cap(ω) ∝ ω^−0.7, and at the fluxonium sweet spot the two differ by 7.1× in T1. That spread is
one-sided rather than symmetric: Nguyen fit their own T1-vs-frequency data with
tan δ_C ∝ ω^0.15, i.e. Q_cap ∝ ω^−0.15. At f01 = 0.363 GHz the flat model sits 1.52× from their
measured exponent while the scqubits default sits 4.7× from it, so the calibration paper puts the
truth much nearer the flat end. The flat value is retained as both the conservative and the
better-supported choice — but it is still a modelling choice, and Nguyen fit that exponent jointly
with thermal stimulated emission, so it is not a free-standing measurement of the loss tangent.

**Both architectures' raw T2 maxima are grid-edge artefacts.** The transmon's 147 µs sits at
EJ=3.03, EC=0.10 — the smallest EC on the grid and the smallest *sampled* EJ with EJ/EC ≥ 20 — at
f01 = 1.45 GHz, below the conventional 3–6 GHz band. The fluxonium's unconstrained maximum, 1398 µs,
sits at EJ=15, EL=0.2 — both grid edges — at 75 MHz with n̄_th = 3.7. In both planes log T2
correlates with log f01 at |r| > 0.98 (−0.9999 and −0.9878), so a raw T2 maximum
is a frequency minimum in disguise. Neither is a design recommendation; quote the
frequency-constrained values together with their constraint.

The transmon's *global* maximum T2 occurs at EJ/EC ≈ 0.5, which is the charge-qubit regime rather
than a physical transmon — headline transmon numbers are therefore reported for EJ/EC ≥ 20.

Absolute values carry substantial uncertainty; comparative trends are the defensible output.

## Citation

This work uses scqubits:

> P. Groszkowski and J. Koch, *scqubits: a Python package for superconducting qubits*,
> Quantum 5, 583 (2021). DOI: 10.22331/q-2021-11-17-583

The physics and reference data come from:

> J. Koch *et al.*, *Charge-insensitive qubit design derived from the Cooper pair box*,
> Phys. Rev. A **76**, 042319 (2007). DOI: 10.1103/PhysRevA.76.042319
>
> P. Krantz *et al.*, *A Quantum Engineer's Guide to Superconducting Qubits*,
> Appl. Phys. Rev. **6**, 021318 (2019). DOI: 10.1063/1.5089550
>
> V. E. Manucharyan *et al.*, *Fluxonium: Single Cooper-Pair Circuit Free of Charge Offsets*,
> Science **326**, 113 (2009). DOI: 10.1126/science.1175552
>
> L. B. Nguyen *et al.*, *High-Coherence Fluxonium Qubit*, Phys. Rev. X **9**, 041041 (2019).
> DOI: 10.1103/PhysRevX.9.041041
>
> G. Ithier *et al.*, *Decoherence in a superconducting quantum bit circuit*,
> Phys. Rev. B **72**, 134519 (2005). DOI: 10.1103/PhysRevB.72.134519

## License

MIT — see [LICENSE](LICENSE).
