# Coherence Sweet Spots Across Superconducting Qubit Architectures

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

| Metric | Transmon | Fluxonium |
|---|---|---|
| Operating point | EJ=13.0, EC=0.26 (EJ/EC=50), f01=4.92 GHz | EJ=8.9, EC=2.5, EL=0.5, f01=0.363 GHz |
| T2 at the operating point | 62.4 µs | 1130 µs |
| Best T2, frequency-constrained | 187 µs (EJ/EC ≥ 20) | 2412 µs (f01 in the published 0.17–0.83 GHz band) |
| Anharmonicity \|α\| | 299 MHz | 8660 MHz (24× its own f01 — see the gate-count caveat) |
| Thermal occupation n̄_th at 15 mK | 1.4×10⁻⁷ | 0.45 |
| Bias window for a 2× T2 loss | 100% of the charge period (T2 varies 1.006×) | 7.2×10⁻⁵ of the flux period |
| Dominant noise channel | dielectric relaxation | relaxation at Φ=0.5; first-order flux noise dominates the rate budget beyond 2.75×10⁻⁵ Φ₀ |
| Sweet-spot falloff exponent | n/a (no sharp knob) | not a constant: −0.45 to −0.96, reaching −0.96 where 1/f theory predicts −1 |
| Practical character | simple, forgiving; no bias to hold | delicate; must hold Φ=0.5 to ~10⁻⁴ Φ₀, and is thermally populated |

**Q3 has no single answer.** The fluxonium/transmon ratio depends on the metric by a factor of
**491×** at the default noise amplitudes (~700× across the Monte Carlo medians) — and the
gate-count metric hides a further **393×** inside itself, in the choice of what bounds the minimum
pulse length.

Monte Carlo over literature noise ranges, N=400 per correlation mode, 90% CI on the median by
bootstrap:

| Metric | Deterministic | MC median [5–95%] | Fluxonium wins |
|---|---|---|---|
| Raw T2, own operating points | 18.1× | 11.3× [3.6, 21.1] | 100% |
| **T2 at matched f01** | 1.07–2.16× | **0.86× [0.27, 1.58]** | **42%** |
| Gate count, matched f01, 1/\|α\| bound | 625–1709× | 606× [191, 1114] | 100% |
| Gate count, own op. points, 1/\|α\| bound | 524× | 326× [104, 611] | 100% |
| Gate count, own op. points, 1/min(\|α\|,f01) | 22.0× | 13.7× [4.4, 25.6] | 100% |
| **Gate count, own op. points, 1/f01 bound** | **1.34×** | **0.83× [0.26, 1.56]** | **40%** |
| **Gate count, matched f01, 1/f01 bound** | 1.07–2.14× | **0.86× [0.27, 1.58]** | **42%** |

The raw coherence advantage is **substantially a frequency effect**: control for f01 and the T2
advantage disappears — across the noise band the transmon is the better device more often than
not. What survives is **anharmonicity per unit coherence** (the gate-count rows) — but only to the
extent that anharmonicity, rather than the qubit's own transition frequency, sets how short a
pulse can be. The fluxonium's \|α\| = 8.7 GHz is 24× its 0.363 GHz transition frequency, so the
1/\|α\| bound assumes a gate 24× shorter than one period of the transition being driven. Bound it
instead by the drive period and the advantage falls to 1.34×, and does not survive the noise band
(40% win rate). **This model does not settle the bound**; `04_analysis.ipynb` §3 computes all
three and the paper must quote the bound alongside the number. The price, on every metric, is
operating-point precision.

> **Superseded.** Earlier versions of this README reported "best T2 = 5578 µs", "bias sensitivity
> 528×", "T2 ∝ detuning^−0.68" and "Monte Carlo median 216×". All four were artefacts of method,
> not physics, and are corrected above; see `04_analysis.ipynb` Appendix A and
> `05_summary.ipynb` §4. In particular 5578 µs was a first-order-only value at a grid-edge design
> at 75 MHz with n̄_th = 3.7, and the 216× comparison used an 18.9 GHz "transmon" held at its
> worst-case bias against that same fluxonium.

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
└── paper/          LaTeX source (IEEE two-column)
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

Check the pipeline first — `pytest tests/ -q` (30 tests, ~3 s) verifies the units, rate-addition
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
It carries 0.63% of the rate budget at the nominal sweet spot, but 43–58% at the interior design
optimum of the (EJ, EL) plane — an optimum that is invisible without it.

Because that channel is the one piece of physics written here rather than taken from scqubits, it
is validated against a published number: run on Nguyen's nine devices at *their* measured
A_Φ = 1.8 µΦ₀, it gives Tφ⁽²⁾ = 15–804 ms, clearing their stated bound of "> 10 ms for all qubits"
(Appendix A 2) on every device. See `01_validation.ipynb` §6.1.

## Limitations

These are model-based estimates using literature noise amplitudes, not measurements.

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
(published version), the model is optimistic by a median factor of **2.0× on T1 and 2.6× on T2**
(worst case 5.7×). The gap is smallest where it matters most: **device I**, the device nearest this
project's nominal design in frequency (0.395 vs 0.363 GHz), is matched to **4% on T2**, and the
model is *pessimistic* on its T1 (402 µs against 500 µs measured).

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

**Both architectures' raw T2 maxima are grid-edge artefacts.** The transmon's 187 µs sits at
EJ=3.0, EC=0.10 — both grid edges — at f01 = 1.45 GHz, outside the 3–6 GHz band real transmons
occupy. The fluxonium's unconstrained maximum sits at EJ=15, EL=0.2 — also both grid edges — at
75 MHz with n̄_th = 3.7. In both planes T2 correlates with f01 at |r| > 0.99, so a raw T2 maximum
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
