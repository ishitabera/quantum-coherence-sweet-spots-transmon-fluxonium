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

| Metric | Transmon | Fluxonium |
|---|---|---|
| Best T2 (real operating regime) | 187 µs | 5578 µs (5.6 ms) |
| Bias-knob sensitivity (T2 max/min) | 1.01× | 528× |
| Dominant noise channel | dielectric relaxation | flux noise off spot; relaxation at Φ=0.5 |
| Sweet-spot scaling | flat (no sharp knob) | T2 ∝ (detuning)^−0.68 |
| Practical character | simple, forgiving | delicate, must hold Φ=0.5 |

Monte Carlo over literature noise ranges (N=60): the fluxonium exceeds the transmon in
**100%** of samples, median ratio **216×**, 90% CI **[52×, 371×]**.

## Repository layout

```
.
├── notebooks/     analysis notebooks, run in numbered order
│   ├── 01_validation.ipynb        reproduce published values (run this first)
│   ├── 02_transmon_sweep.ipynb    800-point EJ/EC sweep, dominant-channel map
│   ├── 03_fluxonium_sweep.ipynb   flux-bias sweep + EJ/EL design map
│   ├── 04_analysis.ipynb          architecture comparison, Monte Carlo robustness
│   └── 05_summary.ipynb           master table + 4-panel summary figure
├── src/           reusable pipeline (qubits, coherence, plotting)
├── data/          generated CSVs + run metadata (JSON)
├── figures/       publication figures (PNG)
├── docs/          lab notebook and findings worksheet
└── paper/         LaTeX source (IEEE two-column)
```

## Reproducing every figure

```bash
conda env create -f environment.yml     # one-time
conda activate qubits
jupyter lab
```

Then run the notebooks in order `01 → 05`. Each notebook begins with a bootstrap cell that
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
| 1/f flux noise | `tphi_1_over_f_flux` | Tφ (fluxonium) |
| Charge (offset) noise | `tphi_1_over_f_ng` | Tφ (transmon) |
| Critical-current noise | `tphi_1_over_f_cc` | Tφ (both) |

## Limitations

These are model-based estimates using literature noise amplitudes, not measurements.
Absolute values carry substantial uncertainty; the comparative conclusions are robust
(see the Monte Carlo result above). The transmon's global maximum T2 occurs at EJ/EC ≈ 0.5,
which is the charge-qubit regime rather than a physical transmon — headline transmon numbers
are therefore reported for EJ/EC ≥ 20.

## Citation

This work uses scqubits:

> P. Groszkowski and J. Koch, *scqubits: a Python package for superconducting qubits*,
> Quantum 5, 583 (2021).

## License

MIT — see [LICENSE](LICENSE).
