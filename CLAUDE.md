# Coherence Sweet Spots — Transmon vs Fluxonium

Parameter-sensitivity study of transmon and fluxonium coherence under realistic noise
channels, built on [scqubits](https://scqubits.readthedocs.io). All numbers are
**model-based estimates**, not measurements. The value is in relative trends across
parameters and architectures — never present absolute coherence times as predictions.

Research questions: **Q1** where are the sweet spots across (EJ, EC) / (EJ, EC, EL, flux);
**Q2** which noise channel dominates where; **Q3** does the fluxonium's flux sweet spot beat
the transmon's charge insensitivity, and at what design cost.

## Non-negotiable conventions

These are the mistakes that have actually cost time in this project.

**Units.** scqubits returns coherence times in **nanoseconds**. Everything in `data/`,
`figures/`, `docs/` and the paper is in **microseconds**. Every conversion is `* 1e-3`
(`NS_TO_US` in `src/coherence.py`). A result that is off by 1000× is a units bug, not physics.

**Rates add, not times.** `1/T2 = 1/(2*T1) + Σ_k 1/Tφ_k`. Never average or min() coherence
times to combine channels. The shortest timescale dominates by construction.

**A plain transmon has no flux loop.** `tphi_1_over_f_flux` is not a valid channel for
`scq.Transmon` and its absence is intentional, not an oversight. Channels are chosen per
architecture:

| Channel | scqubits method | Transmon | Fluxonium |
|---|---|---|---|
| Dielectric loss (T1) | `t1_capacitive` | yes | yes |
| Quasiparticle tunneling (T1) | `t1_quasiparticle_tunneling` | no | yes |
| 1/f flux (Tφ) | `tphi_1_over_f_flux` | **no** | yes |
| Charge / offset (Tφ) | `tphi_1_over_f_ng` | yes | no |
| Critical current (Tφ) | `tphi_1_over_f_cc` | yes | yes |

**Save raw data before plotting.** Every sweep writes `data/*.csv` *before* any figure code
runs, so a plotting error never costs a multi-minute sweep. Each sweep also writes a
`data/*_metadata.json` recording date, grid definition, fixed params, noise amplitudes,
channel list, units string, and output file paths.

**`safe()` returns `inf` on failure.** In `src/coherence.py`, an unsupported or failing
channel returns `np.inf`, which contributes zero rate to the T2 sum. This is deliberate —
but it means a silently broken channel looks like a *good* result. When a T2 jumps
unexpectedly, check whether a channel started returning `inf`.

**EJ/EC ≥ 20 is the real transmon regime.** The global T2 maximum over the swept plane sits
near EJ/EC ≈ 0.5, which is a charge qubit, not a transmon. Headline transmon numbers are
always reported for the EJ/EC ≥ 20 subset. Do not quote the global max as a transmon result.

## Layout and how code is organised

```
src/           qubits.py (constructors) · coherence.py (physics) · plotting.py (figures)
notebooks/     01→05, run in numbered order
data/          generated CSVs + *_metadata.json  (tracked in git — they are the results)
figures/       publication PNGs at dpi=150       (tracked in git)
docs/          lab_notebook.md, Week5_Findings_Worksheet.md
paper/         LaTeX source (IEEE two-column)
results/       scratch only; results/scratch/ is gitignored
```

**`src/` uses flat imports** (`from qubits import make_transmon`), so `src/` must be on
`sys.path` — it is not a package. The bootstrap cell at the top of every notebook handles
this and also chdir's to the project root so `data/` and `figures/` resolve identically
whether Jupyter was launched from the repo root or from `notebooks/`. Copy that cell verbatim
into any new notebook; do not invent a different path scheme.

**Notebooks currently duplicate `src/` logic.** `02_transmon_sweep.ipynb` defines its own
`transmon_coherence` inline rather than importing from `src/coherence.py`. The two are
intended to agree. If you change the physics in one, change it in the other, or the notebook
and the reusable pipeline will silently diverge — the notebook is what produced the committed
CSVs, so treat it as the source of truth for existing results.

## Environment

Conda env `qubits` (Python 3.11): `conda env create -f environment.yml && conda activate qubits`.
`scq.settings.T1_DEFAULT_WARNING = False` is set in `src/qubits.py` to silence a harmless
warning — keep it.

## Nominal designs

```python
TRANSMON_DEFAULTS  = dict(EJ=15.0, EC=0.3, ng=0.0,  ncut=30)     # ng=0   = charge sweet spot
FLUXONIUM_DEFAULTS = dict(EJ=8.9, EC=2.5, EL=0.5, flux=0.5, cutoff=110)  # flux=0.5 = flux sweet spot
DEFAULTS = dict(Q_cap=1e6, A_ng=1e-4, A_flux=1e-6, A_cc=1e-7)    # literature noise amplitudes
```

Changing any of these invalidates every committed CSV and figure. If you change one, say so
explicitly and note which artifacts need regenerating.

## Writing about results

Claims in `docs/` and `paper/` must be traceable to a CSV in `data/`. When stating a headline
number, name the file and the filter that produced it (e.g. "187 µs, `transmon_full_sweep.csv`,
EJ/EC ≥ 20"). Comparative conclusions are robust; absolute values are not — hedge accordingly.
