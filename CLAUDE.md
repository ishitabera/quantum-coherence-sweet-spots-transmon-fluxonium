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
| Second-order 1/f flux (Tφ) | `tphi_flux_second_order` — **our code, not scqubits** | **no** | yes |
| Charge / offset (Tφ) | `tphi_1_over_f_ng` | yes | no |
| Critical current (Tφ) | `tphi_1_over_f_cc` | yes | yes |

The second-order channel is implemented locally in `src/coherence.py` because scqubits
computes 1/f dephasing to **first order only**, so `tphi_1_over_f_flux` diverges to `inf`
exactly at the flux sweet spot (Φ = Φ₀/2) where ∂ω/∂Φ = 0 — which would make the sweet spot
look infinitely coherent. The second-order term, driven by ∂²ω/∂Φ², is what actually limits
the fluxonium there; see Nguyen et al., PRX **9**, 041041 (2019), Eq. (A3), and
Groszkowski & Koch, Quantum **5**, 583 (2021), Sec. 5.1 footnote 4. It applies to the
fluxonium only — a plain transmon has no flux loop, per the rule above.

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

**A raw T2 maximum is a frequency minimum in disguise — in BOTH architectures.** corr(log T2,
log f01) is −0.9999 over the transmon plane (EJ/EC ≥ 20) and −0.9973 over the fluxonium (EJ, EL)
plane, because dielectric loss dominates both and scales with frequency. So the "best T2" point of
any sweep is wherever the frequency constraint was drawn. Both of this project's former headline
maxima were double-grid-edge, out-of-band artefacts: the transmon's 187 µs at EJ=3.0, EC=0.10
(f01 = 1.45 GHz, below the 3–6 GHz band) and the fluxonium's 5427 µs at EJ=15, EL=0.2
(f01 = 75 MHz, n̄_th = 3.7; long quoted as 5578 µs, which was its first-order-only value). **Quote a constrained optimum together with its constraint, never a
raw maximum**, and check whether the optimum sits on a grid edge before reporting it.

**Comparing architectures requires controlling for frequency.** A raw cross-architecture T2 ratio
partly measures the frequency difference, not the architecture. The fluxonium's 18× raw advantage
falls to a median 0.86× at matched f01 across the noise band (`04_analysis.ipynb` §3–4). Report the
metric alongside the ratio; the answer spans ~700× across defensible metrics.

**The gate-count metric needs its pulse bound quoted with it.** `N = T2 / T_gate` is this project's
heuristic, and `T_gate` is set by whatever bounds the pulse bandwidth. For a transmon
|α| = 299 MHz ≪ f01 = 4.9 GHz and every bound agrees. For the nominal fluxonium |α| = 8.7 GHz is
**24× its own** f01 = 0.363 GHz, so the bound decides the answer: the fluxonium/transmon ratio is
**524× under 1/|α|, 22.0× under 1/min(|α|, f01), and 1.34× under 1/f01** — and only the first two
survive the noise band (the drive-period bound gives median 0.83×, 40% win rate). Quoting 524×
alone assumes a gate 24× shorter than one period of the transition being driven.

**`04_analysis.ipynb` §3.1 settles it: 48×, bracketed 36–58×**, from a lab-frame no-RWA pulse
simulation with amplitude, DRAG and Stark detuning calibrated at each duration and leakage counted
in the fidelity. Both naive bounds are excluded. The mechanism they miss is the drive matrix
element: the fluxonium's |⟨1|n̂|2⟩|/|⟨0|n̂|1⟩| = 10.2 against the transmon's 1.37, so the drive is
ten times better at causing leakage than at driving the qubit, which eats most of the
anharmonicity advantage. **The simulated durations are UPPER bounds** — the calibration is coarse
(3 parameters, fixed Gaussian) and the optimiser is not converged at short durations, which shows
as non-monotonicity in the error curve; a monotone envelope is used to read off T_gate and the
bracket is reported. Optimal control would shorten both architectures, not necessarily equally.
All four numbers are in the ledger. Never quote one without its bound.

**The transmon design envelope is Krantz 2019 §II A**: f01 ∈ [3, 6] GHz, |α| ∈ [100, 300] MHz,
EJ/EC ≥ 50. Check any transmon operating point against all three axes before using it — the
superseded Week-4 comparison used an 18.9 GHz qubit with |α| = 1149 MHz and was penalised 3.7× on
T2 for that reason alone.

## Layout and how code is organised

```
src/           qubits.py (constructors) · coherence.py (physics) · plotting.py (figures)
notebooks/     01→05, run in numbered order (03–05 generated by build_nb0N.py)
data/          generated CSVs + *_metadata.json  (tracked in git — they are the results)
figures/       publication PNGs at dpi=300 + .pdf siblings  (tracked in git)
docs/          lab_notebook.md, Week5_Findings_Worksheet.md  (gitignored - local only)
paper/         LaTeX source (IEEE two-column)
results/       scratch only; results/scratch/ is gitignored
```

**`src/` uses flat imports** (`from qubits import make_transmon`), so `src/` must be on
`sys.path` — it is not a package. The bootstrap cell at the top of every notebook handles
this and also chdir's to the project root so `data/` and `figures/` resolve identically
whether Jupyter was launched from the repo root or from `notebooks/`. Copy that cell verbatim
into any new notebook; do not invent a different path scheme.

**No notebook duplicates `src/` logic any more.** Notebooks 02–05 all import from
`src/coherence.py`; the inline `transmon_T2` / `fluxonium_T2` copies in 04 and 05 were removed
on 2026-09-07. They had already diverged: both predated the second-order flux channel, so every
fluxonium number in the old 04/05 came from a superseded pipeline. **Do not reintroduce an inline
copy for convenience** — that divergence is exactly what cost this project its Week-4 headline.

**Notebooks are regenerated from `build_nb0N.py`, not hand-edited.** `build_nb01.py` through
`build_nb05.py` hold the source-of-truth cell definitions for all five notebooks. To change a
notebook, edit the builder and re-run:

```bash
python build_nb04.py && jupyter nbconvert --execute --inplace notebooks/04_analysis.ipynb
```

**Editing an `.ipynb` directly is silently undone by the next build.** The builders for 01 and 02
were reverse-engineered from the executed notebooks and verified to round-trip their cell sources
exactly; if you regenerate them, re-run that check.

**Tests.** `pytest tests/` (30 tests, ~3 s) locks the conventions in this file and the anchor
values every committed CSV depends on:

```bash
pytest tests/ -q
```

`tests/test_conventions.py` covers the rules above — units, rates-add, `T2 <= 2*T1`, the
per-architecture channel lists, `safe()` failure logging, and the second-order flux formula against
Nguyen Eq. (A3). `tests/test_regression.py` pins the validated anchors (54.352 µs transmon,
1136.869 µs / 1129.712 µs fluxonium), the scqubits `NOISE_PARAMS` defaults, and gates the
provenance ledger and the superseded-number audit. **If a regression test fails, every artifact in
`data/` and `figures/` is invalidated** — regenerate them deliberately rather than relaxing the
test.

**Every headline number must survive the ledger.** `05_summary.ipynb` §2 recomputes **32**
quoted values from their named file and filter and asserts each one, so a number that drifts from
its artefact fails a cell rather than reaching the paper. If you add a headline number anywhere,
add it to that ledger; `tests/test_regression.py` fails if the ledger shrinks below 28 rows.

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

**`Q_cap=1e6` is a deliberate modeling choice, not a default.** Passing `Q_cap` explicitly
installs a **flat, frequency-independent** quality factor, overriding scqubits' default
frequency-dependent `Q_cap(ω) = 1e6 * (2π·6 GHz / ω)^0.7` (Smith et al. 2020). Because the
fluxonium's sweet-spot 0→1 frequency is far below 6 GHz, that exponent makes a large
difference to `t1_capacitive` there — the flat value is substantially more pessimistic.
Keep it flat for cross-architecture comparability, but never present the resulting absolute
fluxonium T1 as a literature-calibrated prediction.

That 7.1× is **not** a symmetric ambiguity, and `03_fluxonium_sweep.ipynb` §8 now says so: Nguyen
2019 §III A fit their own T1-vs-frequency data with tan δ_C ∝ ω^0.15, i.e. `Q_cap ∝ ω^-0.15`. At
the fluxonium sweet spot the flat model sits 1.52× from their measured exponent while the scqubits
default sits 4.7× from it. The flat choice is both the conservative *and* the better-supported
one — but Nguyen fit that exponent jointly with thermal stimulated emission, so it is not a
free-standing measurement of the loss tangent.

**Manucharyan 2009 reports three parameter sets, not one.** Main text p. 3 gives EJ=9.0, EC=2.5,
EL=0.52 to **two** significant figures; Table 1 gives two least-squares fits to three —
inductively-shunted **junction** model (8.93, 2.39, 0.52) and **extended** model (8.97, 2.47,
0.52). They span **9.5%** in f01 at the flux sweet spot, because f01 there is exponentially
sensitive to EC. The often-quoted "spectrum validated to 0.012%" is agreement with the *rounded
main-text values*, not with the paper's fits. `01_validation.ipynb` §5 runs the spectrum comparison
at all three: f01 and f02 at Φ_ext=0 are reproduced to <0.7% under every set, and a single
four-parameter fit reproduces all six digitized lines (both flux points) to an RMS of 0.0102 GHz.

**The f03 row at Φ_ext=0 is excluded for a physical reason, not a loose tolerance.** Its residual
is not a misassignment (the measured line disperses +0.046 GHz over Φ 0→0.05, the same sign as
f03 and the *opposite* sign to the two-photon 0→4) and not merely E_C sensitivity (though
∂f03/∂E_C = +1.71 GHz/GHz against −0.19 and −0.10 for f01, f02 — which is why f03 is the sensitive
line). It is **strong hybridization with the junction-array mode** Manucharyan's extended model
puts at ω₂/2π = 10.78 GHz: g₀₃ ≈ 0.53 GHz against a detuning of 0.14–0.36 GHz, i.e. g/|Δ| = 1.5–3.7.
The bare Hamiltonian cannot reproduce that line, which is why the paper needed the extended model.
Do not re-add f03 to the spectrum-validation claim.

**The model has no junction-array modes, and the cost is bounded at sub-percent.** The
superinductance is ideal. `01_validation.ipynb` §5 step (iv) diagonalises the fluxonium coupled to
both of Manucharyan's array modes: f01 shifts −16% at Φ_ext=0 but only **−0.43%** at the flux
sweet spot, where every coherence number here is evaluated. Do not quote the naive g²/Δ estimate at
the sweet spot — it comes out the same size as the transition frequency and is outside its own
regime of validity.

## Writing about results

Claims in `docs/` and `paper/` must be traceable to a CSV in `data/`. When stating a headline
number, name the file and the filter that produced it (e.g. "187 µs, `transmon_full_sweep.csv`,
EJ/EC ≥ 20"). Comparative conclusions are more robust than absolute values — but *"comparative
conclusions are robust"* is itself too strong as a blanket claim: the matched-frequency T2
comparison flips sign under the noise-amplitude uncertainty. Say which comparison, at which metric.

**Absolute coherence times carry a known calibration.** Against the **nine** devices of Nguyen
2019's *published* Table I the model is optimistic by a median 2.0× on T1 and 2.6× on T2 (worst
case 5.7×). Quote that whenever an absolute number appears. Note also that the model runs at
scqubits' T = 15 mK while Nguyen state 20 mK — worth ~1.5× in thermal occupation at these
frequencies, and part of the gap. The gap is **not uniform**: device I,
the device nearest this project's nominal design in frequency (0.395 vs 0.363 GHz), is matched to
4% on T2 with the model *pessimistic* on its T1 — the optimism is concentrated in devices further
from that regime.

**Cite the published Nguyen paper, not the preprint.** Table I in *Phys. Rev. X* **9**, 041041
(2019) has **nine** devices (A–I); arXiv:1810.11006v1 has **eight** and omits device I entirely —
the paper's headline device and the one closest to this project's design. That omission silently
skewed the calibration for weeks and no internal check could catch it, because the transcription was
self-consistent and faithfully reproduced a real table. It was the wrong table. Equation numbering
also differs: published (A1)/(A3) versus preprint (A.1)/(A.3).

**`docs/` is gitignored and local-only.** It is not in a clean clone, so
`05_summary.ipynb` §4's superseded-number audit covers `README.md` alone there. The audit reports
its own coverage (`superseded_audit.docs_dir_present`); a pass with that flag false says nothing
about `docs/`. Run notebook 05 on a checkout that has `docs/` before trusting the audit, and keep
`paper/` claims traceable to the §2 ledger rather than to `docs/` prose.

**`docs/lab_notebook.md` is append-only.** Record corrections in a dated entry with a corrections
table; never retro-edit an earlier entry — supersede it with a new dated row instead. The
interpretation sections were originally reserved for the author and left blank; they were written
with Claude on 2026-09-07 at the author's request and each carries a provenance line saying so.
Keep that line: it is the difference between the author's analysis and a drafted one.
