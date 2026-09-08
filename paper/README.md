# Paper source

IEEE two-column draft of the study. `main.tex` + `refs.bib`.

## Build

No LaTeX toolchain is installed in the `qubits` environment, so this has **not been
compiled**. To build:

```bash
# macOS: brew install --cask mactex-no-gui   (or basictex + tlmgr install ieeetran siunitx)
cd paper && latexmk -pdf main.tex
```

`IEEEtran.cls` ships with TeX Live/MacTeX. `siunitx`, `booktabs` and `hyperref` are also
required.

## Rules for editing this paper

**Every number here must exist in `data/headline_provenance_ledger.csv` first.**
`05_summary.ipynb` §2 recomputes each ledger entry from its named file and filter and
asserts it, and `tests/test_regression.py` fails if the ledger shrinks or any row stops
verifying. The workflow is:

1. add the number to the ledger in `build_nb05.py`
2. re-run `python build_nb05.py && jupyter nbconvert --execute --inplace notebooks/05_summary.ipynb`
3. confirm it verifies
4. only then quote it here

Do not copy numbers out of notebook prose — prose drifts, the ledger does not. The
superseded-number audit in `05_summary.ipynb` §4 scans `README.md` and `docs/` but **not**
`paper/`; if a headline number is ever retired, grep this directory by hand.

## Status

Complete draft, not submitted. Known gaps:

- Figures are included by relative path from `../figures/` (the `.pdf` siblings that
  `src/plotting.save` writes alongside every PNG): `summary_figure`, `transmon_tradeoff`,
  `fluxonium_interior_optimum`, `q3_robustness`, `q3_pulse_bound`. If you move this file,
  fix the paths or copy the PDFs in.
- **Voice pass.** Sections drafted with Claude from the committed artifacts, at the
  author's request, mirroring the notebook interpretation sections. The physics and the
  numbers are checked; the prose has not had an author pass.
- **Author list, affiliation and acknowledgements** are placeholders.
- The pulse-simulation section reports a coarse three-parameter calibration over a fixed
  Gaussian shape, so the gate durations it reports are **upper bounds** and the resulting
  gate-count ratio is an estimate, not a settled value. The optimiser is not reliably
  converged at the shortest durations. A GRAPE/CRAB optimal-control treatment over a free
  pulse shape is the obvious next step, and would shorten both architectures' gates — not
  necessarily equally, which is precisely why it matters.
