# Paper source

IEEE two-column draft of the study. `main.tex` + `refs.bib`.

## Build

```bash
brew install tectonic          # one-time; self-contained, fetches packages on demand
cd paper && tectonic -X compile main.tex
```

Compiles clean to **7 pages**: no overfull boxes, all 14 citations resolved, all cross-
references resolved. Four underfull boxes remain (loose interword spacing in justified
two-column text) — cosmetic. `latexmk -pdf main.tex` also works if you have MacTeX.

`main.pdf` is committed. Build intermediates are gitignored. Figures are copied into
`paper/figures/` so the directory is self-contained for submission; regenerate them with
`cp ../figures/*.pdf figures/` after re-running the notebooks.

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

- **Novelty is not established.** §II states the contribution against the literature I did
  survey, but that survey is not exhaustive and I say so in the text. The frequency scaling
  itself is standard; the claim is about its consequence for design-space studies. Before
  submitting, do a proper search — if someone has already published the metric-sensitivity
  point, the framing needs to change.
- **Venue not chosen.** `IEEEtran` conference class was inherited from the project
  convention, not from a decision. At 7 pages this exceeds a typical 6-page conference
  limit; a journal (PRApplied, IEEE TQE) is probably the better fit and would want a
  different class file.
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
