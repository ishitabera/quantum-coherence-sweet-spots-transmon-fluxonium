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

**Compiles clean. Content-complete. Not submittable until the three items below are closed.**

8 pages, IEEEtran journal mode, 16 references, 5 figures, 3 tables. Zero overfull boxes;
six underfull (loose interword spacing in justified two-column text), one severe. All
citations and cross-references resolve. Every headline number appears in the PDF and
matches `data/headline_provenance_ledger.csv`.

### Must fix before submission

1. **Affiliation is a placeholder.** `\thanks{Affiliation to be completed}` — I left it
   unstated rather than guess your institution. Fill it in.
2. **Venue is a default, not a decision.** Formatted for IEEE Transactions on Quantum
   Engineering (IEEEtran journal mode) because the paper runs past the 6-page limit typical
   of IEEE conferences. Phys. Rev. Applied also suits the content but needs conversion to
   REVTeX 4.2 — a real rewrite of the preamble, floats and bibliography style, not a
   one-line change.
3. **Author voice pass.** The prose is drafted, not written by you. The physics and the
   numbers are checked; the sentences are not yours.

### Known limits of the science, stated in the paper itself

- **Novelty is argued, not proven.** §II-A locates the two closest works (McRae *et al.* on
  experimental reproducibility of coherence comparisons; Lange *et al.* on controlled
  transmon/fluxonium comparison for cross-talk rather than coherence) and §II-B claims only
  "we are not aware of", with the survey's incompleteness stated in the text. A referee may
  know a paper I did not find.
- **The 48× gate-count result rests on a coarse calibration.** Three parameters over a fixed
  Gaussian; durations are upper bounds; the optimiser is not converged at the shortest
  durations. GRAPE/CRAB is the obvious next step and §XI says so.
- **Frequency matching is necessary but not sufficient** as a control — geometry and surface
  participation also set dielectric loss and a Hamiltonian-level model does not capture them.
  Stated as a caveat on our own method in §II-C.

### AI-use disclosure

The Acknowledgements carry a disclosure that a large language model assisted with the
notebook interpretation sections and the manuscript draft, working from the committed
numerical artifacts. Check the wording against your target venue's specific policy — IEEE
and APS both now require disclosure, but the required form differs.
