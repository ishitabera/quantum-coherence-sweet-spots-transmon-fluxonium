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

**Submission-ready except for the author voice pass.**

Formatted for **IEEE Transactions on Quantum Engineering** (IEEEtran, journal mode) —
confirmed by the author. 8 pages, 16 references, 5 vector figures, 3 tables.

Compile and compliance, all verified against the built PDF:

| | |
|---|---|
| overfull boxes | 0 |
| underfull boxes | 5 (loose interword spacing; cosmetic) |
| unresolved citations / cross-refs | 0 |
| abstract | 245 words, single paragraph (IEEE norm 150–250) |
| Index Terms | present |
| affiliation | independent researcher, self-directed project |
| AI-use disclosure | present, cites IEEE policy |
| Reproducibility and Data Availability | present, merged (was duplicated) |
| headline numbers matching the ledger | 28 / 28 |

### The one thing left

**Author voice pass.** The prose was drafted from the committed artifacts, not written by
the author. The physics and every number are checked; the sentences are not yours yet. Read
it line by line before submitting — for a single-author paper this is an authorship matter,
not a style preference.

### Limits stated in the paper itself, not hidden here

- **Novelty is argued, not proven.** §II-A locates the two closest works — McRae *et al.* on
  experimental reproducibility of coherence comparisons, and Lange *et al.* on controlled
  transmon/fluxonium comparison for cross-talk rather than coherence — and §II-B claims only
  "we are not aware of", stating the survey's incompleteness in the body text.
- **The 48× gate-count result rests on a coarse calibration.** Three parameters over a fixed
  Gaussian; durations are upper bounds; the optimiser is not converged at the shortest
  durations. §XI names GRAPE/CRAB as the next step.
- **Frequency matching is necessary but not sufficient** as a control: geometry and surface
  participation also set dielectric loss, and a Hamiltonian-level model does not capture
  them. Stated in §II-C as a caveat on our own method.
- **Model-based throughout.** No new measurement; the calibration set is nine devices from
  one 2019 paper, and §IX notes the field has since moved past them.
