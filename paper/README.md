# Paper source

**This paper is going to arXiv only** (the author's decision, 2026-10-01). There is no journal
submission.

| file | role |
|---|---|
| `preprint.tex` → `preprint.pdf` | **the live source.** Single column, 16 pages. |
| `main.tex` → `main.pdf` | the earlier IEEE TQE version, **frozen**. Not maintained, and not the source of truth. |

Both use `refs.bib` and `figures/`. Edit `preprint.tex` only.

## Build

```bash
brew install tectonic          # one-time; self-contained, fetches packages on demand
cd paper && tectonic -X compile preprint.tex
```

For the arXiv upload, run `python make_arxiv.py` and follow `ARXIV.md`.

### Build notes

- `SciPost.cls` and `SciPost_bibstyle.bst` give the single-column layout and the
  reference style that prints DOIs. They are SciPost's 2021-08 class files, taken from a mirror
  at gitlab.kwant-project.org/qt/skeleton. The paper is **not** a SciPost submission:
  `preprint.tex` loads the class with no journal option, clears the header and leaves out the
  copyright box and line numbers.
- The class's template forbids line breaks inside inline math
  (`\binoppenalty=\relpenalty=10000`). That produced six overfull lines, so `preprint.tex`
  relaxes it to 700/500 and sets `\emergencystretch`.
- The template's `\boldmath` abstract had no enclosing group, which turned every equation in the
  paper bold. Here it is scoped with braces.
- **Micro sign.** Under XeTeX (tectonic) with mathdesign's T1 Charter, siunitx's default `\micro`
  lands on T1 slot 181, which is "ţ", so every "µs" printed as "ţs". The preamble redefines
  `\micro` to TS1 slot 181 (a real upright µ) via `\upmu`, which also works under pdfLaTeX. Use
  `\upmu`, not `\mu`, for unit prefixes such as `\upmu\Phi_0`.
- Clean build: 0 overfull boxes, 0 unresolved citations or cross-references. Two underfull boxes
  in the reference list come from long DOIs and are cosmetic.
- Figures are copied into `paper/figures/` so the directory is self-contained. After re-running
  the notebooks, regenerate them with `cp ../figures/*.pdf figures/`.

## Rules for editing this paper

**Every number here must exist in `data/headline_provenance_ledger.csv` first.**
`05_summary.ipynb` §2 recomputes each ledger entry from its named file and filter and
asserts it, and `tests/test_regression.py` fails if the ledger shrinks or any row stops
verifying. The workflow is:

1. add the number to the ledger in `build_nb05.py`
2. re-run `python build_nb05.py && jupyter nbconvert --execute --inplace notebooks/05_summary.ipynb`
3. confirm it verifies
4. only then quote it here

**American spelling throughout** (the author's choice, 2026-10-01): modeling, optimization,
behavior, analog, Acknowledgments. Reference titles keep their published spelling (e.g. Papič
*et al.*, "…optimisation of transmon-qubit design parameters").

**The arXiv abstract is capped at 1920 characters** and is currently at 1919. Any word added to
the abstract has to be paid for with a cut. `ARXIV.md` holds the form version.

Do not copy numbers out of notebook prose — prose drifts, the ledger does not. The
superseded-number audit in `05_summary.ipynb` §4 scans `README.md` and `docs/` but **not**
`paper/`; if a headline number is ever retired, grep this directory by hand.

## Status

The text is the IEEE draft's body text, re-laid out for a single column, with a new back matter: a
*Use of generative AI* paragraph, a *Funding* paragraph (none) and a code and data availability
section. Every reference except the arXiv-only Ram 2025 has a DOI. All the headline numbers come
from the ledger.

### Before posting — the author has to do these

1. **Voice pass.** The prose was drafted from the committed artifacts, not written by the
   author. The physics and every number are checked, but the sentences are not yours yet. An arXiv
   version is permanent, so read it line by line first.
2. **Make the AI disclosure match what actually happened.** If AI help went beyond the draft text
   and the notebook interpretations, say so.
3. **Check every reference yourself.** Confirm that each DOI opens the right paper and that the
   paper supports the sentence citing it. Start with the 2025 ones: Lange, Bénâtre and Ram. *Yan 2016 was checked on 2026-10-01 against the Nat. Commun. PDF (p. 5):
   A_Φ = 1.4 µΦ₀, γ = 0.9.* *Somoroff 2023 was checked on
   2026-10-01 against the published PRL abstract:* it reports $T_2^* = 1.48 \pm 0.13$ ms and
   gate fidelity 0.99991(1). The draft had said "$T_2$ = 1.43 ms", which was wrong; it is now
   fixed in `preprint.tex`. The frozen `main.tex` still contains the error.
4. **Affiliation:** "Independent researcher". The author confirmed this on 2026-10-01.
5. **Endorsement.** First-time `quant-ph` submitters usually need one. See `ARXIV.md`.
6. *Optional:* archive a tagged release on Zenodo and cite its DOI next to the GitHub link.

### Limits stated in the paper itself, not hidden here

- **Novelty is argued, not proven.** §2.1 locates the two closest works — McRae *et al.* on
  experimental reproducibility of coherence comparisons, and Lange *et al.* on controlled
  transmon/fluxonium comparison for cross-talk rather than coherence — and §2.2 claims only
  "we are not aware of", stating the survey's incompleteness in the body text.
- **The 48× gate-count result rests on a coarse calibration.** Three parameters over a fixed
  Gaussian; durations are upper bounds; the optimiser is not converged at the shortest
  durations. §8.4 and §12 name GRAPE/CRAB as the next step.
- **Frequency matching is necessary but not sufficient** as a control: geometry and surface
  participation also set dielectric loss, and a Hamiltonian-level model does not capture
  them. Stated in §2.3, "A limitation of the present comparison".
- **Model-based throughout.** No new measurement; the calibration set is nine devices from
  one 2019 paper, and §9 notes the field has since moved past them.
