# arXiv submission

The paper is going to arXiv only. There is no journal submission (2026-10-01).

Upload `paper/arxiv.tar.gz`, built by `cd paper && python make_arxiv.py`. Rebuild it after **any**
edit to `preprint.tex`, `refs.bib` or `figures/`. The tarball and `paper/arxiv/` are gitignored
build outputs.

**What the bundle contains:** `preprint.tex` with the comments stripped (arXiv publishes the
source), `preprint.bbl` (arXiv does not run BibTeX), the layout class `SciPost.cls`, its
`SciPost_bibstyle.bst` and the five figure PDFs. The paper carries no journal branding.

**What has been checked:** the bundle compiles using only its own files and the shipped `.bbl`, and
the result matches `preprint.pdf` text-for-text on all 21 pages, with 0 unresolved citations.
**Not checked:** arXiv compiles with pdfLaTeX, and this machine has only XeTeX (tectonic). The
source is pure ASCII and uses standard packages, so a problem is unlikely, but arXiv's compiled
preview is the real test. Read every page of it before you confirm.

## Form fields

| field | value |
|---|---|
| License | your choice. CC BY 4.0 keeps every journal option open later |
| Primary category | `quant-ph` |
| Cross-list | `cond-mat.supr-con` |
| Title | Coherence sweet spots are frequency minima in disguise: Metric sensitivity in transmon and fluxonium design spaces |
| Authors | Ishita Bera |
| Comments | 21 pages, 5 figures, 3 tables. Code and data: https://github.com/ishitabera/quantum-coherence-sweet-spots-transmon-fluxonium (archived: doi:10.5281/zenodo.23988521) |
| Report-no / Journal-ref / DOI | leave blank |

### Abstract (1914 / 1920 characters)

Keep the blank lines between paragraphs.

```
At fixed quality factor, dielectric loss makes a lower-frequency transmon-like qubit relax more slowly. This has a stronger consequence for design studies than is usually acknowledged. Under a common noise model, $\log T_2$ over the transmon $(E_J,E_C)$ plane correlates with $\log f_{01}$ at $r=-0.9999$ within the transmon regime, and over the fluxonium $(E_J,E_L)$ plane at $r=-0.9878$. Because modeled $T_2$ on each plane is largely a function of the transition frequency, a raw $T_2$ maximum primarily identifies the lowest-frequency corner of the swept domain rather than an optimal design. Both our initial headline optima lie at grid edges outside the usual bands.

Controlling for frequency changes these conclusions. For the transmon, we find no interior optimum in $E_J/E_C$ at fixed $f_{01}$: $T_2$ varies by only 3.2% over $E_J/E_C$ = 11-118, while the anharmonicity decreases by a factor of 7.3. Charge insensitivity is therefore bought with anharmonicity rather than coherence, and any reported optimum depends on the imposed regime boundary.

Within this model, the fluxonium has an interior optimum once second-order $1/f$ flux dephasing is included, but at only 5.4-30 MHz, a convention-dependent location, where $\bar n_{\rm th}$ = 10-57 requires active initialization. Across architectures there is no single answer: the fluxonium/transmon ratio spans $\sim$400$\times$ across the metrics we consider, and at matched frequency the transmon is ahead, with a median $T_2$ ratio of 0.37$\times$ at 1 GHz under our noise priors. A lab-frame pulse simulation puts the coherence-to-gate-time ratio at a model-dependent 17-27$\times$, between the two heuristic estimates, partly because the fluxonium drive couples 10.2$\times$ more strongly to its leakage transition than to its qubit transition. Critical-current noise, read as fractional rather than as scqubits applies it, dominates the fluxonium.
```

## Before you press submit

arXiv versions are permanent. You can replace or withdraw a version later, but the original stays
public. So do these first, in this order:

1. **Voice pass.** Read the paper as it will appear and make the sentences yours.
2. **AI disclosure and references.** Make the *Use of generative AI* paragraph match what
   actually happened, and check all 18 references yourself (see `README.md`).
3. **Affiliation:** "Independent researcher". The author confirmed this on 2026-10-01.
4. **Rebuild** with `python make_arxiv.py` if anything changed.
5. **Endorsement.** A first-time `quant-ph` submitter usually needs one. arXiv offers the
   endorsement request once you start the submission, and endorses some users automatically
   based on an institutional e-mail, so a psu.edu address on your arXiv account may help.
6. Then go to arxiv.org → *Start new submission*, upload the tarball, fill in the fields above,
   check the compiled preview, and submit.
