# arXiv submission

The paper is going to arXiv only. There is no journal submission (2026-10-01).

Upload `paper/arxiv.tar.gz`, built by `cd paper && python make_arxiv.py`. Rebuild it after **any**
edit to `preprint.tex`, `refs.bib` or `figures/`. The tarball and `paper/arxiv/` are gitignored
build outputs.

**What the bundle contains:** `preprint.tex` with the comments stripped (arXiv publishes the
source), `preprint.bbl` (arXiv does not run BibTeX), the layout class `SciPost.cls`, its
`SciPost_bibstyle.bst` and the five figure PDFs. The paper carries no journal branding.

**What has been checked:** the bundle compiles using only its own files and the shipped `.bbl`, and
the result matches `preprint.pdf` text-for-text on all 16 pages, with 0 unresolved citations.
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
| Comments | 16 pages, 5 figures, 3 tables. Code and data: https://github.com/ishitabera/quantum-coherence-sweet-spots-transmon-fluxonium |
| Report-no / Journal-ref / DOI | leave blank |

### Abstract (1898 / 1920 characters)

Keep the blank lines between paragraphs.

```
Dielectric loss is well known to scale with transition frequency, implying that a lower-frequency qubit relaxes more slowly at fixed quality factor. We show that this has a stronger consequence for design studies than is usually acknowledged. Under a common noise model, $\log T_2$ over the transmon $(E_J,E_C)$ plane correlates with $\log f_{01}$ at $r=-0.9999$ within the transmon regime, while over the fluxonium $(E_J,E_L)$ plane the correlation is $r=-0.9973$. Because each parameter plane effectively reduces to a single degree of freedom, the transition frequency, a raw $T_2$ maximum primarily identifies the lowest-frequency corner of the swept domain rather than an intrinsic optimal design. Indeed, both of our initial headline optima occur at grid-edge points outside any usable frequency band.

Controlling for frequency changes these conclusions. For the transmon, we find no interior optimum in $E_J/E_C$ at fixed $f_{01}$. $T_2$ varies by only 4.8% across a tenfold range in $E_J/E_C$, while the anharmonicity decreases by a factor of 7.3. Charge insensitivity is therefore obtained primarily at the cost of gate speed rather than through improved coherence, and any reported optimum depends on the imposed regime boundary.

The fluxonium has a genuine interior optimum once second-order $1/f$ flux dephasing is included, but it occurs at frequencies of only 10-73 MHz, where $\bar n_{\rm th}$ = 4-32. For the cross-architecture comparison, there is no single answer: the fluxonium/transmon ratio spans $\sim$500$\times$ across defensible metrics, while at matched frequency the $T_2$ advantage reverses to a median 0.86$\times$ across the literature noise band. A lab-frame pulse simulation yields a 36-58$\times$ advantage in gate count, excluding both analytic bounds: the fluxonium drive couples 10.2$\times$ more strongly to its leakage transition than to its qubit transition.
```

## Before you press submit

arXiv versions are permanent. You can replace or withdraw a version later, but the original stays
public. So do these first, in this order:

1. **Voice pass.** Read the paper as it will appear and make the sentences yours.
2. **AI disclosure and references.** Make the *Use of generative AI* paragraph match what
   actually happened, and check all 16 references yourself (see `README.md`).
3. **Affiliation:** "Independent researcher". The author confirmed this on 2026-10-01.
4. **Rebuild** with `python make_arxiv.py` if anything changed.
5. **Endorsement.** A first-time `quant-ph` submitter usually needs one. arXiv offers the
   endorsement request once you start the submission, and endorses some users automatically
   based on an institutional e-mail, so a psu.edu address on your arXiv account may help.
6. Then go to arxiv.org → *Start new submission*, upload the tarball, fill in the fields above,
   check the compiled preview, and submit.
