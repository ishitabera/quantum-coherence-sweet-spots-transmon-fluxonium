---
name: add-noise-channel
description: Add or modify a noise channel in the coherence pipeline (e.g. inductive loss, Purcell, photon shot noise, second-order flux) so it enters T1/Tphi/T2 correctly for the right architecture. Use when the user wants a new decoherence mechanism modelled, wants to change a noise amplitude, or asks why a channel is missing.
---

# Adding a noise channel

## First check whether scqubits supports it for that qubit

Do not assume. List what the object actually offers:

```python
import scqubits as scq
q = scq.Fluxonium(EJ=8.9, EC=2.5, EL=0.5, flux=0.5, cutoff=110)
print(q.supported_noise_effects())
```

A channel that is unsupported for an architecture must be **omitted**, not stubbed. The
canonical case: `tphi_1_over_f_flux` does not apply to a plain transmon because it has no
flux loop. Omitting it is a physics statement, and CLAUDE.md records it as intentional.

## Wire it in at three places, consistently

1. **`src/coherence.py`** — call it through `safe()` so a failure yields `inf` (zero rate),
   add its amplitude as a keyword argument with a literature default, and add it to the
   returned dict as `<name>_us` after multiplying by `NS_TO_US`.
2. **The T1 or T2 combination** — this is where mistakes happen. Fold it in as a *rate*:
   - a new **T1** channel joins the T1 sum: `1/T1 = 1/T1_a + 1/T1_b + ...`
   - a new **Tφ** channel joins the T2 sum: `1/T2 = 1/(2*T1) + Σ 1/Tφ_k`
   Never combine by taking a min or an average.
3. **`dominant_channel()`** — add its rate to the dict for the right architecture, and update
   the label tuple passed to `plotting.dominant_map()`. `CHANNEL_COLORS` currently holds three
   colors; a fourth channel needs a fourth color or the categorical map will misrender.

If the notebooks define their own inline copy of the coherence function (02 does), update it
there too or note explicitly that the notebook now differs from `src/`.

## Amplitudes

New amplitudes go in `DEFAULTS` in `src/coherence.py` alongside `Q_cap`, `A_ng`, `A_flux`,
`A_cc`. Cite where the value came from in a comment — these are literature ranges and the
Monte Carlo robustness analysis in `04_analysis.ipynb` samples over them, so an
undocumented amplitude weakens the headline uncertainty claim.

## Verifying the change did what you think

A new channel can only ever *reduce* T2 (adding a positive rate). If any T2 goes up, the wiring
is wrong — most likely the channel was added to the wrong sum, or it is returning `inf`.

Check that explicitly: recompute a handful of points before and after and assert
`T2_after <= T2_before` elementwise. Then report how much the new channel actually moved the
result and where it dominates, so the user can judge whether it is worth keeping.

## Consequences to flag

Adding a channel invalidates the committed CSVs, figures, and every headline number in
`README.md`, `docs/`, and `paper/`. Say which artifacts need regenerating; do not quietly
regenerate committed results unless asked.
