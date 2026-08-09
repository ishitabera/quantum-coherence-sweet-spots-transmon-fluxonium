---
name: check-physics
description: Sanity-audit a coherence result, CSV, notebook output, or draft claim against this project's units and physics invariants before it is trusted or written up. Use when a number looks surprising, before quoting a result in docs/ or paper/, after regenerating data, or when the user asks "is this right".
---

# Auditing a coherence result

Work through these in order. Stop and report as soon as something fails — a units bug
invalidates everything downstream, so there is no point checking physics on top of it.

## 1. Units

Every coherence column must be microseconds and named `*_us`. scqubits returns nanoseconds;
the conversion is `* 1e-3`.

Order-of-magnitude sanity: transmon T2 in this study lands in the tens-to-low-hundreds of µs;
fluxonium at the flux sweet spot reaches single-digit milliseconds (thousands of µs). A number
1000× outside its expected band is almost always a missing or doubled `NS_TO_US`.

## 2. The T2 inequality

`T2 ≤ 2·T1` always. Since `1/T2 = 1/(2T1) + Σ 1/Tφ_k` with non-negative rates, T2 can only
approach `2T1` from below. Check it elementwise across the whole CSV:

```python
assert (df.T2_us <= 2 * df.T1_us + 1e-9).all()
```

A violation means rates were combined wrong (times averaged or min'd instead of rates summed).

## 3. Infinities from `safe()`

`safe()` returns `inf` for an unsupported or failing channel, contributing zero rate. That is
intentional, but it means a *broken* channel looks like a *good* result. Count them:

```python
for c in [col for col in df.columns if col.endswith('_us')]:
    n = np.isinf(df[c]).sum()
    if n: print(c, 'inf at', n, 'points')
```

Any `inf` in a channel that should be supported for that architecture is a bug. Confirm each
one against `supported_noise_effects()` for that qubit.

## 4. Channel set matches the architecture

Transmon: dielectric T1, charge Tφ, critical-current Tφ — and **no flux channel**.
Fluxonium: dielectric + quasiparticle T1, flux Tφ, critical-current Tφ — and **no charge channel**.
A transmon result carrying a flux column is wrong; so is a fluxonium result carrying `Tphi_ng_us`.

## 5. Expected qualitative behaviour

These are the trends the study established. A result contradicting one of them needs an
explanation, not a shrug:

- Fluxonium T2 **peaks sharply at flux = 0.5** and falls off as roughly `(0.5 − Φ)^−0.68`.
  A flat fluxonium flux sweep means the sweet spot was missed or `flux` never varied.
- Transmon T2 is **nearly flat** against ng at large EJ/EC (bias sensitivity ≈ 1.01×) —
  that flatness *is* charge insensitivity, not a bug.
- Charge dephasing dominates the transmon only at **low EJ/EC**; dielectric relaxation
  dominates the usable regime.
- Fluxonium beats transmon by roughly two orders of magnitude at the sweet spot.

## 6. Regime honesty

If the number is a transmon headline, confirm it came from the **EJ/EC ≥ 20** subset. The
global maximum over the swept plane sits near EJ/EC ≈ 0.5, which is a charge qubit. Quoting
that as a transmon result is the single most misleading error available in this project.

## Reporting

State which checks passed, which failed, and for failures give the specific rows or parameter
values involved. If everything passes, say so plainly and name the file and filter you
verified, so the claim is traceable when it reaches `docs/` or `paper/`.
