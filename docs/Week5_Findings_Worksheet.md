# Lab Notebook

## Aug 4, 2026 — Finding 1 — Peak coherence
At their respective best case designs, the fluxonium reaches T2 ≈ 5578 µs (5.6 ms) versus the transmon's ≈ 187 µs (0.187 ms) which is about a 30× advantage. The 
origin is the fluxonium's flux sweet spot which is at Φ/Φ₀ = 0.5 the qubit frequency sits at a minimum, so dω/dΦ = 0 and the 1/f flux dephasing rate, which is 
proportional to that slope, vanishes to first order. The transmon has no comparable escape from its dominant limiter, dielectric loss which is a T1 energy-relaxation 
process with no bias point that nulls it so it remains pinned near the dielectric ceiling. This is consistent with the Week 1's observation that both qubits satisfy 
T2 = 2·T1 at their sweet spots, with pure dephasing designed away (1/T2 = 1/2T1 + 1/Tφ, and 1/Tφ → 0), each is relaxation-limited, and the fluxonium simply has the 
higher relaxation-limited ceiling. 

Implication: when the design goal is maximum coherence (ex. deep circuits requiring many gates before decoherence), the fluxonium is the better architecture.


## Aug 4, 2026 — Finding 2 — Robustness vs. fragility

As each qubit's bias knob is swept across its full range, transmon T2 changes by only 1.01× while fluxonium T2 changes by 528× so the two architectures sit at opposite 
extremes of bias tolerance. The reason is a difference in functional form. The transmon's charge sensitivity falls off exponentially with EJ/EC, so at EJ/EC = 50 it 
is negligibly small and offset charge noise has almost nothing to couple to. The fluxonium, by contrast, operates at a quadratic frequency minimum where it detunes 
from Φ = 0.5 the slope dω/dΦ grows linearly, so the flux-dephasing rate which scales as (dω/dΦ)² grows quadratically, and coherence collapses. In short, exponential 
suppression (transmon) versus polynomial blow-up (fluxonium). 

Implication: the transmon is far more tolerant of bias drift, offset fluctuations, and fabrication spread, which is why it is the more forgiving architecture and the 
one IBM and Google build their processors on, even though it gives up the fluxonium's peak coherence to get that robustness.


## Aug 4, 2026 — Finding 3 — Dominant-channel crossover (Q2)

At the charge sweet spot (ng = 0), the transmon is relaxation-limited across the entire EJ/EC plane, because the first-order charge sensitivity dω/dng vanishes there 
by symmetry. The charge dephasing is switched off, leaving dielectric loss as the sole limiter. Moving off the sweet spot to ng = 0.25 revives charge noise so it dominates 
for EJ/EC ≲ 25 and hands control back to relaxation above that. The boundary is a straight line of constant EJ/EC because the strength of charge dephasing is governed 
by the charge dispersion, which scales as exp(−√(8EJ/EC)), so equal-ratio designs share the same charge sensitivity. This is the same exponential-suppression mechanism 
responsible for the transmon's bias robustness in Finding 2. 

Implication: a designer can predict which noise source will limit a given transmon from a single number, the EJ/EC ratio, without running the full simulation.

## Aug 4, 2026 — Finding 4 — Sweet-spot scaling

Near Φ = 0.5, fluxonium T2 falls off as (flux detuning)^−0.68, dropping 10× within just 0.002 flux quanta of the sweet spot. Simple flux-noise theory predicts an exponent 
of −1, since the pure-dephasing time scales as Tφ ∝ 1/|dω/dΦ| ∝ 1/detuning at a quadratic minimum. The measured value is shallower because the fit spans two regimes, far 
from the sweet spot the qubit is flux-limited and approaches the −1 slope, but very close to it the flux-dephasing rate collapses and T2 is instead capped by relaxation 
and critical-current dephasing, which flattens the curve toward slope 0. The reported −0.68 is therefore an effective exponent averaged over the fit window [0.001, 0.05] 
so it characterizes the crossover between these regimes rather than a universal power law, and would shift with the chosen range. 

Implication: the usable coherence window around the sweet spot is extremely narrow which is around roughly ±0.002 flux quanta so the flux bias must be stabilized to 
about one part in 500, a concrete requirement for any fluxonium-based device.

## Aug 4, 2026 — Finding 5 — Robustness to model uncertainty 

Across 60 Monte Carlo draws sampling the noise amplitudes over their literature ranges, the fluxonium's best-case T2 exceeds the transmon's in 100% of samples, with a median 
ratio of 216× and a 90% confidence interval of [52×, 371×]. The absolute T2 values are uncertain because the true noise amplitudes (Q_cap, A_flux, A_cc, A_ng) are not known, 
they are literature estimates spanning orders of magnitude, which is precisely why the ratio's magnitude varies so widely across draws. Its direction, however, never flips. 
The entire distribution sits above parity (ratio = 1) because the fluxonium's advantage originates in the sweet-spot mechanism itself. The flux minimum that structurally nulls 
flux dephasing is a geometric feature of the design, and it is independent of how strong any individual noise source happens to be. 

Implication: the central conclusion "the fluxonium reaches higher peak coherence" is robust to the model's input uncertainty, even though no single absolute coherence 
time can be quoted with confidence. This is the honest way to report a model-based comparison: certain in direction, uncertain in magnitude.