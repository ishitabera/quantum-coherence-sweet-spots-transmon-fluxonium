"""Robustness checks added to the paper on 2026-10-14 in response to review.

    python results/referee_checks.py      (from the repo root, conda env `qubits`)

(a) Transmon quasiparticle-tunneling T1. scqubits.Transmon has no such channel, so the
    scqubits Smith-2020 spectral density is re-implemented here, validated against
    Fluxonium.t1_quasiparticle_tunneling at flux=0.3, and applied with the transmon
    <0|sin(phi/2)|1> from a phase-grid diagonalisation.
(b) Second-order 1/f charge dephasing at the transmon charge sweet spot ng=0, in the
    same Nguyen Eq. (A3) form the paper uses for flux: Gamma = |d2 omega/dng2| A_ng^2.
(c) The frequency confound under the scqubits frequency-dependent Q_cap default
    instead of the flat Q_cap=1e6.
(d) The thermal factor coth(h f01 / 2 kB T) in the total dielectric rate.

Numbers quoted in paper/preprint.tex (Sections 2.3, 3.2, 3.3, 7.3) come from this output.
"""
import numpy as np
import pandas as pd
import scipy.constants as C
import scipy.sparse
import scipy.sparse.linalg
import scipy.special
import scqubits as scq
from scqubits.core.noise import NOISE_PARAMS

scq.settings.T1_DEFAULT_WARNING = False
NS_TO_US = 1e-3
OUT = {}   # written to data/robustness_checks.json; the 05_summary ledger reads it
T, DELTA_EV, RK = NOISE_PARAMS["T"], NOISE_PARAMS["Delta"], C.h / C.e**2


# ---------------------------------------------------------------- (a) transmon QP
def s_qp(omega, EJ_GHz, x_qp):
    """scqubits' quasiparticle spectral density (Smith et al. 2020, Eqs. 19, S23)."""
    x = C.hbar * omega / (C.k * T)
    delta_hz, f_hz, ej_hz = DELTA_EV * C.e / C.h, omega / (2 * np.pi), EJ_GHz * 1e9
    re_y = (np.sqrt(2 / np.pi) * (8 / RK) * (ej_hz / delta_hz) * (2 * delta_hz / f_hz) ** 1.5
            * x_qp * np.sqrt(x / 2) * scipy.special.kv(0, x / 2) * np.sinh(x / 2))
    return 2 * C.hbar * omega / C.e**2 * re_y / np.tanh(x / 2) / (1 + np.exp(-x))


def t1_qp_us(me2, f_GHz, EJ_GHz, x_qp):
    """Total (down + up) T1 in microseconds, as scqubits computes with total=True."""
    w = 2 * np.pi * f_GHz * 1e9
    s = s_qp(w, EJ_GHz, x_qp)
    return 1e6 / (me2 * s * (1 + np.exp(-C.hbar * w / (C.k * T))))


fx = scq.Fluxonium(EJ=8.9, EC=2.5, EL=0.5, flux=0.3, cutoff=110)
ev, es = fx.eigensys(2)
op = fx.sin_phi_operator(alpha=0.5, beta=0.5 * 2 * np.pi * fx.flux)
me2 = abs(es[:, 1].conj() @ op @ es[:, 0]) ** 2
ours, ref = t1_qp_us(me2, ev[1] - ev[0], 8.9, 3e-6), fx.t1_quasiparticle_tunneling() * NS_TO_US
assert abs(ours / ref - 1) < 1e-6, (ours, ref)
print(f"(a) validation, fluxonium flux=0.3: {ours:.2f} us vs scqubits {ref:.2f} us")


def transmon_f01_me2(EJ, EC, n=4001):
    """f01 and |<0|sin(phi/2)|1>|^2 from a finite-difference phase-basis Hamiltonian."""
    phi = np.linspace(-np.pi, np.pi, n)
    d = phi[1] - phi[0]
    lap = scipy.sparse.diags([np.ones(n - 1), -2 * np.ones(n), np.ones(n - 1)], [-1, 0, 1]) / d**2
    H = -4 * EC * lap + scipy.sparse.diags(-EJ * np.cos(phi))
    w, v = scipy.sparse.linalg.eigsh(H.tocsc(), k=2, sigma=-EJ, which="LM")
    o = np.argsort(w)
    w, v = w[o], v[:, o]
    return w[1] - w[0], abs(v[:, 1] @ (np.sin(phi / 2) * v[:, 0])) ** 2


tm = scq.Transmon(EJ=13.0, EC=0.26, ng=0.0, ncut=30)
print(f"    transmon operating point, dielectric T1 = "
      f"{tm.t1_capacitive(i=1, j=0, Q_cap=1e6) * NS_TO_US:.1f} us")
for EJ, EC in [(13.0, 0.26), (6.0, 0.3)]:
    f, m = transmon_f01_me2(EJ, EC)
    row = ", ".join(f"x_qp={x:.0e}: {t1_qp_us(m, f, EJ, x):.3g} us" for x in (3e-6, 1e-7, 1e-8))
    print(f"    transmon EJ/EC={EJ / EC:.0f} (f01={f:.3f} GHz)  T1_qp  {row}")
    OUT.setdefault("transmon_T1_qp_us", {})[f"EJ_EC_{EJ / EC:.0f}"] = {
        f"x_qp_{x:.0e}": float(t1_qp_us(m, f, EJ, x)) for x in (3e-6, 1e-7, 1e-8)}

# ---------------------------------------------------------------- (b) 2nd-order charge
def f01(EJ, EC, ng):
    e = scq.Transmon(EJ=EJ, EC=EC, ng=ng, ncut=30).eigenvals(2)
    return e[1] - e[0]


A_NG, h = 1e-4, 1e-3
for EJ, EC in [(13.0, 0.26), (2.0, 0.1), (4.0, 0.2), (6.0, 0.3)]:
    d2 = (f01(EJ, EC, h) - 2 * f01(EJ, EC, 0.0) + f01(EJ, EC, -h)) / h**2   # GHz per ng^2
    gamma = abs(2 * np.pi * d2 * 1e9) * A_NG**2
    print(f"(b) EJ/EC={EJ / EC:4.0f}: second-order charge Tphi = {1 / gamma:.3g} s")
    OUT.setdefault("second_order_charge_Tphi_s", []).append(
        dict(EJ=EJ, EC=EC, EJ_EC=EJ / EC, Tphi_s=float(1 / gamma)))

# ---------------------------------------------------------------- (c) Q_cap(omega)
t = pd.read_csv("data/transmon_full_sweep.csv")
t = t[t.EJ_EC >= 20].copy()
t["T1_def"] = [scq.Transmon(EJ=r.EJ, EC=r.EC, ng=r.ng, ncut=30).t1_capacitive(i=1, j=0)
               * NS_TO_US for r in t.itertuples()]
t["T2_def"] = 1 / (1 / (2 * t.T1_def) + 1 / t.Tphi_ng_us + 1 / t.Tphi_cc_us)
f = pd.read_csv("data/fluxonium_design_sweep.csv")
t1d = np.array([scq.Fluxonium(EJ=r.EJ, EC=r.EC, EL=r.EL, flux=0.5, cutoff=110)
                .t1_capacitive(i=1, j=0) * NS_TO_US for r in f.itertuples()])
f["T1_def"] = 1 / (1 / t1d + 1 / f.T1_qp_us)
f["T2_def"] = 1 / (1 / (2 * f.T1_def) + 1 / f.Tphi_flux_us + 1 / f.Tphi_flux2_us + 1 / f.Tphi_cc_us)
for name, d in [("transmon, EJ/EC >= 20", t), ("fluxonium (EJ, EL)", f)]:
    for col, label in [("T2_us", "flat Q_cap=1e6"), ("T2_def", "scqubits Q_cap(omega)")]:
        r = np.corrcoef(np.log(d[col]), np.log(d.freq_GHz))[0, 1]
        at_fmin = d[col].idxmax() == d.freq_GHz.idxmin()
        print(f"(c) {name:22s} {label:22s} corr = {r:.4f}   max T2 at lowest f01: {at_fmin}")
        OUT.setdefault("confound_under_Qcap_model", []).append(
            dict(plane=name, Q_cap=label, corr=float(r), max_at_lowest_f01=bool(at_fmin)))

# ---------------------------------------------------------------- (d) thermal factor
for fr in (0.075, 0.363, 4.92):
    x = C.h * fr * 1e9 / (C.k * T)
    print(f"(d) coth(h f / 2 kB T) at {fr} GHz, T = {T * 1e3:.0f} mK: {1 / np.tanh(x / 2):.2f}")

# ======================================================================================
# Second review round (2026-10-14)
# ======================================================================================
import sys
sys.path.insert(0, "src")
from scipy.optimize import minimize_scalar
from coherence import transmon_coherence, fluxonium_coherence
from qubits import make_transmon, make_fluxonium, qubit_frequency

# ---------------------------------------------------------------- (e) analytic QP cross-check
# Catelani et al., PRB 84, 064517 (2011), leading order for a transmon at hf >> kT:
#   Gamma = x_qp * (omega_p / pi) * sqrt(2 Delta / (hbar omega_10)),  omega_p = sqrt(8 EJ EC)/hbar.
# Parity: QP tunneling connects the even and odd charge-parity sectors; at EJ/EC >= 20 the two
# sectors' wavefunctions coincide to the charge dispersion (< 1e-3 relative), which is why the
# extended-phase matrix element above and the analytic result agree.
for EJ, EC in [(13.0, 0.26), (6.0, 0.3)]:
    f10, _ = transmon_f01_me2(EJ, EC)
    wp = 2 * np.pi * np.sqrt(8 * EJ * EC) * 1e9
    w10 = 2 * np.pi * f10 * 1e9
    d_over_hbar = DELTA_EV * C.e / C.hbar
    g = 3e-6 * wp / np.pi * np.sqrt(2 * d_over_hbar / w10)
    print(f"(e) Catelani 2011 analytic, EJ/EC={EJ / EC:.0f}, x_qp=3e-6: T1_qp = {1e6 / g:.2f} us")
    OUT.setdefault("transmon_T1_qp_catelani_us", {})[f"EJ_EC_{EJ / EC:.0f}"] = float(1e6 / g)

# ---------------------------------------------------------------- (f) residuals at fixed f01
for name, d in [("transmon, EJ/EC >= 20", t), ("fluxonium (EJ, EL)", f)]:
    x, y = np.log(d.freq_GHz.values), np.log(d.T2_us.values)
    k, b = np.polyfit(x, y, 1)
    res = y - (k * x + b)
    print(f"(f) {name:22s} log-log slope {k:.3f}; residual T2 scatter about the fit: "
          f"rms {100 * (np.exp(res.std()) - 1):.1f}%, worst {100 * (np.exp(abs(res).max()) - 1):.0f}%")
    OUT.setdefault("residual_T2_scatter_at_fixed_f01", []).append(
        dict(plane=name, slope=float(k), rms_pct=float(100 * (np.exp(res.std()) - 1)),
             worst_pct=float(100 * (np.exp(abs(res).max()) - 1))))
fl = pd.read_csv("data/fluxonium_design_sweep.csv")
for lo, hi in [(0.17, 0.19), (0.30, 0.33), (0.50, 0.55)]:
    s = fl[(fl.freq_GHz >= lo) & (fl.freq_GHz <= hi)]
    if len(s) > 1:
        print(f"(f) fluxonium designs with f01 in [{lo}, {hi}] GHz: n={len(s)}, "
              f"T2 max/min = {s.T2_us.max() / s.T2_us.min():.3f}")

# ---------------------------------------------------------------- (g) comparison under Q_cap(omega)
def fluxonium_at_frequency(f_target, EJ=8.9, EC=2.5, EL_range=(0.2, 6.0), n=140):
    els = np.linspace(*EL_range, n)
    freqs = np.array([qubit_frequency(make_fluxonium(EJ=EJ, EC=EC, EL=float(e), flux=0.5))
                      for e in els])
    i = int(np.abs(freqs - f_target).argmin())
    return dict(EJ=EJ, EC=EC, EL=float(els[i]), flux=0.5, cutoff=110)


def transmon_at_frequency(f_target, Q_cap, min_EJ_EC=50, ncut=60):
    best = None
    for ec in np.linspace(0.01, 0.5, 160):
        ej = (f_target + ec) ** 2 / (8 * ec)
        if ej / ec < min_EJ_EC:
            continue
        r = transmon_coherence(EJ=float(ej), EC=float(ec), ng=0.0, ncut=ncut, Q_cap=Q_cap)
        if best is None or r["T2_us"] > best["T2_us"]:
            best = r
    return best


for label, Q in [("flat Q_cap=1e6", 1e6), ("scqubits Q_cap(omega)", None)]:
    fx0 = fluxonium_coherence(Q_cap=Q)
    tm0 = transmon_coherence(EJ=13.0, EC=0.26, ng=0.0, Q_cap=Q)
    out = [f"raw T2 ratio {fx0['T2_us'] / tm0['T2_us']:.1f}x"]
    rec = OUT.setdefault("comparison_under_Qcap_model", {}).setdefault(label, {})
    rec["raw_T2_ratio"] = float(fx0["T2_us"] / tm0["T2_us"])
    for ft in (0.5, 1.0, 2.0):
        fx = fluxonium_coherence(**fluxonium_at_frequency(ft), Q_cap=Q)
        tm = transmon_at_frequency(ft, Q)
        out.append(f"matched {ft} GHz {fx['T2_us'] / tm['T2_us']:.2f}x")
        rec[f"matched_T2_ratio_{ft}GHz"] = float(fx["T2_us"] / tm["T2_us"])
    print(f"(g) {label:22s} " + ", ".join(out))

# ---------------------------------------------------------------- (h) refined interior optimum
rows = []
for EC in (1.5, 2.5, 3.5):
    for EL in (0.2, 0.5, 1.0):
        g = pd.read_csv("data/fluxonium_EC_slice_scan.csv")
        g = g[(g.EC == EC) & (g.EL == EL)].sort_values("EJ").reset_index(drop=True)
        i = int(g.T2_us.idxmax())
        lo, hi = g.EJ[max(i - 1, 0)], g.EJ[min(i + 1, len(g) - 1)]
        r = minimize_scalar(lambda ej: -fluxonium_coherence(EJ=ej, EC=EC, EL=EL)["T2_us"],
                            bounds=(lo, hi), method="bounded", options={"xatol": 0.02})
        c = fluxonium_coherence(EJ=r.x, EC=EC, EL=EL)
        share = (1 / c["Tphi_flux2_us"]) / (1 / c["T2_us"])
        # A maximum on the last (or first) grid point is the edge of the scan, not an optimum:
        # T2 may still be rising. Such slices are kept in the CSV but left out of the summary.
        interior = 0 < i < len(g) - 1
        rows.append(dict(EC=EC, EL=EL, grid_EJ=g.EJ[i], grid_T2=g.T2_us[i], EJ=r.x,
                         T2=c["T2_us"], f01_MHz=1e3 * c["freq_GHz"], share=share,
                         interior=interior))
        print(f"(h) EC={EC} EL={EL}: grid EJ={g.EJ[i]:5.1f} T2={g.T2_us[i]:7.0f} -> refined "
              f"EJ={r.x:6.2f} T2={c['T2_us']:7.0f} us  f01={1e3 * c['freq_GHz']:5.1f} MHz  "
              f"second-order share of 1/T2 = {100 * share:.1f}%")
h_all = pd.DataFrame(rows)
h = h_all[h_all.interior]
print(f"(h) {len(h)} of {len(h_all)} slices have an interior optimum; "
      f"edge-limited: {h_all[~h_all.interior][['EC', 'EL']].to_dict('records')}")
print(f"(h) refined optima: f01 {h.f01_MHz.min():.1f}-{h.f01_MHz.max():.1f} MHz, "
      f"share {100 * h.share.min():.0f}-{100 * h.share.max():.0f}%, "
      f"T2 gain over grid {100 * (h.T2 / h.grid_T2 - 1).max():.1f}% at most")

h_all.to_csv("data/fluxonium_interior_refined.csv", index=False)
OUT["interior_optimum_refined"] = dict(
    n_interior=int(len(h)), n_slices=int(len(h_all)),
    f01_MHz_min=float(h.f01_MHz.min()), f01_MHz_max=float(h.f01_MHz.max()),
    share_of_total_rate_min=float(h.share.min()), share_of_total_rate_max=float(h.share.max()),
    csv="data/fluxonium_interior_refined.csv")

import json, subprocess, datetime
OUT["_metadata"] = dict(
    generated=datetime.date.today().isoformat(), script="results/referee_checks.py",
    git_revision=subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True,
                                text=True).stdout.strip(),
    scqubits_version=scq.__version__, noise_params={k: NOISE_PARAMS[k] for k in ("T", "x_qp", "Delta")},
    units="times in us unless the key says _s; frequencies as named")
json.dump(OUT, open("data/robustness_checks.json", "w"), indent=1)
print("wrote data/robustness_checks.json and data/fluxonium_interior_refined.csv")
