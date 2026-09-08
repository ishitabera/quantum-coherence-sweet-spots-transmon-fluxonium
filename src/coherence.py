"""Coherence estimation pipeline.

IMPORTANT UNITS: scqubits returns coherence times in NANOSECONDS.
Every function here converts to MICROSECONDS (multiply by 1e-3).

Rates add, not times:  1/T2 = 1/(2*T1) + sum_k 1/Tphi_k

SECOND-ORDER FLUX DEPHASING
---------------------------
scqubits implements 1/f dephasing to FIRST order only: its `tphi_1_over_f`
uses the diagonal matrix element of the noise operator, i.e. d(omega01)/d(lambda).
At a sweet spot that derivative vanishes by symmetry and scqubits returns a
divergent time. This is documented behaviour, not a bug -- see Groszkowski &
Koch, Quantum 5, 583 (2021), Sec. 5.1 footnote 4: "Currently scqubits returns
a value of numpy.inf at sweet spots. Higher order corrections will be added at
a later time."

The physical dephasing at the flux sweet spot is set by the SECOND-order term.
We implement the standard closed form used for fluxonium by

    Nguyen et al., Phys. Rev. X 9, 041041 (2019), Appendix, Eq. (A3):

        Gamma_Phi2 = (d^2 omega01 / d Phi^2) * A_flux^2

    which Nguyen cite to Ithier et al., Phys. Rev. B 72, 134519 (2005).

Both papers use the same 1/f normalisation as scqubits,
S(omega) = 2 pi A^2 / |omega|  (Nguyen Eq. A1; Groszkowski & Koch Eq. 2),
so A_flux carries the same meaning in all three and no rescaling is needed.

CONVENTION WARNING. Eq. (A3) carries NO logarithmic factor. Other treatments
of second-order 1/f dephasing attach a log-dependent prefactor of order
ln(omega_low * t_exp) or its square, which changes the answer by 1-2 orders of
magnitude. `tphi_flux_second_order` exposes this through `log_factor` so the
convention dependence can be reported rather than hidden. The default
(`log_factor="none"`) reproduces Nguyen Eq. (A3) exactly.
"""
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import scqubits as scq
from scqubits.core.noise import NOISE_PARAMS

from qubits import make_transmon, make_fluxonium

NS_TO_US = 1e-3

# Literature-range noise amplitudes used as defaults
DEFAULTS = dict(Q_cap=1e6, A_ng=1e-4, A_flux=1e-6, A_cc=1e-7)

# Physical constants for thermal occupation
H_PLANCK = 6.62607015e-34
K_B = 1.380649e-23


# --------------------------------------------------------------------------
# failure bookkeeping
# --------------------------------------------------------------------------
# `safe()` returning inf contributes ZERO rate to the T2 sum, so a silently
# broken channel looks like a *good* result. Every failure is therefore
# recorded here and notebooks are expected to assert this list is empty.
CHANNEL_FAILURES = []


def reset_channel_failures():
    """Clear the recorded channel-failure log. Call before a sweep."""
    CHANNEL_FAILURES.clear()


def safe(fn, _ctx=None, **kw):
    """Return inf if a noise channel is unsupported/fails -> contributes zero rate.

    The failure is appended to CHANNEL_FAILURES with the channel name, the
    exception, and the caller-supplied context, so it can never pass unnoticed.
    """
    try:
        return fn(**kw)
    except Exception as exc:                      # noqa: BLE001 - deliberate
        CHANNEL_FAILURES.append({
            "channel": getattr(fn, "__name__", repr(fn)),
            "error": f"{type(exc).__name__}: {exc}",
            "context": _ctx or {},
        })
        return np.inf


# --------------------------------------------------------------------------
# provenance
# --------------------------------------------------------------------------
def git_revision():
    """Short SHA of the current commit, plus a dirty flag. 'unknown' if not a repo."""
    try:
        root = Path(__file__).resolve().parent.parent
        sha = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=root,
            stderr=subprocess.DEVNULL).decode().strip()
        dirty = subprocess.check_output(
            ["git", "status", "--porcelain"], cwd=root,
            stderr=subprocess.DEVNULL).decode().strip()
        return sha + ("-dirty" if dirty else "")
    except Exception:                             # noqa: BLE001
        return "unknown"


def noise_param_snapshot():
    """Every scqubits NOISE_PARAMS default that silently enters these results.

    Q_cap / A_flux / A_cc are passed explicitly by this module, but T, x_qp,
    Delta, t_exp, omega_low and omega_high are scqubits defaults that are never
    named anywhere else and yet set the quasiparticle and 1/f results outright.
    """
    return {k: float(v) for k, v in NOISE_PARAMS.items()}


def n_thermal(f_GHz, T=None):
    """Mean thermal photon number of a mode at f_GHz.

    Defaults to the same temperature scqubits uses for its T1 channels
    (NOISE_PARAMS['T'], 15 mK), so the thermal statement and the relaxation
    numbers cannot drift apart.
    """
    T = float(NOISE_PARAMS["T"]) if T is None else T
    x = H_PLANCK * np.asarray(f_GHz, dtype=float) * 1e9 / (K_B * T)
    return 1.0 / (np.exp(x) - 1.0)


# --------------------------------------------------------------------------
# second-order 1/f flux dephasing
# --------------------------------------------------------------------------
def d2f01_dflux2(EJ, EC, EL, flux=0.5, cutoff=110, h=1e-4, richardson=True):
    """d^2 f01 / d(Phi/Phi0)^2 in GHz per Phi0^2, by central difference.

    The plain central second difference has error O(h^2). That is ample near the
    nominal design, but NOT in the deep-well regime reached at large EJ and small
    EL, where f01 collapses toward zero and the curvature grows by orders of
    magnitude: at EJ=40, EL=0.2 a fixed h=1e-4 is still 4% from converged.
    Because the second-order dephasing rate is directly proportional to this
    quantity, and because that rate is what creates the interior design optimum
    in section 6, the derivative is Richardson-extrapolated by default:

        d2 ~= (4 * d2(h/2) - d2(h)) / 3        [error O(h^4)]

    Set richardson=False for the plain two-point-stencil value.
    """
    cache = {}

    def f01(phi):
        key = round(float(phi), 12)
        if key not in cache:
            ev = make_fluxonium(EJ=EJ, EC=EC, EL=EL, flux=float(phi),
                                cutoff=cutoff).eigenvals(evals_count=2)
            cache[key] = ev[1] - ev[0]
        return cache[key]

    def second_diff(step):
        return (f01(flux + step) - 2.0 * f01(flux) + f01(flux - step)) / step ** 2

    d2_h = second_diff(h)
    if not richardson:
        return d2_h
    return (4.0 * second_diff(h / 2.0) - d2_h) / 3.0


def tphi_flux_second_order(d2f_dphi2, A_flux=1e-6, log_factor="none"):
    """Second-order 1/f flux dephasing time in NANOSECONDS (scqubits convention).

    Parameters
    ----------
    d2f_dphi2 : d^2 f01 / d(Phi/Phi0)^2 in GHz / Phi0^2 (ordinary frequency).
    A_flux    : 1/f flux noise amplitude in units of Phi0, same normalisation
                as scqubits and Nguyen 2019 Eq. (A.1).
    log_factor: which convention to apply.
        "none"  -> Gamma = A^2 * |d^2 omega/dPhi^2|          (Nguyen Eq. A3)
        "ln"    -> multiplied by |ln(omega_low * t_exp)|
        "ln2"   -> multiplied by ln(omega_low * t_exp)^2
        The latter two exist only so the notebook can report how much the
        answer depends on a convention the literature does not fix.

    Returns inf when the curvature is zero (no second-order dephasing).
    """
    d2w = 2.0 * np.pi * np.abs(np.asarray(d2f_dphi2, dtype=float))   # rad*GHz/Phi0^2
    ln = abs(np.log(NOISE_PARAMS["omega_low"] * NOISE_PARAMS["t_exp"]))
    factor = {"none": 1.0, "ln": ln, "ln2": ln ** 2}[log_factor]
    # d2w is in rad*1e9/s per Phi0^2; A^2 is dimensionless -> Gamma in 1e9 rad/s
    gamma_per_ns = factor * (A_flux ** 2) * d2w                      # 1/ns
    with np.errstate(divide="ignore"):
        return np.where(gamma_per_ns > 0, 1.0 / gamma_per_ns, np.inf)


# --------------------------------------------------------------------------
# transmon
# --------------------------------------------------------------------------
def transmon_coherence(EJ, EC, ng=0.0, ncut=30,
                       Q_cap=1e6, A_ng=1e-4, A_cc=1e-7):
    """T1, per-channel Tphi and combined T2 for one transmon, in microseconds.

    Channels: t1_capacitive (dielectric), tphi_1_over_f_ng (charge),
    tphi_1_over_f_cc (critical current). A plain transmon has NO flux loop,
    so tphi_1_over_f_flux is not supported and is intentionally absent.
    """
    ctx = dict(kind="transmon", EJ=EJ, EC=EC, ng=ng)
    tm = make_transmon(EJ=EJ, EC=EC, ng=ng, ncut=ncut)
    T1      = safe(tm.t1_capacitive,    _ctx=ctx, i=1, j=0, Q_cap=Q_cap)
    Tphi_ng = safe(tm.tphi_1_over_f_ng, _ctx=ctx, A_noise=A_ng, i=0, j=1)
    Tphi_cc = safe(tm.tphi_1_over_f_cc, _ctx=ctx, A_noise=A_cc, i=0, j=1)
    T2 = 1.0 / (1.0 / (2 * T1) + 1.0 / Tphi_ng + 1.0 / Tphi_cc)
    ev = tm.eigenvals(evals_count=2)
    return {"EJ": EJ, "EC": EC, "EJ_EC": EJ / EC, "ng": ng,
            "freq_GHz": ev[1] - ev[0],
            "T1_us": T1 * NS_TO_US,
            "Tphi_ng_us": Tphi_ng * NS_TO_US,
            "Tphi_cc_us": Tphi_cc * NS_TO_US,
            "T2_us": T2 * NS_TO_US}


# --------------------------------------------------------------------------
# fluxonium
# --------------------------------------------------------------------------
def fluxonium_coherence(EJ=8.9, EC=2.5, EL=0.5, flux=0.5, cutoff=110,
                        Q_cap=1e6, A_flux=1e-6, A_cc=1e-7,
                        second_order_flux=True, log_factor="none",
                        deriv_step=1e-4):
    """T1, per-channel Tphi and combined T2 for one fluxonium, in microseconds.

    T1 channels: dielectric + quasiparticle tunneling.
    Tphi channels: 1/f flux (first order, scqubits) + 1/f flux (second order,
    Nguyen 2019 Eq. A3) + critical current.

    Set second_order_flux=False to reproduce the first-order-only pipeline that
    produced the pre-existing committed results (used by 01_validation.ipynb to
    keep the chain of custody intact).
    """
    ctx = dict(kind="fluxonium", EJ=EJ, EC=EC, EL=EL, flux=flux)
    fx = make_fluxonium(EJ=EJ, EC=EC, EL=EL, flux=flux, cutoff=cutoff)
    T1_diel = safe(fx.t1_capacitive,                _ctx=ctx, i=1, j=0, Q_cap=Q_cap)
    T1_qp   = safe(fx.t1_quasiparticle_tunneling,   _ctx=ctx, i=1, j=0)
    T1 = 1.0 / (1.0 / T1_diel + 1.0 / T1_qp)
    Tphi_flux = safe(fx.tphi_1_over_f_flux, _ctx=ctx, A_noise=A_flux, i=0, j=1)
    Tphi_cc   = safe(fx.tphi_1_over_f_cc,   _ctx=ctx, A_noise=A_cc,   i=0, j=1)

    if second_order_flux:
        d2 = float(d2f01_dflux2(EJ, EC, EL, flux=flux, cutoff=cutoff, h=deriv_step))
        Tphi_flux2 = float(tphi_flux_second_order(d2, A_flux=A_flux,
                                                  log_factor=log_factor))
    else:
        d2, Tphi_flux2 = np.nan, np.inf

    T2 = 1.0 / (1.0 / (2 * T1) + 1.0 / Tphi_flux
                + 1.0 / Tphi_flux2 + 1.0 / Tphi_cc)
    ev = fx.eigenvals(evals_count=2)
    return {"EJ": EJ, "EC": EC, "EL": EL, "flux": flux,
            "freq_GHz": ev[1] - ev[0],
            "T1_us": T1 * NS_TO_US,
            "T1_diel_us": T1_diel * NS_TO_US,
            "T1_qp_us": T1_qp * NS_TO_US,
            "Tphi_flux_us": Tphi_flux * NS_TO_US,
            "Tphi_flux2_us": Tphi_flux2 * NS_TO_US,
            "d2f01_dflux2_GHz": d2,
            "Tphi_cc_us": Tphi_cc * NS_TO_US,
            "T2_us": T2 * NS_TO_US}


def dominant_channel(row, kind="transmon"):
    """Which noise channel has the largest rate at this point."""
    if kind == "transmon":
        rates = {"relaxation(T1)": 1 / (2 * row.T1_us),
                 "charge":         1 / row.Tphi_ng_us,
                 "crit-current":   1 / row.Tphi_cc_us}
    else:
        rates = {"relaxation":     1 / (2 * row.T1_us),
                 "flux":           1 / row.Tphi_flux_us,
                 "flux-2nd":       1 / getattr(row, "Tphi_flux2_us", np.inf),
                 "crit-current":   1 / row.Tphi_cc_us}
    return max(rates, key=rates.get)


def transmon_sweep(EJ_grid, EC_grid, ng=0.0, **kw):
    """Full EJ x EC sweep -> DataFrame. Save with .to_csv() BEFORE plotting."""
    rows = []
    for EJ in EJ_grid:
        for EC in EC_grid:
            try:
                rows.append(transmon_coherence(EJ, EC, ng=ng, **kw))
            except Exception:                     # noqa: BLE001
                rows.append({"EJ": EJ, "EC": EC, "EJ_EC": EJ / EC, "ng": ng,
                             "freq_GHz": np.nan, "T1_us": np.nan,
                             "Tphi_ng_us": np.nan, "Tphi_cc_us": np.nan,
                             "T2_us": np.nan})
    return pd.DataFrame(rows)


def fluxonium_flux_sweep(flux_values, **kw):
    """Sweep flux bias through the sweet spot -> DataFrame."""
    df = pd.DataFrame([fluxonium_coherence(flux=f, **kw) for f in flux_values])
    df["detuning"] = 0.5 - df["flux"]
    return df
