"""Coherence estimation pipeline.

IMPORTANT UNITS: scqubits returns coherence times in NANOSECONDS.
Every function here converts to MICROSECONDS (multiply by 1e-3).

Rates add, not times:  1/T2 = 1/(2*T1) + sum_k 1/Tphi_k
"""
import numpy as np
import pandas as pd

from qubits import make_transmon, make_fluxonium

NS_TO_US = 1e-3

# Literature-range noise amplitudes used as defaults
DEFAULTS = dict(Q_cap=1e6, A_ng=1e-4, A_flux=1e-6, A_cc=1e-7)


def safe(fn, **kw):
    """Return inf if a noise channel is unsupported/fails -> contributes zero rate."""
    try:
        return fn(**kw)
    except Exception:
        return np.inf


def transmon_coherence(EJ, EC, ng=0.0, ncut=30,
                       Q_cap=1e6, A_ng=1e-4, A_cc=1e-7):
    """T1, per-channel Tphi and combined T2 for one transmon, in microseconds.

    Channels: t1_capacitive (dielectric), tphi_1_over_f_ng (charge),
    tphi_1_over_f_cc (critical current). A plain transmon has NO flux loop,
    so tphi_1_over_f_flux is not supported and is intentionally absent.
    """
    tm = make_transmon(EJ=EJ, EC=EC, ng=ng, ncut=ncut)
    T1      = safe(tm.t1_capacitive,     i=1, j=0, Q_cap=Q_cap)
    Tphi_ng = safe(tm.tphi_1_over_f_ng,  A_noise=A_ng, i=0, j=1)
    Tphi_cc = safe(tm.tphi_1_over_f_cc,  A_noise=A_cc, i=0, j=1)
    T2 = 1.0 / (1.0 / (2 * T1) + 1.0 / Tphi_ng + 1.0 / Tphi_cc)
    ev = tm.eigenvals(evals_count=2)
    return {"EJ": EJ, "EC": EC, "EJ_EC": EJ / EC, "ng": ng,
            "freq_GHz": ev[1] - ev[0],
            "T1_us": T1 * NS_TO_US,
            "Tphi_ng_us": Tphi_ng * NS_TO_US,
            "Tphi_cc_us": Tphi_cc * NS_TO_US,
            "T2_us": T2 * NS_TO_US}


def fluxonium_coherence(EJ=8.9, EC=2.5, EL=0.5, flux=0.5, cutoff=110,
                        Q_cap=1e6, A_flux=1e-6, A_cc=1e-7):
    """T1, per-channel Tphi and combined T2 for one fluxonium, in microseconds.

    T1 channels: dielectric + quasiparticle tunneling.
    Tphi channels: 1/f flux + critical current.
    """
    fx = make_fluxonium(EJ=EJ, EC=EC, EL=EL, flux=flux, cutoff=cutoff)
    T1_diel = safe(fx.t1_capacitive, i=1, j=0, Q_cap=Q_cap)
    T1_qp   = safe(fx.t1_quasiparticle_tunneling, i=1, j=0)
    T1 = 1.0 / (1.0 / T1_diel + 1.0 / T1_qp)
    Tphi_flux = safe(fx.tphi_1_over_f_flux, A_noise=A_flux, i=0, j=1)
    Tphi_cc   = safe(fx.tphi_1_over_f_cc,   A_noise=A_cc,   i=0, j=1)
    T2 = 1.0 / (1.0 / (2 * T1) + 1.0 / Tphi_flux + 1.0 / Tphi_cc)
    ev = fx.eigenvals(evals_count=2)
    return {"EJ": EJ, "EC": EC, "EL": EL, "flux": flux,
            "freq_GHz": ev[1] - ev[0],
            "T1_us": T1 * NS_TO_US,
            "T1_diel_us": T1_diel * NS_TO_US,
            "T1_qp_us": T1_qp * NS_TO_US,
            "Tphi_flux_us": Tphi_flux * NS_TO_US,
            "Tphi_cc_us": Tphi_cc * NS_TO_US,
            "T2_us": T2 * NS_TO_US}


def dominant_channel(row, kind="transmon"):
    """Which noise channel has the largest rate at this point."""
    if kind == "transmon":
        rates = {"relaxation(T1)": 1 / (2 * row.T1_us),
                 "charge":         1 / row.Tphi_ng_us,
                 "crit-current":   1 / row.Tphi_cc_us}
    else:
        rates = {"relaxation":   1 / (2 * row.T1_us),
                 "flux":         1 / row.Tphi_flux_us,
                 "crit-current": 1 / row.Tphi_cc_us}
    return max(rates, key=rates.get)


def transmon_sweep(EJ_grid, EC_grid, ng=0.0, **kw):
    """Full EJ x EC sweep -> DataFrame. Save with .to_csv() BEFORE plotting."""
    rows = []
    for EJ in EJ_grid:
        for EC in EC_grid:
            try:
                rows.append(transmon_coherence(EJ, EC, ng=ng, **kw))
            except Exception:
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
