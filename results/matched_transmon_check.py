"""How much does the choice of matched-frequency transmon move the Q3 comparison?

    python results/matched_transmon_check.py      (from the project root; ~20 min)

04_analysis.ipynb matches each fluxonium to the transmon with the HIGHEST T2 at the same
f01, subject to EJ/EC >= 50, scanning EC over [0.01, 0.5] GHz. Modeled transmon T2 rises
slowly with EJ/EC at fixed f01, so that "best" transmon sits at the scan's EC floor
(EC = 10 MHz, EJ/EC = 325-5050, |alpha| ~ 10 MHz): a grid-edge choice, and the most
favorable transmon the scan allows. This script

  1. rebuilds the notebook's matched pairs and REPRODUCES its Monte Carlo exactly
     (same seed, same draw order, same bands) as a check on itself, then
  2. repeats both with the transmon held at the envelope floor, EJ/EC = 50,

and writes data/matched_transmon_check.json. The pairing logic is copied from
build_nb04.py on purpose, so a drift between the two is caught by step 1.
"""
import json
import sys
from datetime import datetime

import numpy as np
import pandas as pd

sys.path.insert(0, "src")
from qubits import make_fluxonium, qubit_frequency          # noqa: E402
from coherence import (transmon_coherence, fluxonium_coherence,  # noqa: E402
                       git_revision, noise_param_snapshot)

MATCH_FREQS = (0.5, 1.0, 2.0)
SEED, N_MC = 20260907, 400
BANDS = {"Q_cap": (5, 7), "A_flux": (-6, -5), "A_cc": (-7, -6), "A_ng": (-4, -3)}


def fluxonium_at_frequency(f_target, EJ=8.9, EC=2.5, EL_range=(0.2, 6.0), n=140):
    els = np.linspace(*EL_range, n)
    freqs = np.array([qubit_frequency(make_fluxonium(EJ=EJ, EC=EC, EL=float(e), flux=0.5))
                      for e in els])
    i = int(np.abs(freqs - f_target).argmin())
    return dict(EJ=EJ, EC=EC, EL=float(els[i]), flux=0.5, cutoff=110)


def transmon_best_T2(f_target, min_EJ_EC=50, ncut=60):
    """04_analysis.ipynb's transmon_at_frequency, verbatim in logic."""
    best = None
    for ec in np.linspace(0.01, 0.5, 160):
        ej = (f_target + ec) ** 2 / (8 * ec)
        if ej / ec < min_EJ_EC:
            continue
        r = transmon_coherence(EJ=float(ej), EC=float(ec), ng=0.0, ncut=ncut)
        if best is None or r["T2_us"] > best[0]["T2_us"]:
            best = (r, dict(EJ=float(ej), EC=float(ec), ncut=ncut))
    return best[1]


def transmon_at_ratio(f_target, ratio=50.0, ncut=60):
    """Same frequency relation as the notebook, EJ/EC fixed: (f + EC)^2/(8 EC) = r EC."""
    ec = f_target / (np.sqrt(8 * ratio) - 1)
    return dict(EJ=float(ratio * ec), EC=float(ec), ncut=ncut)


def draw(rng, mode):
    q_f = 10 ** rng.uniform(*BANDS["Q_cap"])
    c_f = 10 ** rng.uniform(*BANDS["A_cc"])
    q_t, c_t = (q_f, c_f) if mode == "shared" else (10 ** rng.uniform(*BANDS["Q_cap"]),
                                                    10 ** rng.uniform(*BANDS["A_cc"]))
    return dict(Q_cap_fx=q_f, A_cc_fx=c_f, Q_cap_tm=q_t, A_cc_tm=c_t,
                A_flux=10 ** rng.uniform(*BANDS["A_flux"]),
                A_ng=10 ** rng.uniform(*BANDS["A_ng"]))


def main():
    out = {"date": datetime.now().isoformat(timespec="seconds"), "git_revision": git_revision(),
           "noise_params": noise_param_snapshot(), "seed": SEED, "N_MC": N_MC, "bands": BANDS,
           "units": "T2 in microseconds; ratios fluxonium/transmon", "deterministic": []}

    pairs = {}
    for f in MATCH_FREQS:
        fx = fluxonium_at_frequency(f)
        tm_best, tm_50 = transmon_best_T2(f), transmon_at_ratio(f)
        fx_T2 = fluxonium_coherence(**fx)["T2_us"]
        rb = transmon_coherence(**tm_best, ng=0.0)
        r50 = transmon_coherence(**tm_50, ng=0.0)
        pairs[f] = (fx, tm_best, tm_50)
        row = dict(f01_target=f, fx_EL=fx["EL"], fx_T2=fx_T2,
                   tm_best_EJ_EC=tm_best["EJ"] / tm_best["EC"], tm_best_f01=rb["freq_GHz"],
                   tm_best_T2=rb["T2_us"], ratio_best=fx_T2 / rb["T2_us"],
                   tm50_EC=tm_50["EC"], tm50_f01=r50["freq_GHz"], tm50_T2=r50["T2_us"],
                   ratio_50=fx_T2 / r50["T2_us"])
        out["deterministic"].append(row)
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)

    fx, tm_best, tm_50 = pairs[1.0]
    rng = np.random.default_rng(SEED)
    rows = []
    for mode in ("shared", "independent"):
        for k in range(N_MC):
            s = draw(rng, mode)
            fx_m = fluxonium_coherence(**fx, Q_cap=s["Q_cap_fx"], A_flux=s["A_flux"],
                                       A_cc=s["A_cc_fx"])["T2_us"]
            kw = dict(ng=0.0, Q_cap=s["Q_cap_tm"], A_ng=s["A_ng"], A_cc=s["A_cc_tm"])
            tb = transmon_coherence(**tm_best, **kw)["T2_us"]
            t50 = transmon_coherence(**tm_50, **kw)["T2_us"]
            rows.append(dict(mode=mode, ratio_best=fx_m / tb, ratio_50=fx_m / t50))
            if k % 100 == 99:
                print(f"  {mode}: {k + 1}/{N_MC}", flush=True)
    mc = pd.DataFrame(rows)
    out["monte_carlo_1GHz"] = {
        mode: {"median_best": float(g.ratio_best.median()), "win_best": float((g.ratio_best > 1).mean()),
               "median_50": float(g.ratio_50.median()), "win_50": float((g.ratio_50 > 1).mean())}
        for mode, g in mc.groupby("mode")}

    # step 1: the copy must reproduce the committed Monte Carlo before step 2 means anything
    ref = pd.read_csv("data/q3_monte_carlo_summary.csv").set_index(["mode", "metric"])
    want = ref.loc[("shared", "T2 at matched f01 = 1 GHz")]
    got = out["monte_carlo_1GHz"]["shared"]
    out["reproduces_committed_mc"] = bool(abs(got["median_best"] - want["median"]) < 1e-9
                                          and abs(got["win_best"] - want["win_fraction"]) < 1e-12)
    print(json.dumps(out["monte_carlo_1GHz"], indent=2))
    print("reproduces committed MC:", out["reproduces_committed_mc"])
    with open("data/matched_transmon_check.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print("wrote data/matched_transmon_check.json")


if __name__ == "__main__":
    main()
