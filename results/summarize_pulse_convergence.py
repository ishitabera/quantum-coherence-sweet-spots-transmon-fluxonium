"""Build data/pulse_convergence.json from the raw convergence runs and the current T2 values.

    python results/summarize_pulse_convergence.py     (after results/pulse_convergence.py)

results/pulse_convergence.py writes one row per (architecture, case, duration) to
results/pulse_convergence.json. Until 2026-11-04 the summary that 05_summary.ipynb's ledger reads
was written by hand, so its gate-count ratio froze at the T2 values of the day and went stale when
critical-current noise was made fractional. This script derives every summary field from the rows
and from the operating-point T2 in data/q3_analysis_metadata.json, so it cannot drift.
"""
import json
from datetime import date

ROWS = json.load(open("results/pulse_convergence.json"))
Q3 = json.load(open("data/q3_analysis_metadata.json"))["operating_points"]
T2 = {a: Q3[a]["T2_us"] for a in ("transmon", "fluxonium")}


def shortest_pass(arch, case_pred):
    ts = [r["T_ns"] for r in ROWS if r["arch"] == arch and case_pred(r["case"]) and r["passes"]]
    return min(ts) if ts else None


def verdicts(arch, levels):
    """{duration: pass?} after recalibration on `levels` levels with the notebook's starts."""
    tag = f"{levels} levels"
    return {r["T_ns"]: r["passes"] for r in ROWS if r["arch"] == arch and r["case"] == tag
            or (r["arch"] == arch and levels == 6 and r["case"] == "6 levels, notebook starts")}


strong = {a: shortest_pass(a, lambda c: "7 starts" in c) for a in T2}
notebook = {a: shortest_pass(a, lambda c: c == "6 levels, notebook starts") for a in T2}
trunc_ok = all(verdicts(a, 6) == verdicts(a, n) for a in T2 for n in (8, 10)
               if verdicts(a, n))
tight = {(r["arch"], r["T_ns"]): r["passes"] for r in ROWS if "tight ODE" in r["case"]}
base = {(r["arch"], r["T_ns"]): r["passes"] for r in ROWS if r["case"] == "6 levels, notebook starts"}
ode_ok = all(base.get(k) == v for k, v in tight.items())


def ratio(tg):
    return (T2["fluxonium"] / tg["fluxonium"]) / (T2["transmon"] / tg["transmon"])


summary = {
    "truncation_converged": ("pass/fail verdict at every tested duration is unchanged at 6, 8 and 10 "
                             "levels after recalibration") if trunc_ok else "NOT converged - inspect rows",
    "ode_tolerance_converged": ("tighter atol/rtol reproduce every notebook pass/fail verdict"
                                if ode_ok else "NOT converged - inspect rows"),
    "notebook_Tgate_ns": notebook,
    "notebook_ratio": ratio(notebook),
    "stronger_optimizer_Tgate_ns": strong,
    "stronger_optimizer_ratio": ratio(strong),
    "T2_us_used": T2,
    "shorter_durations_not_retested": ("transmon < 6 ns and fluxonium < 2.5 ns were not re-run with "
                                       "the stronger optimizer"),
}
out = {"rows": ROWS, "summary": summary,
       "_metadata": {"script": "results/pulse_convergence.py (rows) + "
                               "results/summarize_pulse_convergence.py (summary)",
                     "generated": date.today().isoformat(),
                     "T2_source": "data/q3_analysis_metadata.json operating_points"}}
json.dump(out, open("data/pulse_convergence.json", "w"), indent=1)
print(json.dumps(summary, indent=1))
assert trunc_ok and ode_ok, "convergence checks failed; see rows"
