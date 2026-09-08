"""Pin the numbers the committed CSVs, figures and paper claims are built on.

If any of these move, every artifact in data/ and figures/ is invalidated and the
headline numbers in README.md, docs/ and paper/ must be regenerated. That is
occasionally the right thing to do -- but it must be a deliberate act, not a
side effect. These tests make it deliberate.

Anchors are sourced from 01_validation.ipynb and the metadata JSONs written by the
notebooks; the tolerances are numerical, not physical.
"""
import json
from pathlib import Path

import numpy as np
import pytest

import coherence as co
from qubits import FLUXONIUM_DEFAULTS, TRANSMON_DEFAULTS

ROOT = Path(__file__).resolve().parent.parent
REL = 1e-9


# ---------------------------------------------------------------- validated anchors
def test_transmon_nominal_matches_the_validation_record():
    """01_validation.ipynb section 6. Notebooks 02 and 04 assert against this too."""
    assert co.transmon_coherence(**TRANSMON_DEFAULTS)["T2_us"] == pytest.approx(
        54.35208346332045, rel=REL)


def test_fluxonium_nominal_first_order_only_matches_the_validation_record():
    """The chain of custody to 01_validation.ipynb runs through the first-order-only
    path, because notebook 01 predates the second-order flux channel."""
    r = co.fluxonium_coherence(**FLUXONIUM_DEFAULTS, second_order_flux=False)
    assert r["T2_us"] == pytest.approx(1136.8688060132743, rel=REL)


def test_fluxonium_nominal_with_second_order_matches_the_committed_pipeline():
    """03_fluxonium_sweep.ipynb section 2: adding the second-order channel moves T2
    from 1136.87 to 1129.71 us, i.e. it carries 0.63% of the rate budget."""
    r = co.fluxonium_coherence(**FLUXONIUM_DEFAULTS)
    assert r["T2_us"] == pytest.approx(1129.7123, rel=1e-6)
    share = (1 / r["Tphi_flux2_us"]) / (1 / r["T2_us"])
    assert share == pytest.approx(0.0063, abs=5e-4)


def test_operating_points_reproduce_notebook_02_and_04():
    tm = co.transmon_coherence(EJ=13.0, EC=0.26, ng=0.0, ncut=30)
    assert tm["freq_GHz"] == pytest.approx(4.925, abs=5e-3)
    assert tm["T2_us"] == pytest.approx(62.4244, rel=1e-5)
    fx = co.fluxonium_coherence(**FLUXONIUM_DEFAULTS)
    assert fx["freq_GHz"] == pytest.approx(0.3634, abs=5e-4)


def test_sweet_spot_flux_curvature():
    """Sets the second-order dephasing rate, so it sets the interior design optimum."""
    geom = {k: FLUXONIUM_DEFAULTS[k] for k in ("EJ", "EC", "EL", "cutoff")}
    assert float(co.d2f01_dflux2(**geom, flux=0.5)) == pytest.approx(886.8, abs=0.5)


# ---------------------------------------------------------------- committed artifacts
def test_committed_csvs_still_reproduce_from_the_pipeline():
    """Spot-check that data/ was produced by this code and has not drifted."""
    pd = pytest.importorskip("pandas")
    df = pd.read_csv(ROOT / "data" / "fluxonium_design_sweep.csv")
    row = df[np.isclose(df.EJ, 8.9, atol=0.4) & np.isclose(df.EL, 0.5, atol=0.05)]
    if not len(row):
        pytest.skip("nominal design not on the committed grid")
    r = co.fluxonium_coherence(EJ=float(row.EJ.iloc[0]), EC=2.5,
                               EL=float(row.EL.iloc[0]), flux=0.5)
    assert r["T2_us"] == pytest.approx(float(row.T2_us.iloc[0]), rel=1e-6)


def test_headline_ledger_is_fully_verified():
    """05_summary.ipynb section 2 asserts this too; here it is a standalone gate."""
    pd = pytest.importorskip("pandas")
    p = ROOT / "data" / "headline_provenance_ledger.csv"
    if not p.exists():
        pytest.skip("ledger not generated yet - run 05_summary.ipynb")
    led = pd.read_csv(p)
    # 28 rows as of the 2026-09-07 audit pass; ten were added for numbers that were quoted in
    # prose but had no ledger row. CLAUDE.md: every headline number goes in the ledger.
    assert len(led) >= 28, f"ledger shrank to {len(led)} rows - a headline number lost its source"
    assert bool(led.verified.all()), led[~led.verified].to_string()


def test_no_superseded_number_is_quoted_as_a_live_claim():
    pd = pytest.importorskip("pandas")
    p = ROOT / "data" / "superseded_number_audit.csv"
    if not p.exists():
        pytest.skip("audit not generated yet - run 05_summary.ipynb")
    aud = pd.read_csv(p)
    live = aud[aud.kind == "LIVE CLAIM"]
    assert len(live) == 0, live.to_string()


# ---------------------------------------------------------------- noise-model provenance
def test_scqubits_noise_defaults_have_not_changed_under_us():
    """Every committed number depends on these. scqubits has changed them between
    releases, which is why requirements.txt pins the version."""
    snap = co.noise_param_snapshot()
    for key, want in [("T", 0.015), ("x_qp", 3e-06), ("Delta", 0.00034),
                      ("t_exp", 10000.0), ("omega_low", 2 * np.pi * 1e-9)]:
        assert snap[key] == pytest.approx(want, rel=1e-6), f"{key}: {snap[key]} != {want}"


def test_project_noise_amplitudes_are_the_documented_ones():
    assert co.DEFAULTS == {"Q_cap": 1e6, "A_ng": 1e-4, "A_flux": 1e-6, "A_cc": 1e-7}


def test_flat_Q_cap_is_a_deliberate_override_of_the_frequency_dependent_default():
    """Passing Q_cap explicitly installs a FLAT quality factor, overriding scqubits'
    Q_cap(w) = 1e6*(2pi*6GHz/w)^0.7. At the fluxonium sweet spot that shortens T1 by
    ~7x, and it is the largest single modelling choice in the project (CLAUDE.md)."""
    from qubits import make_fluxonium
    fx = make_fluxonium(**FLUXONIUM_DEFAULTS)
    flat = fx.t1_capacitive(i=1, j=0, Q_cap=1e6)
    freq_dep = fx.t1_capacitive(i=1, j=0)
    assert freq_dep / flat == pytest.approx(7.1, abs=0.3)


# ---------------------------------------------------------------- published-number anchors
NGUYEN_DESIGNS = [("A", 3.00, 0.84, 1.00), ("B", 4.86, 0.84, 1.14), ("C", 2.20, 0.55, 0.72),
                  ("D", 2.20, 0.83, 0.52), ("E", 1.60, 0.86, 0.50), ("F", 3.40, 0.80, 0.41),
                  ("G", 1.65, 1.14, 0.19), ("H", 4.43, 1.00, 0.79), ("I", 3.43, 1.00, 0.58)]


def test_second_order_flux_clears_nguyens_published_bound():
    """The second-order flux channel is the only physics this project implements itself, so it
    is pinned against a published number rather than against its own formula.

    Nguyen et al., PRX 9, 041041 (2019), Appendix A 2: having measured A_flux = 1.8 uPhi_0 from
    off-sweet-spot echo, they state Eq. (A3) gives "an upper Tphi limit > 10 ms for all qubits",
    and Sec. III A gives "10-100 ms for the typical devices". Run our channel on their nine
    devices at their amplitude. Mirrors 01_validation.ipynb section 6.1.
    """
    A_NGUYEN = 1.8e-6
    times_ms = {}
    for dev, EJ, EC, EL in NGUYEN_DESIGNS:
        d2 = float(co.d2f01_dflux2(EJ=EJ, EC=EC, EL=EL, flux=0.5, cutoff=110))
        times_ms[dev] = float(co.tphi_flux_second_order(d2, A_flux=A_NGUYEN)) * co.NS_TO_US / 1e3
    worst = min(times_ms, key=times_ms.get)
    assert times_ms[worst] > 10.0, f"device {worst}: {times_ms[worst]:.1f} ms < Nguyen's 10 ms bound"
    assert times_ms[worst] == pytest.approx(15.4, abs=1.0)      # device B, the tightest


def test_manucharyan_published_parameter_sets_span_the_recorded_sweet_spot_range():
    """Manucharyan 2009 reports THREE parameter sets and they are not interchangeable: the
    rounded main-text values, and Table 1's junction- and extended-model least-squares fits.
    They span ~9.5% in f01 at the flux sweet spot, so "validated to 0.012%" is agreement with a
    rounding, not with the paper's fits. 01_validation.ipynb section 5 records this.
    """
    from qubits import make_fluxonium, qubit_frequency
    sets = {"main_text": dict(EJ=9.00, EC=2.50, EL=0.52),
            "table1_junction": dict(EJ=8.93, EC=2.39, EL=0.52),
            "table1_extended": dict(EJ=8.97, EC=2.47, EL=0.52)}
    f = {k: qubit_frequency(make_fluxonium(**v, flux=0.5, cutoff=110)) for k, v in sets.items()}
    assert f["main_text"] == pytest.approx(0.363378, abs=5e-5)
    assert f["table1_junction"] == pytest.approx(0.330056, abs=5e-5)
    assert f["table1_extended"] == pytest.approx(0.355097, abs=5e-5)
    spread = 100 * (max(f.values()) - min(f.values())) / np.mean(list(f.values()))
    assert spread == pytest.approx(9.53, abs=0.2)


def test_gate_count_pulse_bound_changes_the_q3_answer_by_two_orders_of_magnitude():
    """The gate-count metric N = T2 / T_gate needs a bound on T_gate, and the choice is worth
    ~390x for this fluxonium because |alpha| = 8.7 GHz exceeds its own f01 = 0.363 GHz by 24x.
    04_analysis.ipynb section 3 reports all three; this pins them so the headline cannot quietly
    revert to the |alpha| bound alone.
    """
    from qubits import anharmonicity, make_fluxonium, make_transmon
    tm = co.transmon_coherence(EJ=13.0, EC=0.26, ng=0.0, ncut=30)
    fx = co.fluxonium_coherence(**FLUXONIUM_DEFAULTS)
    a_tm = abs(anharmonicity(make_transmon(EJ=13.0, EC=0.26, ng=0.0, ncut=30)))
    a_fx = abs(anharmonicity(make_fluxonium(**FLUXONIUM_DEFAULTS)))

    def N(T2, alpha, f01, bound):
        bw = {"alpha": alpha, "both": min(alpha, f01), "f01": f01}[bound]
        return 1e3 * T2 * bw

    ratios = {b: (N(fx["T2_us"], a_fx, fx["freq_GHz"], b)
                  / N(tm["T2_us"], a_tm, tm["freq_GHz"], b))
              for b in ("alpha", "both", "f01")}
    assert ratios["alpha"] == pytest.approx(524.5, rel=1e-3)
    assert ratios["both"] == pytest.approx(22.01, rel=1e-3)
    assert ratios["f01"] == pytest.approx(1.335, rel=1e-3)
    # the bound is worth as much as the entire spread across metrics 1-3
    assert ratios["alpha"] / ratios["f01"] > 300
    # and for a transmon the bounds coincide, which is why this was easy to miss
    assert a_tm < tm["freq_GHz"] and a_fx > fx["freq_GHz"]
