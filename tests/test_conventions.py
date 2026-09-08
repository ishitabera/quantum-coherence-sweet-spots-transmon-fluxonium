"""Lock down the conventions in CLAUDE.md that have actually cost this project time.

Each test here corresponds to a rule stated in CLAUDE.md under "Non-negotiable
conventions". They are cheap, and they exist so that a future edit to src/ cannot
silently break an invariant that every committed CSV depends on.
"""
import numpy as np
import pytest
import scqubits as scq

import coherence as co
from qubits import (FLUXONIUM_DEFAULTS, TRANSMON_DEFAULTS, anharmonicity,
                    make_fluxonium, make_transmon, qubit_frequency)


# ---------------------------------------------------------------- units
def test_scqubits_returns_nanoseconds_and_we_convert_to_microseconds():
    """The single easiest mistake in this project: a 1000x units bug."""
    tm = make_transmon(**{k: v for k, v in TRANSMON_DEFAULTS.items()})
    raw_ns = tm.t1_capacitive(i=1, j=0, Q_cap=1e6)
    got_us = co.transmon_coherence(**TRANSMON_DEFAULTS)["T1_us"]
    assert co.NS_TO_US == 1e-3
    assert got_us == pytest.approx(raw_ns * 1e-3, rel=1e-12)
    # sanity: a transmon T1 in microseconds is tens, not tens of thousands
    assert 1.0 < got_us < 1e4


# ---------------------------------------------------------------- rates add
@pytest.mark.parametrize("kind", ["transmon", "fluxonium"])
def test_rates_add_never_times(kind):
    """1/T2 = 1/(2*T1) + sum_k 1/Tphi_k. Never average or min() coherence times."""
    if kind == "transmon":
        r = co.transmon_coherence(**TRANSMON_DEFAULTS)
        rate = 1 / (2 * r["T1_us"]) + 1 / r["Tphi_ng_us"] + 1 / r["Tphi_cc_us"]
    else:
        r = co.fluxonium_coherence(**FLUXONIUM_DEFAULTS)
        rate = (1 / (2 * r["T1_us"]) + 1 / r["Tphi_flux_us"]
                + 1 / r["Tphi_flux2_us"] + 1 / r["Tphi_cc_us"])
    assert 1 / r["T2_us"] == pytest.approx(rate, rel=1e-12)


@pytest.mark.parametrize("kind", ["transmon", "fluxonium"])
def test_T2_never_exceeds_2T1(kind):
    """T2 <= 2*T1 always; equality only if every pure-dephasing channel is infinite."""
    r = (co.transmon_coherence(**TRANSMON_DEFAULTS) if kind == "transmon"
         else co.fluxonium_coherence(**FLUXONIUM_DEFAULTS))
    assert r["T2_us"] <= 2 * r["T1_us"] * (1 + 1e-12)


def test_fluxonium_T1_combines_dielectric_and_quasiparticle_as_rates():
    r = co.fluxonium_coherence(**FLUXONIUM_DEFAULTS)
    assert 1 / r["T1_us"] == pytest.approx(1 / r["T1_diel_us"] + 1 / r["T1_qp_us"], rel=1e-12)


# ---------------------------------------------------------------- channels per architecture
def test_a_plain_transmon_has_no_flux_loop():
    """tphi_1_over_f_flux is not a valid channel for scq.Transmon. Its absence is
    intentional, not an oversight (CLAUDE.md)."""
    keys = co.transmon_coherence(**TRANSMON_DEFAULTS)
    assert "Tphi_flux_us" not in keys and "Tphi_flux2_us" not in keys
    assert "tphi_1_over_f_flux" not in make_transmon().supported_noise_channels()


def test_fluxonium_supports_flux_and_quasiparticle_channels():
    ch = make_fluxonium().supported_noise_channels()
    assert "tphi_1_over_f_flux" in ch and "t1_quasiparticle_tunneling" in ch
    r = co.fluxonium_coherence(**FLUXONIUM_DEFAULTS)
    assert {"Tphi_flux_us", "Tphi_flux2_us", "T1_qp_us"} <= set(r)


# ---------------------------------------------------------------- safe() bookkeeping
def test_safe_returns_inf_and_records_the_failure():
    """`safe()` returning inf contributes ZERO rate, so a broken channel looks like a
    GOOD result. Every failure must therefore be logged (CLAUDE.md)."""
    co.reset_channel_failures()
    assert co.CHANNEL_FAILURES == []

    def boom(**kw):
        raise RuntimeError("channel exploded")

    out = co.safe(boom, _ctx={"where": "unit test"})
    assert np.isinf(out)
    assert len(co.CHANNEL_FAILURES) == 1
    rec = co.CHANNEL_FAILURES[0]
    assert "RuntimeError" in rec["error"] and rec["context"] == {"where": "unit test"}
    co.reset_channel_failures()
    assert co.CHANNEL_FAILURES == []


def test_a_clean_sweep_logs_no_failures():
    co.reset_channel_failures()
    co.fluxonium_coherence(**FLUXONIUM_DEFAULTS)
    co.transmon_coherence(**TRANSMON_DEFAULTS)
    assert co.CHANNEL_FAILURES == [], co.CHANNEL_FAILURES


# ---------------------------------------------------------------- second-order flux dephasing
def test_second_order_flux_matches_nguyen_eq_A3_analytically():
    """Gamma = A^2 * d^2(omega01)/dPhi^2, with omega = 2*pi*f. No log factor by default."""
    d2f, A = 886.8, 1e-6                      # GHz/Phi0^2, Phi0
    expected_ns = 1.0 / (A ** 2 * 2 * np.pi * d2f)
    assert float(co.tphi_flux_second_order(d2f, A_flux=A)) == pytest.approx(expected_ns, rel=1e-12)


def test_second_order_log_conventions_are_ordered_and_documented():
    """The literature does not fix this convention; the code must expose it, not hide it."""
    d2f = 886.8
    none = float(co.tphi_flux_second_order(d2f, log_factor="none"))
    ln = float(co.tphi_flux_second_order(d2f, log_factor="ln"))
    ln2 = float(co.tphi_flux_second_order(d2f, log_factor="ln2"))
    assert none > ln > ln2                     # bigger prefactor -> shorter time
    with pytest.raises(KeyError):
        co.tphi_flux_second_order(d2f, log_factor="bogus")


def test_zero_curvature_gives_no_second_order_dephasing():
    assert np.isinf(float(co.tphi_flux_second_order(0.0)))


def test_richardson_second_derivative_is_step_independent():
    """The plain central difference is not converged in the deep-well regime, and the
    second-order rate is directly proportional to this derivative."""
    geom = {k: FLUXONIUM_DEFAULTS[k] for k in ("EJ", "EC", "EL", "cutoff")}
    vals = [float(co.d2f01_dflux2(**geom, flux=0.5, h=h)) for h in (1e-3, 5e-4, 2e-4, 1e-4)]
    assert (max(vals) - min(vals)) / abs(np.mean(vals)) < 1e-5


def test_second_order_can_be_switched_off_for_the_validated_code_path():
    """01_validation.ipynb's chain of custody runs through the first-order-only path."""
    r = co.fluxonium_coherence(**FLUXONIUM_DEFAULTS, second_order_flux=False)
    assert np.isinf(r["Tphi_flux2_us"]) and np.isnan(r["d2f01_dflux2_GHz"])


# ---------------------------------------------------------------- thermal occupation
def test_n_thermal_matches_bose_einstein_and_uses_the_scqubits_temperature():
    from scqubits.core.noise import NOISE_PARAMS
    h, kB = 6.62607015e-34, 1.380649e-23
    f, T = 0.363, float(NOISE_PARAMS["T"])
    expected = 1.0 / (np.exp(h * f * 1e9 / (kB * T)) - 1.0)
    assert co.n_thermal(f) == pytest.approx(expected, rel=1e-12)
    assert co.n_thermal(5.0) < co.n_thermal(0.5)          # monotone decreasing in f


# ---------------------------------------------------------------- Manucharyan f03 diagnosis
def test_manucharyan_third_line_is_f03_not_the_two_photon_transition():
    """01_validation.ipynb section 5 resolves the f03 residual at Phi_ext = 0 by flux
    dispersion: the measured third line moves +0.046 GHz from Phi = 0 to 0.05, f03 moves the
    same way and the two-photon 0-4 transition moves the OTHER way. Sign, not tolerance."""
    def levels(flux, n=5):
        e = make_fluxonium(EJ=8.93, EC=2.39, EL=0.52, flux=flux,
                           cutoff=110).eigenvals(evals_count=n + 1)
        return e[1:] - e[0]
    a, b = levels(0.0), levels(0.05)
    measured = 11.236 - 11.190                       # digitized Figs. 3 and 4B
    assert (b[2] - a[2]) * measured > 0, "f03 no longer disperses with the measured line"
    assert (b[3] - a[3]) / 2 * measured < 0, "two-photon 0-4 no longer excluded by its sign"


def test_manucharyan_f03_is_strongly_hybridised_with_the_array_mode():
    """The residual is not a model error: at Phi_ext = 0 the bare f03 sits within 0.14-0.36 GHz
    of the junction-array mode Manucharyan's extended model puts at 10.78 GHz, with a 0-3
    coupling near 0.53 GHz. g > detuning, so the bare Hamiltonian cannot apply and the f03 row
    is excluded from the spectrum-validation claim (CLAUDE.md)."""
    W2, G2 = 10.78, 0.324
    for EJ, EC, EL in [(9.00, 2.50, 0.52), (8.93, 2.39, 0.52), (8.97, 2.47, 0.52)]:
        fx = make_fluxonium(EJ=EJ, EC=EC, EL=EL, flux=0.0, cutoff=110)
        ev = fx.eigenvals(evals_count=4)
        phi = fx.phi_operator(energy_esys=True)
        g03 = G2 * abs(phi[0, 3]) / abs(phi[0, 1])
        detuning = (ev[3] - ev[0]) - W2
        assert g03 / abs(detuning) > 1.0, f"({EJ},{EC},{EL}): no longer strongly hybridised"
    # and the 0->2 transition is parity-forbidden at Phi = 0, so the digitized f02 is a
    # theory-curve reading rather than a measured peak
    assert abs(phi[0, 2]) < 1e-9
