"""Convergence checks for the lab-frame pulse simulation of 04_analysis.ipynb section 3.1.

    python results/pulse_convergence.py           (repo root, conda env `qubits`)

The simulation functions are copied verbatim from build_nb04.py (GATE_LEVELS, _envelope,
gate_error, calibrate_gate), with the truncation, ODE tolerances and optimiser starts exposed
as parameters. Re-calibrates the gate at the durations that decide T_gate (transmon 6, 8 ns;
fluxonium 2.5, 3 ns) and checks whether the < 1e-3 verdict changes with
  (i) 6 / 8 / 10 levels, (ii) more optimiser starts and evaluations, (iii) tighter ODE tolerances.
"""
import sys, time, json
sys.path.insert(0, "src")
import numpy as np
import qutip as qt
from scipy.integrate import quad
from scipy.optimize import minimize, minimize_scalar
from qubits import make_transmon, make_fluxonium

ERR_TARGET = 1e-3
TOL = {"atol": 1e-10, "rtol": 1e-8}


def drive_system(qubit, nlev):
    qubit.truncated_dim = nlev
    esys = qubit.eigensys(evals_count=nlev)
    f = esys[0] - esys[0][0]
    n = np.array(qubit.n_operator(energy_esys=esys))[:nlev, :nlev]
    return f, n - np.diag(np.diag(n))


def envelope(t, T):
    sig = T / 4.0
    g = np.exp(-((t - T / 2) ** 2) / (2 * sig ** 2))
    g0 = np.exp(-((T / 2) ** 2) / (2 * sig ** 2))
    return (g - g0) / (1 - g0)


def gate_error(f, n, T, amp, drag, det, nlev, tol=TOL):
    wd = f[1] + det
    alpha = (f[2] - f[1]) - f[1]
    H0 = qt.Qobj(2 * np.pi * np.diag(f)); N = qt.Qobj(2 * np.pi * n)
    sig = T / 4.0
    g0 = np.exp(-((T / 2) ** 2) / (2 * sig ** 2))
    def eps(t, *a, **k):
        g = np.exp(-((t - T / 2) ** 2) / (2 * sig ** 2))
        I = amp * (g - g0) / (1 - g0)
        dI = amp * (-(t - T / 2) / sig ** 2) * g / (1 - g0)
        Q = drag * dI / (2 * np.pi * alpha)
        return I * np.cos(2 * np.pi * wd * t) + Q * np.sin(2 * np.pi * wd * t)
    U = qt.sesolve([H0, [N, eps]], qt.qeye(nlev), [0, T],
                   options={**tol, "nsteps": 10**7}).states[-1].full()
    R = np.diag(np.exp(1j * 2 * np.pi * wd * np.arange(nlev) * T))
    Uq = (R @ U)[:2, :2]
    Xt = np.array([[0, 1], [1, 0]], dtype=complex)
    def err(th):
        V = np.diag([1, np.exp(1j * th)]) @ Uq
        F = (abs(np.trace(Xt.conj().T @ V)) ** 2 + np.trace(V.conj().T @ V).real) / 6
        return 1 - F.real
    return minimize_scalar(err, bounds=(-np.pi, np.pi), method="bounded").fun


def calibrate(f, n, T, nlev, starts, maxfev, tol=TOL):
    area = quad(lambda t: envelope(t, T), 0, T)[0]
    a0 = 1.0 / (2 * abs(n[0, 1]) * area)
    best = None
    for (sa, d0, det0) in starts:
        r = minimize(lambda x: np.log10(max(gate_error(f, n, T, x[0], x[1], x[2], nlev, tol), 1e-14)),
                     np.array([a0 * sa, d0, det0]), method="Nelder-Mead",
                     options={"xatol": 1e-6, "fatol": 1e-3, "maxfev": maxfev})
        if best is None or r.fun < best.fun:
            best = r
    return 10 ** best.fun, best.x


QUBITS = {"transmon": make_transmon(EJ=13.0, EC=0.26, ng=0.0),
          "fluxonium": make_fluxonium(EJ=8.9, EC=2.5, EL=0.5, flux=0.5)}
DUR = {"transmon": [6.0, 8.0], "fluxonium": [2.5, 3.0]}
ORIG = [(1.0, 0.0, 0.0), (1.0, -1.0, 0.0), (1.0, 1.0, 0.0)]            # the notebook's starts
MORE = ORIG + [(0.9, -0.5, 0.01), (1.1, -0.5, -0.01), (1.0, -2.0, -0.03), (1.0, -20.0, -0.03)]

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "full"
    if mode == "time":
        f, n = drive_system(QUBITS["fluxonium"], 6)
        t0 = time.time(); e = gate_error(f, n, 3.0, 7.2856, -19.486, -0.0301, 6)
        print(f"one fluxonium evaluation: {time.time()-t0:.2f} s, error {e:.2e} (notebook 7.7e-07)")
        f, n = drive_system(QUBITS["fluxonium"], 10)
        t0 = time.time(); e = gate_error(f, n, 3.0, 7.2856, -19.486, -0.0301, 10)
        print(f"same pulse, 10 levels: {time.time()-t0:.2f} s, error {e:.2e}")
        sys.exit()
    rows = []
    cases = [("6 levels, notebook starts", 6, ORIG, 140, TOL),
             ("6 levels, 7 starts, 400 evals", 6, MORE, 400, TOL),
             ("6 levels, tight ODE tolerance", 6, ORIG, 140, {"atol": 1e-12, "rtol": 1e-10}),
             ("8 levels", 8, ORIG, 140, TOL),
             ("10 levels", 10, ORIG, 140, TOL),
             ("10 levels, 7 starts, 400 evals", 10, MORE, 400, TOL)]
    for lab, q in QUBITS.items():
        for name, nlev, starts, maxfev, tol in cases:
            f, n = drive_system(q, nlev)
            for T in DUR[lab]:
                t0 = time.time()
                e, x = calibrate(f, n, T, nlev, starts, maxfev, tol)
                rows.append(dict(arch=lab, case=name, T_ns=T, error=float(e), amp=float(x[0]),
                                 drag=float(x[1]), det_GHz=float(x[2]), passes=bool(e < ERR_TARGET)))
                print(f"{lab:9s} {name:32s} T={T:4.1f} ns  error={e:.2e}  "
                      f"{'PASS' if e < ERR_TARGET else 'fail'}  amp={x[0]:.3f}  ({time.time()-t0:.0f} s)",
                      flush=True)
    json.dump(rows, open("results/pulse_convergence.json", "w"), indent=1)
