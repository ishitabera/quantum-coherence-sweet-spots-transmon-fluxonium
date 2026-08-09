"""Qubit builders for the coherence sweet-spot study.

Thin wrappers around scqubits constructors so every notebook builds
transmons and fluxoniums with identical, documented defaults.
"""
import scqubits as scq

scq.settings.T1_DEFAULT_WARNING = False

# Nominal designs used throughout the study
TRANSMON_DEFAULTS  = dict(EJ=15.0, EC=0.3, ng=0.0, ncut=30)
FLUXONIUM_DEFAULTS = dict(EJ=8.9, EC=2.5, EL=0.5, flux=0.5, cutoff=110)


def make_transmon(EJ=15.0, EC=0.3, ng=0.0, ncut=30):
    """Cooper-pair-box transmon. ng=0 is the charge sweet spot."""
    return scq.Transmon(EJ=EJ, EC=EC, ng=ng, ncut=ncut)


def make_fluxonium(EJ=8.9, EC=2.5, EL=0.5, flux=0.5, cutoff=110):
    """Fluxonium. flux=0.5 (half flux quantum) is the flux sweet spot."""
    return scq.Fluxonium(EJ=EJ, EC=EC, EL=EL, flux=flux, cutoff=cutoff)


def qubit_frequency(qubit):
    """0->1 transition frequency in GHz."""
    ev = qubit.eigenvals(evals_count=2)
    return ev[1] - ev[0]


def anharmonicity(qubit):
    """(E2-E1)-(E1-E0) in GHz. Negative for a transmon."""
    ev = qubit.eigenvals(evals_count=3)
    return (ev[2] - ev[1]) - (ev[1] - ev[0])
