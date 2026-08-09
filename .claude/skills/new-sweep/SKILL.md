---
name: new-sweep
description: Run a new parameter sweep (transmon or fluxonium) end to end following this project's conventions — grid definition, CSV-before-plot, metadata JSON, unit conversion, sweet-spot reporting. Use whenever the user asks to sweep, scan, or map a parameter (EJ, EC, EL, flux, ng, a noise amplitude), add points to an existing sweep, or re-run one at different resolution.
---

# Running a parameter sweep

## Before writing any code

Establish four things, asking the user only about what you cannot infer:

1. **Which qubit** — transmon or fluxonium. This fixes the channel set (see CLAUDE.md table).
2. **Which axes and ranges**, and how many points per axis. Existing sweeps for reference:
   transmon 40×20 = 800 points over EJ ∈ [1, 80], EC ∈ [0.1, 2.0] (~2 min);
   fluxonium flux sweep through 0.5.
3. **What is held fixed** — every parameter not swept, including noise amplitudes. Default to
   `TRANSMON_DEFAULTS` / `FLUXONIUM_DEFAULTS` / `DEFAULTS` in `src/`, and say so out loud.
4. **Output names** — `data/<name>.csv` and `data/<name>_metadata.json`. Never overwrite an
   existing CSV unless the user explicitly asked to regenerate it; those files are the study's
   committed results.

Estimate the runtime before launching (points × ~0.15 s for a transmon, more for a fluxonium
at cutoff=110). If it exceeds a couple of minutes, tell the user the estimate first.

## Structure every sweep the same way

```python
# 1. bootstrap (copy verbatim from an existing notebook — chdir to root, sys.path += src/)
# 2. grid definition, as explicit np.linspace calls
# 3. loop with try/except appending a NaN row on failure — never let one bad point kill the run
# 4. df.to_csv('data/<name>.csv', index=False)      # BEFORE any plotting
# 5. metadata JSON written next to it
# 6. only then: figures, summary prints
```

The try/except-to-NaN pattern matters: scqubits raises on some parameter combinations
(diagonalization failures at extreme ratios), and losing 800 points to one exception is the
failure mode this guards against.

## Metadata JSON

Mirror the shape of `data/transmon_sweep_metadata.json` exactly:

```json
{
  "date": "YYYY-MM-DD",
  "qubit": "Transmon",
  "sweep": "EJ [1,80]x40  x  EC [0.1,2.0]x20  = 800 points",
  "fixed_params": {"ncut": 30, "Q_cap": 1e6, "A_ng": 1e-4, "A_cc": 1e-7},
  "noise_channels": ["t1_capacitive", "tphi_1_over_f_ng", "tphi_1_over_f_cc"],
  "units": "coherence times in microseconds (scqubits returns ns; multiply by 1e-3)",
  "files": ["data/transmon_full_sweep.csv"]
}
```

Get the date from the environment, not from memory.

## Reporting the result

Always report both the global optimum **and** the optimum inside the physically meaningful
region — for the transmon that is EJ/EC ≥ 20, because the global max sits in the charge-qubit
regime and quoting it as a transmon result is wrong. State the parameters at the optimum, not
just the T2 value.

Then run the dominant-channel breakdown (`dominant_channel` in `src/coherence.py`) and report
where each channel wins — that is Q2 and it is cheap once the sweep exists.

## Before you call it done

- Spot-check one point against `01_validation.ipynb`'s reproduced literature values.
- Confirm the T2 column is in microseconds and in a plausible range (µs to low ms). A value
  near 10⁵ µs usually means a missing `* 1e-3`; a suspiciously large T2 can mean a channel
  silently returned `inf` via `safe()`.
- Confirm the NaN count is small and explain it if it is not.
