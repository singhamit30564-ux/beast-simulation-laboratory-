# Priority pack v2 validation — 2026-09-29

- `python -m py_compile app.py $(find beastlab -name '*.py')`: passed.
- `python -m pytest -q`: **38 passed**, about 10 seconds in this sandbox.
- AppTest: all eight render functions, normal and Lite; root navigation shell,
  persistent Lite toggle, clear session, sickle case execution, guide ranking/input
  invalidation, offline and mocked GRCh38 scan, ML inference/retraining/loss panel,
  RL benchmark and export preparation.
- Unit checks: input alphabet/FASTA/limits; canonical export hash and PDF magic;
  both-strand coordinates, PAM adjacency, overlap, mismatch/seed rules; API failure
  fallback; toy case success/failure; six SVG frame traces/controls; 3 MLP tests and
  3 environment tests; torch absent from ordinary page imports.
- Live Ensembl request from sandbox: **URLError**. Bundled phiX174 fallback returned
  5,386 bases. Successful online scan path tested with a fixture, not claimed as a
  successful live integration check.
- Single-process CPU model load + DQN/random/greedy comparison: **0.136 s** after
  module imports; peak RSS **367.7 MiB** including imports (Linux ru_maxrss).
  This is not a multiuser load test or low-end physical phone measurement.
- Default `ACGTACGTACGTACGTACGT`, seed 42, shipped MLP, fixed NGG:

  | Policy | Before (%) | After (%) |
  | --- | ---: | ---: |
  | DQN | 71.68 | 68.51 |
  | random | 71.68 | 71.81 |
  | greedy | 71.68 | 75.78 |

  Actual measured toy scores; DQN loses here. One guide is not a statistical
  benchmark, and greedy evaluates 80 candidates per step.

Remaining scope limits: window-only GRCh38 coverage, no biological confidence
calibration or target-preserving RL constraints, synthetic cure-rule demonstration,
linear phiX174 scan, no physical-device/browser animation performance testing.
CI is supplied but remote run status must be checked on the PR.
