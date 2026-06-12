# Hickit fork for phase3

This submodule keeps the Hickit core plus the cleaned P9016 blind-diploid softall baseline.

The active P9016 entry point is:

```text
run_blind_p9016_minimal.c
```

The active baseline mode is:

```text
raw_expected_soft_all
```

In this mode, each raw contact contributes to all four copy-state edges with the E-step posterior weights. The runner rejects old experiment switches for hard phase locks, top-k thinning, pcut/oracle modes, trans-specific force scaling, outlier models, and annealing schedules.

## Build

Run from this directory inside the repository-required `analysis` Conda environment:

```bash
make run_blind_p9016_minimal.bin audit_blind_p9016_full_cpu_output.bin
make smoke_blind_p9016_minimal
```

## Files kept for the active baseline

```text
run_blind_p9016_minimal.c
audit_blind_p9016_full_cpu_output.c
blind.c
hickit.h
test_blind_*.c
testdata/
```

The Hickit command-line core (`main.c`, `io.c`, `pair.c`, `bin.c`, `fdg.c`, and related headers) is retained for compile compatibility. Some internal enum values and helper routines remain in `blind.c`/`hickit.h`; they are not active P9016 baseline entry points.
