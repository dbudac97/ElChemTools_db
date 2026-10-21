# Extraction validation — 2026-10-21

Validated locally with Julia 1.10.10:

- **93 numerical/IO/provenance checks passed** in `test/runtests.jl`.
- The suite runs without loading Python, PyPlot, Plots or a physical model.
- Synthetic RC/EEC spectra check recovery of known R/C values, fixed parameters, alpha and
  user bounds, finite errors at zero crossings, and meaningful local uncertainty handling.
- DRT checks cover nonnegative reconstruction, inductance recovery, Gold iteration validation,
  endpoint peak bookkeeping, resistance conservation, Gaussian allocation and empty crops.
- Input validation, CSV/ZView round trips, configurable decimal/column conventions, explicit
  outlier removal, nonmutation of source data and retained operation records are exercised.
- The optional PyPlot extension passed `test/plotting_smoke.jl` in a separate environment.
  The exported Nyquist/Bode figure was visually inspected. PyPlot 2.11.6, Python 3.12 and
  Matplotlib 3.10.8 were used. Legacy animation rendering was not exercised.

Actual numerical command:

```sh
julia --compiled-modules=no -O0 --project=. test/runtests.jl
```

The model extraction separately passed **15 source-integrity/structure checks**. Its five
numerical source files match the upstream baseline SHA-256 hashes. Model runtime equations
were neither changed nor physically revalidated; the model repository documents deferred work.

CI is configured for package tests on Julia 1.10 and the current stable 1.x, plus an independent
optional plotting job. Local success does not imply those remote jobs have completed.
