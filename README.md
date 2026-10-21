# ElChemTools

General-purpose Julia tools for **analysis and processing of electrochemical impedance spectra**.
No material, electrode chemistry or physical transport model is assumed. Input is a table with
`f` in Hz and complex `Z` in Ω. Arbitrary experiment metadata can be attached to results.

The LSM|YSZ gPNP model has been extracted into the independent **GPNPTPB** project.
The split is based on `dbudac97/ElChemTools_db@e47fdf753be44a09de29b10b73360a82ecf97e16`.
See [migration and ownership](docs/migration.md) before updating older scripts.
Original copyright and MIT license are retained.

## Install

Julia 1.10 or newer:

```julia
using Pkg
Pkg.develop(path="/path/to/ElChemTools_db")
using ElChemTools
```

Numerical analysis needs neither Python nor VoronoiFVM. Plotting is optional:
`Pkg.add("PyPlot")`, followed by `using PyPlot`, then call `plot_eis` or `plot_drt`.
For legacy animations install and load `Plots` and `LaTeXStrings`, then obtain the animation
module with `animations = ElChemTools.enable_animations!()` and construct controls through
`animations.gif_maker_control(...)`.

## Analyze EIS

```julia
using ElChemTools
spectrum = read_eis("measurement.z") # or CSV with f, ReZ, ImZ
analysis = analyze_eis(spectrum;
    metadata=(system="any electrochemical system", sample="A"),
    control=DRT_control_struct(alg="Gold", k_Gold=10000))
analysis.result.h          # discrete resistance weights [Ω], not a log-time density
analysis.result.peaks_df   # tau_c, R, C from thresholded contiguous groups
write_eis("reconstructed.csv", analysis.reconstructed)

circuit = fit_eec(spectrum; structure="R-RCPE",
    initial=[2.0, 5.0, 0.002, 0.9], fixed=Dict("R1"=>2.0),
    alpha_low=0.2, alpha_upp=1.0)
```

Circuits are series combinations of R, C, L, RC, RL, RCPE and RCPE(tilt).
`initial` in `fit_eec` contains physical parameters only, without the legacy `Err` field.
Fixed parameters keep their names; no implicit frequency-based block swapping is performed.
`with_errors=true` uses bounded nonlinear least squares and reports estimable local standard
errors; unavailable uncertainties are NaN, never copied fitted values.

## Process and reconstruct explicitly

```julia
processed = preprocess_eis(spectrum; drop_indices=[5, 19], smooth=:drt)
processed.original          # unchanged original points
processed.removed_indices   # indices into original, before any sort
processed.operations        # record of corrections and reconstructions
write_eis("presentation.csv", processed.processed)
```

Default preprocessing preserves all points and performs no reconstruction. Automated outlier
filtering must be enabled explicitly through `EIS_preprocessing_control(outlayers_threshold=...)`.
EEC reconstruction is also supported: `preprocess_eis(spectrum; smooth=:eec, circuit=circuit)`.
The reconstructed curve is a model representation, not additional measurements.

See [numerical conventions and limitations](docs/numerics.md), the runnable examples and
`julia --project=. -e 'using Pkg; Pkg.test()'`.
