# Numerical conventions and repair notes

## EIS and circuits

Frequencies must be positive, finite and unique; impedances must be finite. Import sorts by
increasing frequency without changing the caller's table. Resampling interpolates linearly in
frequency and rejects extrapolation. CSV imaginary values are Im(Z), not -Im(Z).

- R: Ω; C and RC capacitance: F.
- A standalone L preserves the historical µH convention. The RL block uses H.
  This difference is explicitly retained to avoid changing existing circuit values by 10^6.
- RCPE: `R/(1+(jωRC)^α)`; C defines τ=RC and is not the CPE coefficient Q.
- RCPE(tilt): `R/(1+RQ(jω)^α)`; the field historically named C is Q [S·s^α].
- R+Rpol remains an endpoint resistance estimate with a fixed illustrative relaxation time.
  It is not a general nonlinear fit or a measurement of capacitance.

The default optimizer is a scaled Nelder–Mead fit with log/logistic bound transforms, with a squared residual objective.
All fits enforce nonnegative physical values and alpha/user bounds. Legacy Err is metadata and
is excluded from optimization and covariance. The modern API requires explicit initial values;
the older batch initial estimator supports R-L-RC and R-L-RCPE[-RCPE...][-C].
Block names are not automatically permuted after fitting. Local covariance does not establish
physical uniqueness or confidence under model misspecification.

Normalized fit error is now complex RMS divided by max |Zexp|; relative error is complex RMS
weighted by pointwise |Zexp| with a finite floor. Both are finite at component zero crossings.
This replaces the historical separate real/imag weighting and sqrt(sum)/N convention, so old
success thresholds and numerical error values need recalibration.

## DRT

`Z = Rohm + jωL + Σ h[j]/(1+jωτ[j])`. Each h[j] is a discrete resistance weight [Ω].
Comparisons of densities across grids require dividing by the log-time bin widths.

Tikhonov uses nonnegative least squares on [A; λI], retaining the historical zero-order penalty
on all coefficients. Column equilibration improves solver scaling. The former undocumented
extra constraint transferring the shortest-time weight into Rohm was removed.

Gold uses only RC and Rohm; L=NaN means unestimated. Reconstruction uses no inductive term
for this algorithm. Gold requires positive iteration counts, λ=0 and nonnegative kernel/data;
inductive data must be explicitly corrected or cropped first. H=A'A is computed once, followed
by H*x per iteration, preserving the multiplicative update.

Thresholded peak tables retain final open peaks; below-threshold weights are intentionally not
part of these groups. Advanced `divide_R_peaks="valley sum"` or `"curv sum"` accounts for every
included grid weight, including the last point. Gaussian partitioning is available with `gauss`.
Peak segmentation is heuristic and does not identify physical processes.

Automatic edge-artifact cropping is now opt-in through `crop_artifacts=true`; the DRT itself is
never modified by peak extraction. Cropping may remove genuine edge contributions, so compare
against the uncropped distribution. No independent KK validity test is implemented.

Preprocessing never subtracts an unestimated Gold inductance, and filtering/reconstruction are
not default operations. `preprocess_eis` keeps the original spectrum and records removals and
operations separately. Optional manual modifications to tables outside this workflow cannot be
tracked by the package.

Legacy (CPE) plotting aliases are rejected because the stored table does not identify which
CPE convention was used. Supply an explicit expression appropriate to the circuit instead.
