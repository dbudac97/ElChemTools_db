# Separation from the original combined toolkit

Baseline: https://github.com/dbudac97/ElChemTools_db/commit/e47fdf753be44a09de29b10b73360a82ecf97e16

The existing repository keeps its package identity and UUID (`ElChemTools`, version 0.1).
The physical model has its own package identity (`GPNPTPB`, a different UUID).
Neither numerical core imports the other. Original history remains available on both split branches.

| Original area | Destination |
|---|---|
| EEC, DRT, peak analysis, EIS preprocessing | ElChemTools/src |
| Generic EIS plotting, batch circuit fits and parameter exports | ElChemTools/src |
| Dataset-specific file naming, CV import and physical-model calibration scripts | GPNPTPB/legacy (preserved, not loaded) |
| YSZ material models, `ysz_experiments`, time-domain impedance | GPNPTPB/src and prototypes |
| Physical model CV/CAP descriptors, fitting and parameter-study code | GPNPTPB/legacy (preserved, not loaded) |
| Physical-model documentation | GPNPTPB/txt |

## Script migration

- Use `read_eis(path)` / `write_eis(path, df)` for generic ZView or CSV IO.
  For PROKOP-style or other custom tables, specify `f_column`, `re_column`, `im_column`
  as names or column indices, `imaginary_sign=-1` for exported -Im(Z), and CSV options such
  as `decimal=','` / `delim=';'`. Column positions are no longer inferred from a filename.
- Dataset labels such as MONO_110 no longer imply a laboratory directory tree. Pass a file path
  or explicitly register a keyword-based loader using `register_eis_importer!`.
- No `custom_import_experimental_data.jl` is created or executed at package load.
- `EIS_simulation(TC,pO2,bias;...)` remains as a compatibility descriptor for EIS plots and batch
  metadata; model, mesh and time-domain simulation keywords have moved out of this package.
- `get_DRT`, `view_EEC`, `run_EEC_fitting` and existing peak plots remain available.
- Prefer `analyze_eis`, `fit_eec`, `preprocess_eis` and `reconstruct_eis` for new scripts.
- The old `simple_run`, `CV_simulation`, `CAP_simulation`, physical fitting and cluster drivers
  are not part of ElChemTools. The model repository retains their source for later migration.
- `test_DRT` is a small compatibility demo returning a DRT result, rather than an unconditional
  Boolean. It accepts experimental input or an EEC; physical mode was removed.
- Models can export ordinary `f,ReZ,ImZ` CSV files for analysis here. No circular dependency or
  separate shared repository is needed.

This is a breaking 0.1 extraction. It does not silently claim compatibility with old model scripts.
The historical material-specific entry point is intentionally not an analytical package API.
