# GPNPTPB — gPNP model of the LSM|YSZ interface

Independent Julia project extracted from ElChemTools. It contains the physical YSZ transport,
space-charge and surface reaction models and the EIS/CV numerical experiment driver.
Generic impedance analysis and preprocessing belong in **ElChemTools**.

The numerical model equations and time-domain solver are preserved byte for byte from
`dbudac97/ElChemTools_db@e47fdf753be44a09de29b10b73360a82ecf97e16`.
They have not been repaired or physically revalidated as part of the extraction.
Original MIT license and authorship are retained. This is a model-development baseline.

## Entry points

```julia
using GPNPTPB
# These can be expensive; set the physical parameters for your intended experiment first.
# eis = simulate_eis(TC=800, pO2=100, bias=0.0,
#     prms_names_in=["A.r_C"], prms_values_in=[24.0])
# cv = simulate_cv(TC=800, pO2=100, voltrate=0.01)
```

`simulate_eis` and `simulate_cv` take temperature in °C and nominal O2 fraction in percent at
1 bar total pressure. As in the historical wrapper, nominal zero is regularized to 10⁻⁷ bar.
`run_experiment` exposes the historical driver: T is kelvin and pO2 follows its bar convention.
The active model modules are exported as `TemperatureModel` and `GenericModel`.

## Status

- Active: `ysz_model_GAS_LoMA_Temperature`, `ysz_model_GAS_LoMA_generic`.
- The ARO file is a fully commented historical design, not an active module.
- EIS, CV, steady I-V and analytical capacitance code are retained. Numerical correctness
  and compatibility with current dependencies still need dedicated model validation.
- STEP / chronoamperometry is not implemented; the old driver contains a placeholder.
- Generic-model potential coupling still uses a hard-coded mesh node. This needs repair
  before mesh or geometry studies.
- Legacy fitting, CV/CAP descriptors and parameter-study scripts are retained in `legacy/`.
  They are not imported by GPNPTPB and still require migration to the independent APIs.
- `legacy/Project.toml` and `legacy/Manifest.toml` preserve the original combined environment
  as provenance, not as a second active ElChemTools installation.

There is no dependency on ElChemTools. Export generated EIS as f,ReZ,ImZ CSV and analyze it
there. Model fitting across EIS/CV and future chronoamperometry remains work for this project.
