# Provenance and extraction

Source: https://github.com/dbudac97/ElChemTools_db
Commit: e47fdf753be44a09de29b10b73360a82ecf97e16
Upstream author: Vojtěch Miloš / Masicko, MIT license (LICENSE retained).

Unmodified numerical files:
- src/models/ysz_model_GAS_LoMA_Temperature.jl
- src/models/ysz_model_GAS_LoMA_generic.jl
- src/models/YSZ_LSM_model_ARO.jl (commented prototype)
- src/ysz_experiments.jl
- prototypes/timedomain_impedance.jl

New code is only the GPNPTPB package wrapper and explicit imports for the generic support
functions. Fitting and batch orchestration are preserved under legacy for later migration.
No EEC fitting, DRT engine, EIS preprocessing or laboratory dataset selector is loaded by the
model package. Existing model bugs are explicitly deferred, rather than silently changed.
