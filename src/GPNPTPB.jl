module GPNPTPB
include("ysz_experiments.jl")
const TemperatureModel = ysz_experiments.ysz_model_GAS_LoMA_Temperature
const GenericModel = ysz_experiments.ysz_model_GAS_LoMA_generic
export TemperatureModel, GenericModel, simulate_eis, simulate_cv, run_experiment

"""Run the historical model driver. T is in kelvin and pO2 in the model's bar convention."""
function run_experiment(; physical_model_name="ysz_model_GAS_LoMA_Temperature", kwargs...)
    physical_model_name in ("ysz_model_GAS_LoMA_Temperature", "ysz_model_GAS_LoMA_generic") ||
        throw(ArgumentError("Only Temperature and generic are active; ARO is commented out upstream."))
    ysz_experiments.run_new(; physical_model_name, kwargs...)
end

"""Generate EIS; TC is °C, pO2 is percent O2 at a nominal total pressure of 1 bar."""
function simulate_eis(; TC=800.0, pO2=100.0, bias=0.0, timedomain=false, kwargs...)
    0 <= pO2 <= 100 || throw(ArgumentError("pO2 must be 0–100 percent"))
    TC > -273.15 || throw(ArgumentError("temperature must exceed absolute zero"))
    run_experiment(; T=TC+273.15, pO2=max(pO2, 1e-5)/100,
                   bias, EIS_IS=!timedomain, EIS_TDS=timedomain, out_df_bool=true, kwargs...)
end

"""Generate a cyclic voltammogram using the historical model driver."""
function simulate_cv(; TC=800.0, pO2=100.0, fast=false, kwargs...)
    0 <= pO2 <= 100 || throw(ArgumentError("pO2 must be 0–100 percent"))
    TC > -273.15 || throw(ArgumentError("temperature must exceed absolute zero"))
    run_experiment(; T=TC+273.15, pO2=max(pO2, 1e-5)/100,
                   voltammetry=true, fast_CV_mode=fast, out_df_bool=true, kwargs...)
end
end
