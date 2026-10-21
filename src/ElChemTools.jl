module ElChemTools
using CSV, DataFrames, LinearAlgebra, NonNegLeastSquares, Optim, LsqFit, Printf, ProgressMeter, Serialization, Random
import Base: string
include("plotting.jl")
include("general_supporting_stuff.jl")
include("data.jl")
include("DRT.jl")
include("EIS_preprocessing.jl")
include("EIS_data.jl")
include("EEC_module.jl")
include("workflows.jl")
export read_eis, write_eis, validate_eis, resample_eis, register_eis_importer!
export EIS_simulation, DRT_control_struct, EIS_preprocessing_control
export get_DRT, reconstruct_eis, analyze_eis, preprocess_eis, EIS_preprocessing
export get_EEC, evaluate_EEC, get_EIS_from_EEC, fit_eec, run_EEC_fitting, view_EEC
export plot_eis, plot_drt, plot_DRT_h, plot_DRT_RC, plot_DRT_Rtau, plot_DRT_Rf
end
