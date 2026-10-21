
using Printf
using PyPlot
#using Plots
using DataFrames
#using LeastSquaresOptim
#using Optim
#using BlackBoxOptim
using LinearAlgebra
using NonNegLeastSquares
using Revise
using Serialization

import Base.string

include("../src/general_supporting_stuff.jl")
include("../src/import_experimental_data.jl")
include("../src/export_simulated_data.jl")

include("../src/DRT.jl")
include("../src/EIS_preprocessing.jl")


include("../src/ysz_experiments.jl")

include("../src/simulations/general_simulation.jl")
include("../src/simulations/CV_simulation.jl")
include("../src/simulations/EIS_simulation.jl")
include("../src/simulations/CAP_simulation.jl")

#include("../src/EIS_preprocessing_tool.jl")

#include("../src/par_study.jl")
include("../src/EEC_module.jl")


function construct_SIM_list(;TC=Nothing, pO2=Nothing, bias=0.0, data_set=Nothing, simulations=Array{String}(undef, 0),
                    fitness_factors=Nothing, physical_model_name=Nothing)
  SIM_list = []
  try
    if fitness_factors == Nothing
      aux_array = zeros(length(simulations))
      aux_array .= 1
      fitness_factors = aux_array
    end
    SIM_list = get_SIM_list_rectangle(TC, pO2, bias, data_set, simulations, fitness_factors, physical_model_name)
  catch e
    println(e)
    println("ERROR: please define TC, pO2, bias, data_set, simulations ... OR ... define SIM_fitting")
    return throw(Exception)
  end
  return SIM_list
end

function simple_run(SIM_list=Nothing; TC=Nothing, pO2=Nothing, bias=0.0, data_set=Nothing, simulations=Array{String}(undef, 0),
                    fitness_factors=Nothing, physical_model_name=Nothing, 
                        pyplot=0, use_experiment=true, prms_names=[], prms_values=[], 
                        test=false, save_files=false, save_dir="default",
                        line_color_idx = 1, evaluate_temp=false
                        )
  # here starts the true body
  if SIM_list==Nothing
    SIM_list = construct_SIM_list(TC=TC, pO2=pO2, bias=bias, data_set=data_set, simulations=simulations, 
                                    fitness_factors=fitness_factors, physical_model_name=physical_model_name)
  end
  
  save_path="../data/simple_run/"*save_dir
  
  if evaluate_temp
    plot_temp_parameters(prms_values=prms_values, prms_names=prms_names, label=save_dir, 
                          save_file=save_files, save_dir=save_path*"/", file_name="temp_prms",
                          line_color_idx=line_color_idx)  
  end
  
  if test
    test_result = 0
  end
  res = DataFrame()
  for SIM in SIM_list
    
    
    function recursive_simple_run_call(output_prms, plot_names, plot_values, active_idx)
      if active_idx > size(prms_names,1)

        ###################################################################
        # for each combination of output_prms from prms_lists perform this #######
        prms_names_in=[]
        prms_values_in=[]
        
        append!(prms_names_in, prms_names)
        append!(prms_values_in, output_prms)
      
        
        if check_equal_size(prms_names, prms_values)  
          SIM_raw = typical_run_simulation(SIM, prms_names_in, prms_values_in, pyplot)        
        end    
        if size(plot_names,1) < 1
          plot_prms_string = ""
        else
          plot_prms_string = " $(string(plot_names)) = $(plot_values)"          
        end
        SIM_sim = apply_checknodes(SIM, SIM_raw, SIM.checknodes)
        
        ## TODO !!!
        res = SIM_sim
        ###
        
        if pyplot > 0
            typical_plot_sim(SIM, SIM_sim, plot_prms_string)
        end
        if save_files
          save_file_prms(SIM, SIM_sim, save_path, plot_values, plot_names, [], mode="sim")
        end
        if use_experiment
          if test
            test_result+=fitnessFunction(SIM, SIM_sim, SIM_exp)
          else
            fitness_error_report(SIM, plot_prms_string, SIM_exp, SIM_sim)
          end
        end
        return
        ###################################################################
        ###################################################################
      
      end
      if size(prms_values[active_idx],1)>1
        for i in prms_values[active_idx]          
          recursive_simple_run_call(
            push!(deepcopy(output_prms),i),
            push!(deepcopy(plot_names),prms_names[active_idx]),
            push!(deepcopy(plot_values),i),
            active_idx + 1)
        end
      else
        recursive_simple_run_call(
          push!(deepcopy(output_prms),prms_values[active_idx][1]),
          plot_names,
          plot_values,
          active_idx + 1)
      end
    end

    # here the true body continues
    if use_experiment
      SIM_exp = apply_checknodes(SIM, import_data_to_DataFrame(SIM), SIM.checknodes)
      if (pyplot > 0)
        typical_plot_exp(SIM, SIM_exp)
      end
      if save_files
        save_file_prms(SIM, SIM_exp, save_path, [], [], [], mode="exp")
      end
    end
    if prms_names!=Nothing
      recursive_simple_run_call([], Array{String}(undef,(0)), Array{Any}(undef,(0)), 1)
    end
  end
  
  if test
    return test_result
  else
    return res
  end
end



###########################################################
###########################################################
#### Fitness_function_interpolation #######################
###########################################################
###########################################################




###########################################################
###########################################################
####     D R T     ########################################
###########################################################
###########################################################


