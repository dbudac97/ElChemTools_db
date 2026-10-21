using CSV
using DataFrames

#using Revise




const EIS_standard_figure_num = 6


######################################


mutable struct EIS_simulation <: abstract_simulation
  TC::Float64
  pO2::Float64
  bias::Float64
  data_set::String
  #
  #
  f_range::Tuple
  f_interval::Tuple
  fig_size::Tuple
  fig_num::Int16
  #
  checknodes::Any
  fitness_factor::Float64
  #
  #
  use_DRT::Bool
  DRT_control::DRT_control_struct
  DRT_draw_semicircles::Bool
  DRT_backward_check::Bool
  #
  plot_option::String
  plot_legend::Bool
  #
  name::String
  ID::Int16

  EIS_simulation() = new(0.,0.,0.,"",(1.1,9500.,1.2),(-Inf,Inf),(9,6),6,
    EIS_get_checknodes_geometrical(1.1,9500.,1.2),1.,true,DRT_control_struct(),false,false,
    "Nyq Bode DRT RC leg",false,"EIS",2)
end

function string(SIM::EIS_simulation)
  return "$(SIM.name)_TC_$(SIM.TC)_pO2_$(SIM.pO2)_bias_$(SIM.bias)_data_set_$(SIM.data_set)"
end

function EIS_simulation(TC, pO2, bias=0.0; data_set="", f_range=EIS_get_shared_f_range(), f_interval=(-Inf, Inf), fig_size=(9, 6), fig_num=-1, fitness_factor=1.0, use_DRT=true, DRT_control=DRT_control_struct(), DRT_draw_semicircles=false, DRT_backward_check=false, plot_option="Nyq Bode DRT RC leg", plot_legend=false)
  output = Array{abstract_simulation}(undef,0)
  for TC_item in TC
    for pO2_item in pO2
      for bias_item in bias
        this = EIS_simulation()
        #
        this.TC = TC_item
        this.pO2 = pO2_item
        this.bias = bias_item
        this.data_set = data_set
        #
        #
        this.f_range = f_range
        this.f_interval = f_interval
        this.fig_size = fig_size
        this.fig_num = (fig_num >= 0 ? fig_num : EIS_standard_figure_num)
        #
        this.checknodes = EIS_get_checknodes_geometrical(f_range...)
        this.fitness_factor = fitness_factor
        #
        #
        this.use_DRT = use_DRT
        this.DRT_control = deepcopy(DRT_control)
        this.DRT_draw_semicircles = DRT_draw_semicircles
        this.DRT_backward_check = DRT_backward_check
        #
        this.plot_option = plot_option
        this.plot_legend = plot_legend
        #
        this.name = "EIS"
        this.ID = 2

        push!(output, this)
      end
    end
  end
  return output
end

#######################################
# potentially changable stuff #########

function EIS_get_shared_f_range()
    # experimental range is 0.1 - 65000 Hz
    # fs = (w0, w1, w_fac)
    return f_range = (1.1, 9500, 1.2)
end

function get_shared_checknodes(SIM::EIS_simulation)
    return copy(SIM.checknodes)
end

function setting_legend(SIM::EIS_simulation; latex=true)
  if latex
    return "\$\\theta=$(SIM.TC)\$°C \$\\mathrm{O}_2=$(SIM.pO2)\\% \$ \$\\mathrm{bias}=$(SIM.bias)\$ $(SIM.data_set)"
  else
    return "TC=$(SIM.TC)°C pO2=$(SIM.pO2)% bias=$(SIM.bias)V $(SIM.data_set)"
  end
end

function typical_plot_general(SIM::EIS_simulation, EIS_df, my_label, additional_string="", to_standard_figure=true; marker_style="x-", color=Nothing)
  _require_plotting()
  if to_standard_figure
    figure(SIM.fig_num, figsize=SIM.fig_size)
  end

  if  (occursin("Rf", SIM.plot_option)) ||
      (occursin("RC", SIM.plot_option)) ||
      (occursin("DRT", SIM.plot_option)) ||
      (occursin("Rtau", SIM.plot_option))
    num_rows = 300
  else
    num_rows = 200
  end

  _py(:subplots_adjust, bottom=0.07, top=0.95)

  if occursin("leg", SIM.plot_option) || SIM.plot_legend
    sLEG = subplot(num_rows + 23)
    _py(:plot, [], marker_style, label = my_label)
    _py(:axis, "off")
    legend(loc="best")
  end

  if occursin("Nyq", SIM.plot_option)
    s1 = subplot(num_rows + 21)
    title("Nyquist plot")
    xlabel("Re\$(Z) \\ [\\Omega]\$")
    ylabel("-Im\$(Z) \\ [\\Omega]\$")
    _py(:plot, real.(EIS_df.Z), -imag.(EIS_df.Z), marker_style, label = my_label)
#     if !(my_label == "") && SIM.plot_legend
#         legend(loc="best")
#     end
    grid(true)
    s1.set_aspect(1.0)
  end



  if occursin("Bode", SIM.plot_option)
    s2 = subplot(num_rows + 22)
    title("Bode plot: Re vs log f \$(\$Hz\$)\$")
    ylabel("Re\$(Z) \\ [\\Omega]\$")
    _py(:plot, log10.(EIS_df.f), real.(EIS_df.Z), marker_style, label = my_label)

    s3 = subplot(num_rows + 24)
    title("Bode plot: - Im vs log f \$(\$Hz\$)\$")
    ylabel("-Im\$(Z) \\ [\\Omega]\$")
    _py(:plot, log10.(EIS_df.f), -imag.(EIS_df.Z), marker_style, label = my_label)
  end

  if SIM.use_DRT
    DRT_actual = get_DRT(EIS_df, SIM.DRT_control)

    if occursin("DRT", SIM.plot_option)
      s5 = subplot(num_rows + 25)
      plot_DRT_h(DRT_actual, false)
    elseif occursin("Rtau", SIM.plot_option)
      s5 = subplot(num_rows + 25)
      plot_DRT_Rtau(DRT_actual, false)
    elseif occursin("Rf", SIM.plot_option)
      s5 = subplot(num_rows + 25)
      plot_DRT_Rf(DRT_actual, false)
    end

    if occursin("RC", SIM.plot_option)
      s4 = subplot(num_rows + 26)
      plot_DRT_RC(DRT_actual, false)
    end

    if SIM.DRT_draw_semicircles
      R_ohm_auxilary = DRT_actual.R_ohm
      SIM_aux = deepcopy(SIM)
      SIM_aux.use_DRT = false
      for i in 1:size(DRT_actual.peaks_df, 1)
        EIS_RC = EIS_get_RC_CPE_elements(DRT_actual.peaks_df.R[i], DRT_actual.peaks_df.C[i], 0, 1, 1, R_ohm_auxilary, f_range=SIM.f_range)
        typical_plot_general(SIM_aux, EIS_RC, "", "", true, marker_style=marker_style)
        R_ohm_auxilary += DRT_actual.peaks_df.R[i]
      end
    end

    if SIM.DRT_backward_check
      typical_plot_sim(EIS_simulation(800, 80, 0.0, use_DRT=false, fig_num=SIM.fig_num, plot_option=SIM.plot_option)..., DRT_actual.EIS_df, "! DRT_backward_check")
    end
    return DRT_actual
  else
    if occursin("DRT", SIM.plot_option)
      s5 = subplot(num_rows + 25)
      _py(:plot, [])
    elseif occursin("Rtau", SIM.plot_option)
      s5 = subplot(num_rows + 25)
      _py(:plot, [])
    elseif occursin("Rf", SIM.plot_option)
      s5 = subplot(num_rows + 25)
      _py(:plot, [])
    end

    if occursin("RC", SIM.plot_option)
      s4 = subplot(num_rows + 26)
      _py(:plot, [])
    end
    return
  end
end

function typical_plot_sim(SIM::EIS_simulation, EIS_df, additional_string="", to_standard_figure=true; color=Nothing)
  if length(additional_string)>0 && additional_string[1]=='!'
    my_label = additional_string[2:end]
  else
    my_label = "sim $(setting_legend(SIM))$(additional_string)"
  end
  typical_plot_general(SIM, EIS_df, my_label, additional_string, to_standard_figure, marker_style="x-", color=color)
end

function typical_plot_exp(SIM::EIS_simulation, EIS_df, additional_string="", to_standard_figure=true; color=Nothing)
  if length(additional_string)>0 && additional_string[1]=='!'
    my_label = additional_string[2:end]
  else
    my_label = "exp $(setting_legend(SIM))$(additional_string)"
  end
  typical_plot_general(SIM, EIS_df, my_label, additional_string, to_standard_figure, marker_style="x:", color=color)
end

function fitness_error_report(SIM::EIS_simulation, plot_prms_string, EIS_exp, EIS_sim)
  println("EIS fitness error $(setting_legend(SIM, latex=false)) $(plot_prms_string)  => ", fitnessFunction(SIM, EIS_sim, EIS_exp))
end


#######################################
# relatively static part ##############

function apply_checknodes(SIM::EIS_simulation, EIS_in, checknodes)
    DataFrame( f = checknodes, Z = EIS_get_Z_values(EIS_in, checknodes))
end

function load_file_prms(SIM::EIS_simulation; save_dir, prms, prms_names=("A0", "R0", "DGA", "DGR", "betaR", "SR"), scripted_tuple=(0, 0, 0, 0, 0, 0), throw_exception=true, verbose=true)

    (out_path, out_name) = filename_format_prms( save_dir=save_dir, prefix="EIS", prms=prms, prms_names=prms_names, scripted_tuple=scripted_tuple)

    df_out = DataFrame()
    try
      df = CSV.read(string(out_path, out_name), DataFrame)
      df_out = DataFrame(f = df.f, Z = (df.ReZ .+ df.ImZ .*im))
    catch e
      if e isa InterruptException
        rethrow(e)
      else
        if verbose
          println("file not found: $(out_path)$(out_name)")
        end
        if throw_exception
          rethrow(e)
        end
      end
    end
    return df_out
end

function save_file_prms(SIM::EIS_simulation, df_out, save_dir, prms, prms_names, scripted_tuple; mode="")
    if mode=="exp"
      prefix = "$(string(SIM))_experiment"
      scripted_tuple = Array{Int16}(undef, length(prms_names))
      scripted_tuple .= 0
      prms = []
      prms_names = []
    elseif mode=="sim"
      prefix = "$(string(SIM))_"
      scripted_tuple = Array{Int16}(undef, length(prms_names))
      scripted_tuple .= 0
    else
      prefix="EIS"
    end

    (out_path, out_name) = filename_format_prms( save_dir=save_dir, prefix=prefix, prms=prms, prms_names=prms_names, scripted_tuple=scripted_tuple)
    mkpath(out_path)

    CSV.write(
        string(out_path,out_name),
        DataFrame(f = df_out.f, ReZ = real.(df_out.Z), ImZ = imag.(df_out.Z))
    )
    return
end


function import_data_to_DataFrame(SIM::EIS_simulation)
  import_EIStoDataFrame(TC=SIM.TC, pO2=SIM.pO2, bias=SIM.bias, data_set=SIM.data_set)
end

function EIS_test_checknodes_range(f_range=EIS_get_shared_f_range())
  Nyquist_plot(apply_checknodes(EIS_simulation(), import_EIStoDataFrame(TC=800,pO2=100,bias=0.0),EIS_get_checknodes_geometrical(f_range...)))
end










function EIS_view_experimental_data(;TC=1, pO2=1, bias=1, data_set="",
                  plot_option="Nyq Bode Rtau RC", plot_legend=true, fig_num=12,
                  backward_check=true,
                  DRT_control=DRT_control_struct(), use_DRT=true,
                  use_checknodes=false, save_to_folder=Nothing,
                  EIS_preprocessing_control = EIS_preprocessing_control(
                                  f_interval=Nothing,
                                  add_inductance=0,
                                  trim_inductance=false,
                                  outlayers_threshold=Inf,
                                  use_DRT=false, DRT_control=DRT_control_struct()
                  ),
                  graph_3D_control = nothing,
                  gif_maker_control = nothing
          )
  _require_plotting()
  postprocessing_bool = (
                    (graph_3D_control !== nothing) ||
                    (gif_maker_control !== nothing)
                  )
  if postprocessing_bool
    enable_animations!()
    x_ranges_storage = []
    y_ranges_storage = []
    z_ranges_storage = []
    legend_ranges_storage = []
  end

  save_dir = save_to_folder == Nothing ? "" : String(save_to_folder)
  EIS_exp = DataFrame()
  EIS_exp_NEW = DataFrame()
  for TC_item in TC
    for pO2_item in pO2
      for data_set_item in (typeof(data_set)==String ? [data_set] : data_set)
        for bias_item in bias
          if use_checknodes
            SIM = EIS_simulation()
            checknodes =  get_shared_checknodes(SIM)
            EIS_exp = apply_checknodes(SIM, import_EIStoDataFrame(TC=TC_item, pO2=pO2_item, bias=bias_item, data_set=data_set_item), checknodes)
          else
            EIS_exp = import_EIStoDataFrame(TC=TC_item, pO2=pO2_item, bias=bias_item, data_set=data_set_item)
          end

          #typical_plot_exp(EIS_simulation(TC_item, pO2_item, bias_item, fig_num=fig_num, data_set=data_set_item, use_DRT=use_DRT, DRT_backward_check=true, DRT_control=DRT_control, plot_option=plot_option, plot_legend=plot_legend)..., EIS_exp, "")

          # NEW branch
          EIS_exp_NEW = EIS_preprocessing(EIS_exp, EIS_preprocessing_control)
          #
          loc_SIM = EIS_simulation(TC_item, pO2_item, bias_item, fig_num=fig_num, data_set=data_set_item,
                                use_DRT=use_DRT, DRT_control=DRT_control, DRT_backward_check=backward_check,
                                plot_option=plot_option, plot_legend=plot_legend
                              )[1]
          if postprocessing_bool
            push!(x_ranges_storage, deepcopy(real.(EIS_exp_NEW[!, :Z]) ))
            push!(y_ranges_storage, deepcopy((-1).*imag.(EIS_exp_NEW[!, :Z]) ))
            push!(z_ranges_storage, bias_item)
            push!(legend_ranges_storage, "TC=$(TC_item)°C pO2=$(pO2_item)% data_set=$(data_set_item)")
          else
            typical_plot_exp(loc_SIM, EIS_exp_NEW, "")
            #
            if save_to_folder!=Nothing
              save_file_prms(loc_SIM, EIS_exp_NEW, save_dir, [], [], []; mode="exp")
            end
          end
        end
      end
    end
  end
  if (graph_3D_control !== nothing)
    graph_3D_control.z_ranges = z_ranges_storage
    graph_3D_control.legend_ranges =  legend_ranges_storage
    graph_3D_control.save_dir = save_dir
    #
    _animation_backend[].plot_2D_sequence_in_3D(x_ranges_storage, y_ranges_storage, graph_3D_control=graph_3D_control
    )
  end
  if (gif_maker_control !== nothing)
    gif_maker_control.x_ranges = x_ranges_storage
    gif_maker_control.y_ranges = y_ranges_storage
    gif_maker_control.z_ranges = z_ranges_storage
    gif_maker_control.legend_ranges = legend_ranges_storage
    gif_maker_control.save_dir = save_dir
    _animation_backend[].make_GIF(gif_maker_control)
  end
  return EIS_exp_NEW
end

# function EIS_capacitance_test()
#
# end

function EIS_get_Z_values(df, checknodes)
    resample_eis(df, checknodes).Z
end

function linComb(SIM::EIS_simulation, Q1,Q2,lambda)
    if Q1[!, 1] == Q2[!, 1]
        return DataFrame(f = Q1.f, Z = Q1.Z.*(1-lambda) .+ Q2.Z.*(lambda))
    else
        println("ERROR: linComb_EIS: shape or value mismatch:\n",Q1.Z,"\n",Q2.Z)
        throw(ArgumentError("incompatible spectra or interpolation outside measured range"))
    end
end

function biliComb(SIM::EIS_simulation, Q11,Q12,Q21,Q22,x,y)
    if Q11[!, 1] == Q12[!, 1] && Q11[!, 1] == Q21[!, 1] && Q11[!, 1] ==Q22[!, 1]
        return DataFrame(f = Q11[!, 1],
             Z = Q11[!, 2].*(1-x).*(1-y) .+ Q21[!, 2].*x.*(1-y) .+ Q12[!, 2].*(1-x).*y .+ Q22[!, 2].*x.*y
            )
    else
        println("ERROR: biliComb_EIS: shape mismatch or different *.f values")
        throw(ArgumentError("incompatible spectra or interpolation outside measured range"))
    end
end

function real_Bode_get_normalization_factor(a, b, A, B=1)
  # a, A is LF value of experimental and the reference Bode, respectively.
  # b, B is HF value of experimental and the reference Bode, respectively.
  # the following formula comes from the problem:
  #
  #               min_alpha ( (alpha*a - A)^2 + (alpha*b - B)^2 )
  return (a*A + b*B)/(a*a + b*b)
end

function fitnessFunction(::EIS_simulation, sim::DataFrame, exp::DataFrame; error_type="normalized")
    validate_eis(sim); validate_eis(exp)
    sim.f == exp.f || throw(ArgumentError("frequency grids differ"))
    floor = max(maximum(abs.(exp.Z))*1e-12, eps(Float64))
    if error_type == "relative"
        scale = max.(abs.(exp.Z), floor)
        return sqrt(sum(abs2, (sim.Z-exp.Z)./scale)/nrow(exp))
    elseif error_type == "normalized"
        scale = max(maximum(abs.(exp.Z)),floor)
        return sqrt(sum(abs2, sim.Z-exp.Z)/nrow(exp))/scale
    else
        throw(ArgumentError("error_type must be normalized or relative"))
    end
end

# function monotony_test_function(SIM::EIS_simulation, exp_EIS::DataFrame, sim_EIS::DataFrame)
#         err = 0.0
#         if  exp_EIS.f == sim_EIS.f
#
#                 for row = 1:size(exp_EIS,1)
#                         err +=( (real(exp_EIS.Z[row]) - real(sim_EIS.Z[row]))^2
#                                  + (imag(exp_EIS.Z[row]) - imag(sim_EIS.Z[row]))^2)
#                 end
#         else
#                 println("ERROR: EIS_fitnesFunction: shape mismatch or different *.f values")
#                 throw(ArgumentError("incompatible spectra or interpolation outside measured range"))
#         end
#         # returns average error per checknode
#         return sqrt(err)/Float32(size(exp_EIS,1))
# end


function EIS_get_checknodes_short()
    # frequency nodes to compare
    # must be sorted upwards
    checknodes = (
        1.0e+0,
        1.0e+1,
        1.0e+2,
        1.0e+3
    )
end

function initialize_trend_tuples(SIM::EIS_simulation, EIS_ref::DataFrame)
  trend_tuples = DataFrame(prm_value=[])
  for i in 1:size(EIS_ref, 1)
    trend_tuples[!, Symbol(string(i))] = Array{Complex}(undef, 0)
  end
  return trend_tuples
end

function get_trend_tuple(SIM::EIS_simulation, EIS_ref::DataFrame, EIS_test::DataFrame)
  # plot signed deviation from referent simulation_curve
  trend_tuple = deepcopy(EIS_test.Z - EIS_ref.Z)
  return trend_tuple
#   return [
#     real(trend_tuple),
#     imag(trend_tuple),
#     abs(trend_tuple),
#     angle(trend_tuple)
#     ]
end

# function plot_trend_tuples(SIM::EIS_simulation, trend_tuples)
#
# end


