using Printf

#using Plots
using DataFrames
#using LeastSquaresOptim

using LinearAlgebra


#include("../src/simulations/EIS_simulation.jl")


using Random
using Base: @kwdef

const DRT_standard_figure = 33

@kwdef mutable struct DRT_control_struct
  lambda::Float64=0.0
  tau_min_fac::Float64=1.0
  tau_max_fac::Float64=1.0
  tau_range_fac::Float64=1.0
  peak_merge_tol::Float64=0.0
  specified_f_range = Nothing
  k_Gold::Int64=10000
  alg::String="Tikhonov"
  divide_R_peaks::String=""
  R_peaks_plot_bool::Bool=false
  export_in_f::Bool=false
  crop_artifacts::Bool=false
end

mutable struct DRT_struct
  EIS_df::DataFrame
  tau_range::Array{Float64}
  h::Array{Float64}
  R_ohm::Float64
  L::Float64
  peaks_df::DataFrame
  R_peak_list::Array{Any}
  control::DRT_control_struct
end

include("peak_analysis.jl")

# 
# function perform_KK_test(EIS_df::DataFrame, DRT::DRT_struct)
#   EIS_copy = deepcopy(EIS_df)
#   
#   L = DRT.L
#   # subtract inductance
#   for (i, f) in enumerate(EIS_copy.f)
#     EIS_copy.Z[i] += - im*2*pi*f*L
#   end
#   
# end

function evaluate_RC_peaks_from_DRT(DRT::DRT_struct)
  
  h_max = maximum(DRT.h)
  threshold = h_max/50.0
  
  function peak_assesement(temp_peak_df)
    R = 0.0
    tau_c = 0.0
    for i in 1:size(temp_peak_df, 1)
      R += temp_peak_df.h[i]
      # only an intermediate step ... weighted average of peak points
      tau_c += temp_peak_df.h[i] * temp_peak_df.tau[i]
    end
    tau_c = tau_c/R
    C = tau_c/R
    return tau_c, R, C
  end
  
  #peak finding
  peaks = DataFrame(tau_c = [], R = [], C = [])
  
  temp_peak_df = DataFrame(tau = [], h = [])
  in_peak_bool = false
  for (i, tau) in enumerate(DRT.tau_range)
    if DRT.h[i] > threshold
      push!(temp_peak_df, (tau, DRT.h[i]))
    else
      if size(temp_peak_df, 1) > 0
        push!(peaks, peak_assesement(temp_peak_df))
        
        temp_peak_df = DataFrame(tau = [], h = [])
      end
    end
  end
  
  if nrow(temp_peak_df) > 0
    push!(peaks, peak_assesement(temp_peak_df))
  end
  DRT.peaks_df = peaks
  return
end

# function plot_DRT(DRT::DRT_struct, to_standard_figure=true)
#   if to_standard_figure
#     figure(DRT_standard_figure)
#   end
#   title("DRT")
#   plot(log10.(DRT.tau_range), DRT.h, "-x", label="l=$(DRT.lambda)")
#   xlabel("log10(\$\\tau\$ [s])")
#   ylabel("\$h(\\tau)\$ [Ohm]")
#   if to_standard_figure
#     legend()
#   end
#   peaks_df = get_RC_df_from_DRT(DRT.tau_range, DRT.h)
#   println("DRT_parameters:  R_ohm = $(DRT.R_ohm)  L = $(DRT.L)")
#   for i in 1:size(peaks_df,1)
#     println(">> R$i = $(peaks_df.R[i])   C$i = $(peaks_df.C[i])    ... (tau_c$i = $(peaks_df.tau_c[i]))")
#   end
# end



function plot_DRT_h(DRT::DRT_struct, to_standard_figure=true, print_bool=false, plot_lambda=true)
  if to_standard_figure
    figure(DRT_standard_figure)
  end
  
  title("DRT")
  plot(log10.(DRT.tau_range), DRT.h, "-x", label="\$\\lambda\$=$(DRT.control.lambda)")
  xlabel("log10(\$\\tau\$ [s])")
  ylabel("\$h(\\tau)\$ [Ohm]")
  if (to_standard_figure || plot_lambda)
    legend()
  end
  if print_bool
    println("non-DRT_parameters:  R_ohm = $(DRT.R_ohm)  L = $(DRT.L)")  
  end
end
  
function plot_DRT_RC(DRT::DRT_struct, to_standard_figure=true, print_bool=true)
  if to_standard_figure
    figure(DRT_standard_figure)
  end
  peaks_df = DRT.peaks_df
#   
  title("RC diagram")
  xlabel("C [F]")
  ylabel("R [Ohm]")
  plot(peaks_df.C, peaks_df.R, "x")
  grid(true)
  
  if print_bool
    println("non-DRT_parameters:  R_ohm = $(DRT.R_ohm)  L = $(DRT.L)") 
    for i in 1:size(peaks_df,1)
      println(">> f$i = $(1/(2*pi*peaks_df.C[i]*peaks_df.R[i]))    R$i = $(peaks_df.R[i])   C$i = $(peaks_df.C[i])   ... (tau_c$i = $(peaks_df.tau_c[i]))")
    end
  end
end

function plot_DRT_Rtau(DRT::DRT_struct, to_standard_figure=true, print_bool=false)
  if to_standard_figure
    figure(DRT_standard_figure)
  end
  peaks_df = DRT.peaks_df
  
  title("Rtau diagram")
  xlabel("log10(\$\\tau\$ [s])")
  ylabel("R [Ohm]")
  xlim(log10.([DRT.tau_range[1], DRT.tau_range[end]])...)
  
  
  previos_maximum = (_py(:gca)).get_ylim()[2]
  
  ylim(0, max(
    (length(peaks_df.R) > 0 ? maximum(peaks_df.R)*1.1 : 0.0001), 
    previos_maximum)  )
  plot(log10.(peaks_df.C.*peaks_df.R), peaks_df.R, "o")
  grid(true)

end

function plot_DRT_Rf(DRT::DRT_struct, to_standard_figure=true, print_bool=false)
  if to_standard_figure
    figure(DRT_standard_figure)
  end
  peaks_df = DRT.peaks_df
  
  title("Rf diagram")
  xlabel("log10(\$f\$ [Hz])")
  ylabel("R [Ohm]")
  xlim(log10.(1 ./(2*pi*[DRT.tau_range[end], DRT.tau_range[1]]))...)
  plot(log10.(1 ./(2*pi*peaks_df.C.*peaks_df.R)), peaks_df.R, "o")
  grid(true)
  
  if print_bool
    for i in 1:size(peaks_df,1)
      println(">> f$i = $(1/(2*pi*peaks_df.C[i]*peaks_df.R[i]))    R$i = $(peaks_df.R[i])    ... (C$i = $(peaks_df.C[i]))")
    end
  end
end


function get_expspace(A,B,n)
  isfinite(A) && isfinite(B) && 0 < A < B || throw(ArgumentError("invalid positive time range"))
  n >= 2 || throw(ArgumentError("at least two relaxation nodes required"))
  exp.(range(log(A),log(B);length=round(Int,n)))
end

function construct_A_matrix_and_b(f_nodes, tau_range, Z, lambda, control)
  N_f = size(f_nodes, 1)
  N_tau = size(tau_range, 1)
  
  # h, R_ohm, L
  if control.alg == "Gold"
    n_cols = N_tau + 1  
  else
    n_cols = N_tau + 2  
  end
  
  # real(Z), imag(Z), regularization
  if lambda == 0.0
    n_rows = 2*N_f
  else
    n_rows = 2*N_f + n_cols
  end
#   @show tau_range
#   @show f_nodes
    
  A = Matrix{Float64}(undef, n_rows, n_cols)
  b = Vector{Float64}(undef, n_rows)
  
  im_fac = -1.0

  # assemble A and b
  for (i, f) in enumerate(f_nodes)
    #RC
    for (j, tau) in enumerate(tau_range)
      A[i, j]       = real(1/(1 + im*(2*pi*f)*tau))
      A[N_f + i, j] = im_fac*imag(1/(1 + im*(2*pi*f)*tau))
    end
    # R_ohm
    A[i, N_tau + 1] = 1
    A[N_f + i, N_tau + 1] = im_fac*0
    # L
    if control.alg != "Gold"
      A[i, N_tau + 2] = 0
      A[N_f + i, N_tau + 2] = im_fac*2*pi*f
    end
    
    # b
    if Z != Nothing
      b[i] = real(Z[i])
      b[N_f + i] = im_fac*imag(Z[i])
    end
  end
  
  if lambda != 0.0
    # assemble "regularization" part of A and b
    A[2*N_f + 1 : end, :] .= Diagonal([lambda for i in 1:n_cols])
    b[2*N_f + 1 : end] .= 0
  end
  
#   @show A[N_f,:]
#   @show b[N_f]
#   @show 2222222
#   @show A[2*N_f,:]
#   @show b[2*N_f]

  return (A, b, N_f, N_tau)
end


function Gold_solve(A,b,k)
  k > 0 || throw(ArgumentError("Gold iterations must be positive"))
  all(A .>= 0) && all(b .>= 0) || throw(ArgumentError("Gold requires nonnegative kernel and data; remove inductive data first"))
  H=transpose(A)*A
  q=transpose(A)*b
  x=ones(size(A,2))
  denom=similar(x)
  for _ in 1:k
    mul!(denom,H,x)
    for i in eachindex(x)
      x[i] = denom[i] > 0 ? x[i]*q[i]/denom[i] : 0.0
    end
  end
  x
end

function _scaled_nnls(A,b)
  scales=max.(sqrt.(vec(sum(abs2,A;dims=1))),eps(Float64))
  nonneg_lsq(A ./ transpose(scales),b;alg=:nnls)./scales
end

function get_DRT(EIS_df::DataFrame, control::DRT_control_struct, debug_mode=false)
  validate_eis(EIS_df)
  nrow(EIS_df) >= 3 || throw(ArgumentError("DRT needs at least three frequencies"))
  control.alg in ("Tikhonov","Gold") || throw(ArgumentError("unknown DRT algorithm"))
  all(x -> isfinite(x) && x > 0, (control.tau_min_fac,control.tau_max_fac,control.tau_range_fac)) || throw(ArgumentError("invalid DRT grid factors"))
  isfinite(control.lambda) && control.lambda >= 0 || throw(ArgumentError("lambda must be nonnegative"))
  if control.alg == "Gold"
    control.k_Gold > 0 || throw(ArgumentError("Gold iterations must be positive"))
    control.lambda == 0 || throw(ArgumentError("Gold does not accept lambda"))
  end
  #println(control.lambda)
  #@show tau_max_fac
  tau_min = 1.0/(2*pi*maximum(EIS_df.f)) / control.tau_min_fac
  tau_max = 1.0/(2*pi*minimum(EIS_df.f)) * control.tau_max_fac
  tau_range = get_expspace(tau_min, tau_max, control.tau_range_fac*size(EIS_df.f, 1)-2)
  #tau_range = [0.001]
  
  (A, b, N_f, N_tau) = construct_A_matrix_and_b(EIS_df.f, tau_range, EIS_df.Z, control.lambda, control)
  
 # @show A b
  if control.alg == "Tikhonov"
    solution= _scaled_nnls(A,b)
  elseif control.alg == "Gold"
    if control.lambda != 0.0
      println("ERROR: when using alg=\"Gold\", lambda must be set to 0.0")
      return throw(ArgumentError("invalid input; see preceding message"))
    end
    solution = Gold_solve(A, b, control.k_Gold)
  else
    println("ERROR: DRT >>> wrong algoritm chosen: $(control.alg)")
    throw(ArgumentError("invalid input; see preceding message"))
  end
  
  # reconstructin EIS from DRT
  # getting rid of the corner effects
  # if control.lambda == 0.0 && false
  #   solution[end-7 : end-2] .= 0
  #   solution[1 : 3] .=0
  # end
  
  if control.specified_f_range != Nothing  
    f_nodes = EIS_get_checknodes_geometrical(control.specified_f_range...)
        
    # assemble new A w.r.t. specified_f_range
    (A, b, N_f, N_tau) = construct_A_matrix_and_b(f_nodes, tau_range, Nothing, 0.0, control)
  
    EIS_post = A*solution
    
    Z_specified = Array{Complex}(undef, N_f)
    for i in 1:N_f
      Z_specified[i] = Complex(EIS_post[i], -EIS_post[N_f + i])
    end
    
    EIS_new = DataFrame(f = f_nodes, Z = Z_specified)
  else
    EIS_post = A*solution
    
    EIS_new = deepcopy(EIS_df)
    for i in 1:N_f
      EIS_new.Z[i] = EIS_post[i] - im*EIS_post[N_f + i]
    end
  end

  if control.alg == "Gold"
    L_out = NaN
  else
    L_out = solution[N_tau+2]
  end
  
  DRT_out = DRT_struct(
    EIS_new, tau_range, solution[1:N_tau], solution[N_tau+1], L_out, DataFrame(), [], control)
  
  evaluate_RC_peaks_from_DRT(DRT_out)
  (control.divide_R_peaks != "") && divide_and_evaluate_R_peaks(DRT_out)

  # the following do NOT change DRT function, only changes peaks_df
  if control.peak_merge_tol > 0.0
    i = 2
    peaks = DRT_out.peaks_df
    while i <= nrow(peaks)
      if abs(log10(peaks.tau_c[i]/peaks.tau_c[i-1])) < control.peak_merge_tol
        R = peaks.R[i-1]+peaks.R[i]
        tau = (peaks.R[i-1]*peaks.tau_c[i-1]+peaks.R[i]*peaks.tau_c[i])/R
        peaks.R[i-1]=R; peaks.tau_c[i-1]=tau; peaks.C[i-1]=tau/R
        deleteat!(peaks,i)
      else
        i += 1
      end
    end
  end
  
  return DRT_out
end



function EIS_get_RC_CPE_elements(R1, C1, R2, C2, alpha, Rohm=0; f_range=Nothing)
  EIS_RC = DataFrame( f = [], Z = [])
  if f_range==Nothing
    checknodes = EIS_get_checknodes_geometrical(EIS_get_shared_f_range()...)
  else
    checknodes = EIS_get_checknodes_geometrical(f_range...)
  end
  for f in checknodes
    push!(
      EIS_RC, 
      (f, Rohm + R1/(1 + im*2*pi*f*R1*C1) + R2/(1 + ((im*2*pi*f*R2*C2)^alpha))  )
    ) 
  end
  return EIS_RC  
end


