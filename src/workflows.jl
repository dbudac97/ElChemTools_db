"""Fit a series circuit to EIS. `initial` contains physical parameters, without Err.
Bounds are dictionaries by parameter name; fixed parameters retain their block identity.
L blocks use µH; RL blocks retain the historical H convention.
"""
function fit_eec(df; structure="R-L-RCPE-RCPE", initial, fixed=Dict(),
                 lower=Dict(), upper=Dict(), alpha_low=0.0, alpha_upp=1.0,
                 with_errors=false, error_type="normalized")
    if structure=="R+Rpol" && (!isempty(fixed) || !isempty(lower) || !isempty(upper))
        throw(ArgumentError("R+Rpol is an endpoint estimator, not a bounded fit"))
    end
    data=sort(deepcopy(validate_eis(df)),:f)
    eec=get_EEC(structure)
    names=filter(!=("Err"),eec.prms_names)
    length(initial)==length(names) || throw(ArgumentError("initial needs $(length(names)) values: $names"))
    eec.prms_values=Float64.(structure=="R+Rpol" ? initial : vcat(collect(initial),0.0))
    for dict in (fixed,lower,upper)
        all(k -> k in names,keys(dict)) || throw(ArgumentError("unknown parameter name"))
    end
    mask=trues(length(eec.prms_names))
    for (i,name) in enumerate(eec.prms_names)
        if haskey(fixed,name)
            eec.prms_values[i]=fixed[name];mask[i]=false
        end
        eec.lower_limits_for_fitting[i]=get(lower,name,-Inf)
        eec.upper_limits_for_fitting[i]=get(upper,name,Inf)
    end
    EEC_find_fit!(eec,data;mask,alpha_low,alpha_upp,with_errors,error_type)
    eec
end

"""Evaluate a fitted representation on positive frequencies; does not alter original data."""
function reconstruct_eis(drt::DRT_struct; frequencies=drt.EIS_df.f)
    f=Float64.(collect(frequencies))
    all(x -> isfinite(x) && x>0,f) && allunique(f) || throw(ArgumentError("invalid reconstruction frequencies"))
    # Gold models RC + Rohm only. Its NaN L means unestimated, not an inductive contribution.
    L=drt.control.alg=="Gold" ? 0.0 : drt.L
    z=[drt.R_ohm+im*2pi*x*L+sum(drt.h./(1 .+ im*2pi*x.*drt.tau_range)) for x in f]
    DataFrame(f=f,Z=ComplexF64.(z))
end
function reconstruct_eis(eec::EEC_data_struct; frequencies)
    f=Float64.(collect(frequencies))
    all(x -> isfinite(x) && x>0,f) && allunique(f) || throw(ArgumentError("invalid reconstruction frequencies"))
    get_EIS_from_EEC(eec;f_range=f)
end

"""Analyze any EIS spectrum. Metadata is opaque and does not affect the calculation."""
function analyze_eis(df; method=:drt, control=DRT_control_struct(), metadata=NamedTuple(), kwargs...)
    original=deepcopy(validate_eis(df))
    if method==:drt
        isempty(kwargs) || throw(ArgumentError("use control for DRT options"))
        result=get_DRT(original,control)
        reconstructed=reconstruct_eis(result)
    elseif method==:eec
        result=fit_eec(original;kwargs...)
        reconstructed=reconstruct_eis(result;frequencies=original.f)
    else
        throw(ArgumentError("method must be :drt or :eec"))
    end
    (;original,result,reconstructed,metadata)
end
analyze_eis(spectra::AbstractVector;kwargs...) = [analyze_eis(df;kwargs...) for df in spectra]

"""Return original, processed spectrum, removed original rows and an operation log.
All removals and reconstructions are explicit. No source DataFrame is modified.
"""
function preprocess_eis(df; drop_indices=Int[], control=EIS_preprocessing_control(),
                        smooth=:none, drt_control=DRT_control_struct(), circuit=nothing)
    original=deepcopy(validate_eis(df))
    all(i -> i isa Integer && 1<=i<=nrow(df),drop_indices) || throw(ArgumentError("invalid point indices"))
    keep=setdiff(1:nrow(original),drop_indices)
    isempty(keep) && throw(ArgumentError("cannot remove all points"))
    processed=EIS_preprocessing(original[keep,:],control)
    nrow(processed)>0 || throw(ArgumentError("preprocessing removed all points"))
    removed_indices=findall(f -> !(f in processed.f),original.f)
    removed=original[removed_indices,:]
    operations=String[]
    isempty(drop_indices) || push!(operations,"manual removal: $(sort(unique(drop_indices)))")
    control.f_interval==Nothing || push!(operations,"frequency crop: $(control.f_interval)")
    isfinite(control.outlayers_threshold) && push!(operations,"outlier filter: threshold $(control.outlayers_threshold)")
    control.scale_factor==Nothing || push!(operations,"scale: $(control.scale_factor)")
    control.subtract_inductance_fac==0 || push!(operations,"subtract fitted inductance: factor $(control.subtract_inductance_fac)")
    control.add_inductance==0 || push!(operations,"add inductance [H]: $(control.add_inductance)")
    control.use_DRT && push!(operations,"DRT reconstruction in legacy preprocessing")
    control.DAN_prepro && push!(operations,"synthetic boundary points added")
    representation=nothing
    if smooth==:drt
        representation=get_DRT(processed,drt_control)
        processed=reconstruct_eis(representation;frequencies=processed.f)
        push!(operations,"DRT reconstruction for presentation")
    elseif smooth==:eec
        circuit isa EEC_data_struct || throw(ArgumentError("provide a fitted circuit"))
        representation=deepcopy(circuit)
        processed=reconstruct_eis(representation;frequencies=processed.f)
        push!(operations,"EEC reconstruction for presentation")
    elseif smooth!=:none
        throw(ArgumentError("smooth must be :none, :drt or :eec"))
    end
    (;original,processed,removed,removed_indices,operations,representation)
end

# Small compatibility wrapper for the old DRT demo entry point.
function test_DRT(; EIS_input=nothing, mode="EEC", data_set="", EEC_structure="R-RC", EEC_prms=[1.0,2.0,0.01,0.0],
                    alg="Tikhonov",lambda=0.0,k_Gold=10000,tau_min_fac=10.0,tau_max_fac=10.0,tau_range_fac=2.0,
                    divide_R_peaks="",peak_merge_tol=0.0,plot_bool=false,f_range=(0.1,1e5,1.2))
    control=DRT_control_struct(;alg,lambda,k_Gold,tau_min_fac,tau_max_fac,tau_range_fac,divide_R_peaks,peak_merge_tol)
    df = if EIS_input !== nothing
        EIS_input
    elseif mode=="exp"
        read_eis(data_set)
    elseif mode=="EEC"
        eec=view_EEC(;EEC_structure,prms_values=EEC_prms,pyplot=false)
        get_EIS_from_EEC(eec;f_range=get_f_range_from_triplet(f_range))
    else
        throw(ArgumentError("physical simulation moved to GPNPTPB; provide EIS_input"))
    end
    result=get_DRT(df,control)
    plot_bool && plot_drt(result)
    result
end
