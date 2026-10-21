abstract type abstract_simulation end

"""Validate finite EIS data with positive, unique frequencies. Does not change its order."""
function validate_eis(df)
    all(n -> n in propertynames(df), (:f, :Z)) || throw(ArgumentError("EIS needs columns f and Z"))
    nrow(df) > 0 || throw(ArgumentError("empty EIS spectrum"))
    all(f -> f isa Real && isfinite(f) && f > 0, df.f) || throw(ArgumentError("frequencies must be finite and positive"))
    all(z -> z isa Number && isfinite(z), df.Z) || throw(ArgumentError("impedances must be finite numbers"))
    allunique(df.f) || throw(ArgumentError("duplicate frequencies"))
    df
end

"""Read tabular CSV (f,ReZ,ImZ) or a ZView-style .z file; return increasing f."""
function read_eis(path::AbstractString; f_column=:f, re_column=:ReZ, im_column=:ImZ,
                  imaginary_sign=1, kwargs...)
    imaginary_sign in (-1,1) || throw(ArgumentError("imaginary_sign must be +1 or -1"))
    if endswith(lowercase(path), ".z")
        rows = Tuple{Float64,ComplexF64}[]
        started = false
        for line in eachline(path)
            if strip(line) == "End Comments"
                started = true
            elseif started && !isempty(strip(line))
                fields = split(line)
                length(fields) >= 6 || throw(ArgumentError("expected at least six columns after End Comments"))
                push!(rows, (parse(Float64,fields[1]), complex(parse(Float64,fields[5]),parse(Float64,fields[6]))))
            end
        end
        started || throw(ArgumentError("missing End Comments header"))
        df = DataFrame(f=first.(rows), Z=last.(rows))
    else
        raw = CSV.read(path, DataFrame; kwargs...)
        columns = try
            (raw[!,f_column],raw[!,re_column],raw[!,im_column])
        catch err
            err isa InterruptException && rethrow()
            throw(ArgumentError("CSV needs the specified frequency, real and imaginary columns"))
        end
        df = DataFrame(f=Float64.(columns[1]), Z=complex.(Float64.(columns[2]),imaginary_sign.*Float64.(columns[3])))
    end
    sort!(validate_eis(df), :f)
end

function write_eis(path::AbstractString, df)
    validate_eis(df)
    if endswith(lowercase(path), ".z")
        open(path, "w") do io
            println(io, "End Comments")
            for i in sortperm(df.f; rev=true)
                println(io, join((df.f[i],0,0,0,real(df.Z[i]),imag(df.Z[i])), '\t'))
            end
        end
    else
        CSV.write(path, DataFrame(f=df.f, ReZ=real.(df.Z), ImZ=imag.(df.Z)))
    end
    path
end
const _eis_importer = Ref{Any}(nothing)
"""Register a function accepting TC,pO2,bias,data_set keywords. Nothing restores file import."""
function register_eis_importer!(f)
    (f === nothing || f isa Function) || throw(ArgumentError("expected function or nothing"))
    _eis_importer[] = f
end
function import_EIStoDataFrame(;TC=0,pO2=0,bias=0,data_set="",kwargs...)
    f = _eis_importer[]
    result = f === nothing ? read_eis(data_set) : f(;TC,pO2,bias,data_set,kwargs...)
    sort!(deepcopy(validate_eis(result)), :f)
end
import_EIStoDataFrame_path(path) = read_eis(path)
export_DataFrametoEIS_path(df,path) = write_eis(path,df)
export_DataFrametoCSV(df,path) = write_eis(path,df)
Nyquist_plot(df,label="") = plot_eis(df;label)

"""Interpolate real and imaginary components linearly in f, only inside measured range."""
function resample_eis(df, frequencies)
    d = sort(deepcopy(validate_eis(df)), :f)
    f = Float64.(collect(frequencies))
    all(isfinite, f) && all(f .> 0) && allunique(f) || throw(ArgumentError("invalid target frequencies"))
    z = ComplexF64[]
    for x in f
        first(d.f) <= x <= last(d.f) || throw(ArgumentError("target frequency outside measured range"))
        i = searchsortedlast(d.f,x)
        if i == length(d.f) || d.f[i] == x
            push!(z,d.Z[i])
        else
            a=(x-d.f[i])/(d.f[i+1]-d.f[i])
            push!(z,(1-a)*d.Z[i]+a*d.Z[i+1])
        end
    end
    DataFrame(f=f,Z=z)
end
