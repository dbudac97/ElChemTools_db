# Backends are supplied by Julia package extensions, never imported by numerical analysis.
const _plot_backend = Ref{Any}(nothing)
const _animation_backend = Ref{Any}(nothing)
function _require_plotting()
    _plot_backend[] === nothing && throw(ArgumentError("Install PyPlot and load it with `using PyPlot` before plotting."))
    nothing
end
function _py(name,args...;kwargs...)
    _require_plotting()
    Base.invokelatest(getproperty(_plot_backend[],name),args...;kwargs...)
end
for name in (:figure, :subplot, :plot, :title, :xlabel, :ylabel, :xlim, :ylim,
             :legend, :grid, :clf, :gcf, :pause, :axis, :close)
    @eval $name(args...;kwargs...) = _py($(QuoteNode(name)),args...;kwargs...)
end
"""Return the optional animation module after loading Plots and LaTeXStrings."""
function enable_animations!()
    _animation_backend[] === nothing && throw(ArgumentError("Install and load Plots and LaTeXStrings before enabling animations."))
    _animation_backend[]
end
function plot_eis(df; label="", kwargs...)
    validate_eis(df)
    plot(real.(df.Z), -imag.(df.Z); label, kwargs...)
    xlabel("Re Z [Ω]"); ylabel("−Im Z [Ω]")
end
plot_drt(drt; kwargs...) = (plot(log10.(drt.tau_range), drt.h; kwargs...); xlabel("log10 τ [s]"); ylabel("resistance weight [Ω]"))
