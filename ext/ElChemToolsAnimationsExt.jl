module ElChemToolsAnimationsExt
using ElChemTools, Plots, LaTeXStrings
include("../src/my_GIF_maker.jl")
function __init__()
    ElChemTools._animation_backend[] = my_GIF_maker
end
end
