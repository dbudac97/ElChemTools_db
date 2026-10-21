module ElChemToolsPyPlotExt
using ElChemTools, PyPlot
function __init__()
    ElChemTools._plot_backend[] = PyPlot
end
end
