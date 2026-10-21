using ElChemTools, PyPlot
const DataFrame=ElChemTools.DataFrame
f=10.0.^range(-1,5;length=41)
d=DataFrame(f=f,Z=2 .+ 5 ./ (1 .+ im*2pi.*f.*0.01))
ElChemTools.plot_eis(d)
plot_drt(get_DRT(d,DRT_control_struct(;lambda=0.001)))
ElChemTools.typical_plot_exp(EIS_simulation(20,0;use_DRT=false,plot_option="Bode Nyq")[1],d)
ElChemTools._py(:savefig,joinpath(ENV["ELCHEM_PLOT_OUTPUT"],"eis.png"))
println("Optional plotting smoke test passed")
