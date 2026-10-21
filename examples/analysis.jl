using ElChemTools, DataFrames
f=10.0.^range(-1,5;length=61)
data=DataFrame(f=f,Z=2 .+ 5 ./ (1 .+ im*2pi.*f.*0.01))
result=analyze_eis(data;control=DRT_control_struct(lambda=0.001))
circuit=fit_eec(data;structure="R-RC",initial=[1.,4.,0.003])
@show circuit.prms_values result.result.peaks_df
processed=preprocess_eis(data;drop_indices=[5],smooth=:drt)
@show processed.removed_indices processed.operations
