using Test, ElChemTools, DataFrames, LinearAlgebra, Random
const E=ElChemTools

function spectrum(f; R0=2.0,R=5.0,C=0.002,L=0.0)
    DataFrame(f=collect(f),Z=ComplexF64[R0+R/(1+im*2pi*x*R*C)+im*2pi*x*L for x in f])
end
const freq=10.0.^range(-1,5;length=61)

@testset "Independent EIS package and input contract" begin
    @test !isdefined(E,:ysz_experiments)
    @test !isdefined(E,:CV_simulation)
    @test !isdefined(E,:PyPlot)
    @test !isdefined(E,:my_GIF_maker)
    @test_throws ArgumentError validate_eis(DataFrame(f=Float64[],Z=ComplexF64[]))
    for f in ([0.0,1.0],[-1.0,1.0],[1.0,1.0],[NaN,1.0])
        @test_throws ArgumentError validate_eis(DataFrame(f=f,Z=[1im,2im]))
    end
    @test_throws ArgumentError validate_eis(DataFrame(f=[1.0],Z=[NaN+im]))
    @test_throws ArgumentError validate_eis(DataFrame(f=[1.0],ReZ=[1.0]))
    @test_throws ArgumentError E.EIS_get_checknodes_geometrical(0.,100.,1.2)
    @test_throws ArgumentError E.EIS_get_checknodes_geometrical(1.,100.,1.)
    mktempdir() do dir
        before=readdir(dir)
        cd(dir) do
            df=spectrum(freq)
            get_DRT(df,DRT_control_struct(;lambda=0.001))
            @test readdir(dir)==before
        end
    end
end

@testset "CSV / ZView IO and generic importer" begin
    d=spectrum(freq)
    mktempdir() do dir
        for extension in ("csv","z")
            path=joinpath(dir,"spectrum.$extension")
            write_eis(path,d)
            @test read_eis(path).f≈d.f
            @test read_eis(path).Z≈d.Z rtol=1e-12
        end
        path=joinpath(dir,"custom.txt")
        write(path,"Hz;real;minus_imag\n1,0;2,0;3,0\n2,0;4,0;5,0\n")
        custom=read_eis(path;f_column=1,re_column=2,im_column=3,imaginary_sign=-1,decimal=',',delim=';')
        @test custom.f==[1.,2.]
        @test custom.Z==[2-3im,4-5im]
        path=joinpath(dir,"bad.z");write(path,"1 2 3\n")
        @test_throws ArgumentError read_eis(path)
    end
    register_eis_importer!((;kwargs...)->d)
    @test E.import_EIStoDataFrame(;data_set="arbitrary system")==d
    register_eis_importer!(nothing)
    d2=DataFrame(f=[1.,3.],Z=[2.0 + 3im,6.0 + 7im])
    @test resample_eis(d2,[2.]).Z≈[4+5im]
    @test_throws ArgumentError resample_eis(d2,[0.5])
end

@testset "Equivalent circuits and bounded physical parameters" begin
    c=get_EEC("R-RC");c.prms_values=[2.,5.,0.002,0.]
    d=spectrum(freq)
    @test get_EIS_from_EEC(c;f_range=freq).Z≈d.Z
    fitted=fit_eec(d;structure="R-RC",initial=[1.,4.,0.003],with_errors=false)
    @test fitted.prms_values[1:3]≈[2.,5.,0.002] rtol=1e-3
    @test all(isnan,fitted.standard_deviations)
    @test E.fitnessFunction(EIS_simulation(),reconstruct_eis(fitted;frequencies=freq),d)<1e-5
    bounded=fit_eec(d;structure="R-RC",initial=[2.,2.5,0.002],fixed=Dict("R1"=>2.),upper=Dict("R2"=>3.))
    @test bounded.prms_values[1]==2.
    @test bounded.prms_values[2]<=3.
    @test_throws ArgumentError fit_eec(d;structure="R-RC",initial=[1.,4.,0.003],fixed=Dict("typo"=>2.))
    @test_throws ArgumentError fit_eec(d;structure="R-RC",initial=[1.,4.])
    # Confidence estimates must exclude Err and retain fixed block identities.
    noisy=deepcopy(d);Random.seed!(71);noisy.Z .+= 1e-4.*(randn(length(freq)).+im.*randn(length(freq)))
    ls=fit_eec(noisy;structure="R-RC",initial=[2.,4.,0.003],fixed=Dict("R1"=>2.),with_errors=true)
    @test ls.prms_values[2:3]≈[5.,0.002] rtol=0.01
    @test isnan(ls.standard_deviations[1])
    @test isnan(ls.standard_deviations[end])
    @test all(isfinite,ls.standard_deviations[2:3])
    # alpha bounds must constrain the default optimizer, not just the covariance path.
    c=get_EEC("R-RCPE");c.prms_values=[2.,5.,0.002,0.65,0.]
    dc=get_EIS_from_EEC(c;f_range=freq)
    a=fit_eec(dc;structure="R-RCPE",initial=[2.,5.,0.002,0.85],fixed=Dict("R1"=>2.),alpha_low=0.8,alpha_upp=0.9)
    @test 0.8<=a.prms_values[4]<=0.9
    @test isfinite(E.fitnessFunction(EIS_simulation(),DataFrame(f=[1.,2.],Z=[2.0 + 0im,2.0 + 1im]),DataFrame(f=[1.,2.],Z=[2.0 + 0im,2.0 + 1im]);error_type="relative"))
    @test_throws ArgumentError E.fitnessFunction(EIS_simulation(),d,d;error_type="bad")
end

@testset "DRT reconstruction and physical limits" begin
    d=spectrum(freq;L=2e-6)
    result=get_DRT(d,DRT_control_struct(;tau_min_fac=10.,tau_max_fac=10.,tau_range_fac=3.))
    @test all(result.h.>=0)
    @test E.fitnessFunction(EIS_simulation(),result.EIS_df,d)<0.005
    @test result.L≈2e-6 rtol=0.03
    @test sum(result.h)≈5. rtol=0.05
    @test result.R_ohm≈2. rtol=0.05
    @test reconstruct_eis(result).Z≈result.EIS_df.Z rtol=1e-10
    gold=get_DRT(spectrum(freq),DRT_control_struct(;alg="Gold",k_Gold=3000,tau_range_fac=2.))
    @test E.fitnessFunction(EIS_simulation(),gold.EIS_df,spectrum(freq))<0.03
    @test all(isfinite,reconstruct_eis(gold).Z)
    @test isnan(gold.L)
    for k in (0,-10)
        @test_throws ArgumentError get_DRT(d,DRT_control_struct(;alg="Gold",k_Gold=k))
    end
    @test_throws ArgumentError get_DRT(d,DRT_control_struct(;alg="Gold"))
    @test_throws ArgumentError get_DRT(d,DRT_control_struct(;lambda=-1))
    @test_throws ArgumentError get_DRT(d,DRT_control_struct(;tau_range_fac=0.))
    @test_throws ArgumentError get_DRT(d,DRT_control_struct(;alg="bad"))
    @test_throws ArgumentError get_DRT(d[1:2,:],DRT_control_struct())
end

@testset "DRT peak bookkeeping conserves resistance" begin
    tau=10.0.^range(-4,1;length=6)
    d=spectrum(freq)
    for h in ([0.,0.,0.,1.,2.,3.], [3.,2.,1.,0.,1.,2.],zeros(6))
        control=DRT_control_struct(;divide_R_peaks="valley sum")
        result=E.DRT_struct(d,tau,copy(h),0.,0.,DataFrame(),[],control)
        E.evaluate_RC_peaks_from_DRT(result)
        E.divide_and_evaluate_R_peaks(result)
        @test sum(last.(result.R_peak_list);init=0.0)≈sum(h)
        if h[end]>0
            @test nrow(result.peaks_df)>=1
            @test sum(result.peaks_df.R)≈sum(h)
        end
    end
    empty=E.DRT_struct(d,tau,zeros(6),0.,0.,DataFrame(),[],DRT_control_struct(;divide_R_peaks="valley sum",crop_artifacts=true))
    E.divide_and_evaluate_R_peaks(empty)
    @test isempty(empty.R_peak_list)
    # Several adjacent peaks merged in one pass: no skipped rows and resistance preserved.
    merged=get_DRT(spectrum(freq),DRT_control_struct(;peak_merge_tol=10.))
    @test nrow(merged.peaks_df)==1
    @test merged.peaks_df.R[1]>0
end

@testset "Explicit preprocessing and provenance" begin
    d=spectrum(freq);d.Z[20]+=100+100im
    original=deepcopy(d)
    noop=preprocess_eis(d)
    @test noop.processed==d
    @test isempty(noop.removed_indices)
    @test isempty(noop.operations)
    manual=preprocess_eis(d;drop_indices=[20])
    @test manual.removed_indices==[20]
    @test manual.removed.Z==d.Z[[20]]
    @test d==original
    @test manual.original==original
    @test length(manual.operations)==1
    filtered=preprocess_eis(d;control=EIS_preprocessing_control(;outlayers_threshold=1.5))
    @test 20 in filtered.removed_indices
    @test d==original
    @test_throws ArgumentError preprocess_eis(d;drop_indices=[0])
    @test_throws ArgumentError preprocess_eis(d;smooth=:bad)
    @test_throws ArgumentError preprocess_eis(d;control=EIS_preprocessing_control(;subtract_inductance_fac=1.,DRT_control=DRT_control_struct(;alg="Gold")))
    smoothed=preprocess_eis(spectrum(freq);smooth=:drt)
    @test !isempty(smoothed.operations)
    @test smoothed.original==spectrum(freq)
    @test smoothed.processed.f==freq
    c=fit_eec(spectrum(freq);structure="R-RC",initial=[2.,4.,0.003])
    ec=preprocess_eis(spectrum(freq);smooth=:eec,circuit=c)
    @test ec.processed.Z≈spectrum(freq).Z rtol=1e-3
    analyzed=analyze_eis(spectrum(freq);metadata=(system="battery",sample="A"))
    @test analyzed.metadata.system=="battery"
    @test length(analyze_eis([spectrum(freq),spectrum(freq)]))==2
    @test E.test_DRT() isa E.DRT_struct
end

@testset "Legacy batch fits and Gaussian resistance allocation" begin
    d=spectrum(freq)
    batch=E.run_EEC_fitting(;TC=[20.],pO2=[0.],bias=[0.],data_set=["generic"],
        EIS_input=d,EEC_structure="R-L-RC",fixed_prms_names=["R1","L2"],fixed_prms_values=[2.,0.],
        save_file_bool=false,plot_bool=false,with_errors=true)
    @test batch.prms_names==["R1","L2","R3","C3","Err"]
    values=vec(batch.data)
    @test values[1]==2.
    @test values[2]==0.
    @test values[3]≈5. rtol=0.01
    @test values[4]≈0.002 rtol=0.01
    tau=10.0.^range(-4,2;length=101)
    h=[exp(-((i-35)/6)^2/2)+0.5*exp(-((i-70)/5)^2/2) for i in 1:101]
    drt=E.DRT_struct(d,tau,h,0.,0.,DataFrame(),[],DRT_control_struct(;divide_R_peaks="curv gauss"))
    E.divide_and_evaluate_R_peaks(drt)
    @test length(drt.R_peak_list)==2
    @test sum(last.(drt.R_peak_list))≈sum(h) rtol=1e-10
    @test all(p -> all(isfinite,p),drt.R_peak_list)
end
