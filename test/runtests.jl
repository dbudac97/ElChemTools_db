using Test, GPNPTPB
@testset "Independent physical-model package" begin
    @test nameof(GPNPTPB.TemperatureModel)==:ysz_model_GAS_LoMA_Temperature
    @test nameof(GPNPTPB.GenericModel)==:ysz_model_GAS_LoMA_generic
    @test !isdefined(GPNPTPB,:ElChemTools)
    @test_throws ArgumentError run_experiment(physical_model_name="YSZ_LSM_model_ARO")
    @test_throws ArgumentError simulate_eis(pO2=-1)
    @test_throws ArgumentError simulate_cv(TC=-300)
end
