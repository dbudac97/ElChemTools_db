using Test, TOML, SHA
root=dirname(@__DIR__)
@testset "Model extraction preserves numerical sources" begin
    project=TOML.parsefile(joinpath(root,"Project.toml"))
    @test project["name"]=="GPNPTPB"
    @test project["uuid"]!="fd0e0154-a3e7-4b3b-98df-ed047f11e543"
    @test !haskey(project["deps"],"ElChemTools")
    for (path,hash) in TOML.parsefile(joinpath(root,"docs/source_hashes.toml"))
        @test bytes2hex(sha256(read(joinpath(root,path)) ))==hash
    end
    for dir in ("src","prototypes")
        for (base,_,names) in walkdir(joinpath(root,dir)), name in names
            endswith(name,".jl") || continue
            parsed=Meta.parseall(read(joinpath(base,name),String))
            @test parsed.head != :error
        end
    end
end
