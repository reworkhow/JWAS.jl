# Use this same script with --project pointing to untouched upstream or this branch.
# Arguments: output_directory upstream|candidate.
using Test,JWAS,DataFrames,CSV,Random,Statistics,LinearAlgebra
include(joinpath(@__DIR__,"..","test","unit","bayesc_memory_fixture.jl"))
BLAS.set_num_threads(1)
outdir=abspath(ARGS[1]); upstream=ARGS[2]=="upstream";mkpath(outdir)
println("Source: ",pathof(JWAS),"; Julia: ",VERSION);flush(stdout)
rows=NamedTuple[]
for traits in (1,2), sampler in (traits==1 ? (:I,) : (:I,:II)),
    precision in (Float32,Float64), mode in ("ordinary","sequential","independent"),
    seed in (73,123,456)
    name="traits$(traits)_sampler$(sampler)_$(precision)_$(mode)_seed$(seed)"
    result=redirect_stdout(devnull) do
        bayesc_memory_case(;traits,sampler,precision,mode,seed,upstream,
            n=40,t=28,p=19,outer=500,burnin=150,freq=10,folder=joinpath(outdir,name))
    end
    CSV.write(joinpath(outdir,name*"_animals.csv"),result.ebv)
    CSV.write(joinpath(outdir,name*"_markers.csv"),result.marker)
    push!(rows,(;name,traits,sampler=string(sampler),precision=string(precision),mode,seed,annotated=false,
        animals=traits*12,markers=traits*19))
    println("Completed ",name);flush(stdout)
end
for traits in (1,2), sampler in (traits==1 ? (:I,) : (:I,:II)),
    precision in (Float32,Float64), seed in (73,123,456)
    mode="sequential"
    name="traits$(traits)_sampler$(sampler)_$(precision)_annotated_seed$(seed)"
    result=redirect_stdout(devnull) do
        bayesc_memory_case(;traits,sampler,precision,mode,seed,upstream,annotated=true,
            n=40,t=28,p=19,outer=500,burnin=150,freq=10,folder=joinpath(outdir,name))
    end
    CSV.write(joinpath(outdir,name*"_animals.csv"),result.ebv)
    CSV.write(joinpath(outdir,name*"_markers.csv"),result.marker)
    push!(rows,(;name,traits,sampler=string(sampler),precision=string(precision),mode,seed,annotated=true,
        animals=traits*12,markers=traits*19))
    println("Completed ",name);flush(stdout)
end
CSV.write(joinpath(outdir,"cases.csv"),DataFrame(rows))
println("All ",length(rows)," cases completed.")
