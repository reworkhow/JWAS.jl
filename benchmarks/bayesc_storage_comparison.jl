using Test,JWAS,DataFrames,CSV,Random,Statistics,LinearAlgebra
include(joinpath(@__DIR__,"..","test","unit","bayesc_memory_fixture.jl"))
BLAS.set_num_threads(1)
rows=NamedTuple[]
for traits in (1,2),mode in ("ordinary","sequential"),seed in (73,123)
    a=redirect_stdout(devnull) do
        bayesc_memory_case(;traits,mode,seed,copying=true,n=128,t=96,p=257,
            prediction_all=true,outer=200,burnin=50)
    end
    b=redirect_stdout(devnull) do
        bayesc_memory_case(;traits,mode,seed,n=128,t=96,p=257,
            prediction_all=true,outer=200,burnin=50)
    end
    compare_bayesc_memory(a,b)
    push!(rows,(;traits,mode,seed,copied_bytes=a.storage.backing_matrix_bytes,
        shared_bytes=b.storage.backing_matrix_bytes,copied_parents=a.storage.backing_matrix_count,
        shared_parents=b.storage.backing_matrix_count))
end
CSV.write(joinpath(@__DIR__,"reports","2026-10-31-bayesc-only","storage.csv"),DataFrame(rows))
show(stdout,MIME("text/plain"),DataFrame(rows));println()
