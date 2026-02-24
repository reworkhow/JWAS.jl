# Unit tests for pedigree algorithms (AInverse, inbreeding, get_info)
using Test, JWAS, DataFrames, CSV, JWAS.Datasets
using JWAS.PedModule
using SparseArrays, LinearAlgebra

pedfile = Datasets.dataset("pedigree.txt", dataset_name="demo_7animals")

@testset "Pedigree Algorithms" begin
    ped = get_pedigree(pedfile, separator=",", header=true)

    @testset "AInverse" begin
        Ai = PedModule.AInverse(ped)
        n = length(ped.idMap)
        @test size(Ai) == (n, n)
        @test issparse(Ai)
        # A-inverse should be symmetric
        @test Ai ≈ Ai'
        # Diagonal elements should be positive
        @test all(diag(Ai) .> 0)
    end

    @testset "AInverse vs AInverseSlow" begin
        Ai_fast = PedModule.AInverse(ped)
        Ai_slow = PedModule.AInverseSlow(ped)
        # Both methods should produce the same result
        @test Matrix(Ai_fast) ≈ Matrix(Ai_slow) atol=1e-6
    end

    @testset "Inbreeding coefficients" begin
        inbreeding = PedModule.getInbreeding(ped)
        @test length(inbreeding) == length(ped.idMap)
        # All inbreeding coefficients should be between 0 and 1
        @test all(0 .<= inbreeding .<= 1)
        # Founders should have inbreeding = 0
        for (id, node) in ped.idMap
            if node.sire == "missing" && node.dam == "missing"
                @test node.f == 0.0
            end
        end
    end

    @testset "getIDs" begin
        ids = PedModule.getIDs(ped)
        @test length(ids) == length(ped.idMap)
        @test length(unique(ids)) == length(ids)  # all unique
    end

    @testset "get_info with Ai=true" begin
        IDs, Ai, inbreeding = get_info(ped, Ai=true)
        @test length(IDs) == length(ped.idMap)
        @test size(Ai) == (length(IDs), length(IDs))
        @test length(inbreeding) == length(IDs)
    end

    @testset "Pedigree from DataFrame" begin
        df = DataFrame(ID=["a","b","c"], Sire=["missing","missing","a"], Dam=["missing","missing","b"])
        ped2 = get_pedigree(df)
        @test length(ped2.idMap) == 3
        @test ped2.idMap["c"].sire == "a"
        @test ped2.idMap["c"].dam == "b"
    end
end
