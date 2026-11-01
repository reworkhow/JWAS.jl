include("bayesc_memory_fixture.jl")
@testset "BayesC shared memory, ST/MT and ordinary/fast blocks" begin
    for traits in (1,2), sampler in (traits==1 ? (:I,) : (:I,:II)),
        precision in (Float32,Float64), mode in ("ordinary","sequential","independent")
        copied=bayesc_memory_case(;traits,sampler,precision,mode,copying=true)
        shared=bayesc_memory_case(;traits,sampler,precision,mode)
        compare_bayesc_memory(copied,shared)
    end
    for sampler in (:I,:II), precision in (Float32,Float64)
        copied=bayesc_memory_case(;traits=2,sampler,precision,mode="sequential",annotated=true,copying=true)
        shared=bayesc_memory_case(;traits=2,sampler,precision,mode="sequential",annotated=true)
        compare_bayesc_memory(copied,shared)
    end
end
@testset "BayesC-only defaults and storage scope" begin
    for method in ("BayesB","BayesR","RR-BLUP")
        X=Float32[0 1 2;2 1 0;1 2 1]
        g=get_genotypes(X;method,quality_control=false,center=false)
        @test g.genotypes !== X
        @test g.genotypes==X
    end
    X=reshape(Float32.(1:84),12,7)
    shared=JWAS.GibbsMats(view(X,1:8,:),ones(Float32,8);fast_blocks=[1,4,7],use_views=true)
    copied=JWAS.GibbsMats(X,ones(Float32,12);fast_blocks=[1,4,7])
    @test all(b->Base.mightalias(b,X),shared.XArray)
    @test all(b->!Base.mightalias(b,X),copied.XArray)
    @test_throws ErrorException JWAS.contiguous_genotype_rows(["a","c"],["a","b","c"])
end
