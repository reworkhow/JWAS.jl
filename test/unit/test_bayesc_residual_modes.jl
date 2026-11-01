using Test,JWAS,Random,LinearAlgebra
@testset "BayesC residual reconstruction ST and MT" begin
    for T in (Float32,Float64), independent in (false,true),seed in (73,123,456)
        rng=MersenneTwister(seed)
        parent_X=randn(rng,T,41,11);original=copy(parent_X)
        X=view(parent_X,1:29,:)
        blocks=[view(X,:,s:min(s+3,11)) for s in 1:4:11]
        grams=[b'b for b in blocks]
        diagonals=T[dot(view(X,:,j),view(X,:,j)) for j in 1:11]
        for sampler in (:single,:I,:II)
            Random.seed!(seed)
            alpha=[randn(rng,T,11) for _ in 1:(sampler==:single ? 1 : 2)]
            beta=deepcopy(alpha);delta=[ones(T,11) for _ in alpha]
            residual=[randn(rng,T,29) for _ in alpha]
            before=deepcopy(residual);oldalpha=deepcopy(alpha)
            if sampler==:single
                JWAS.BayesABC_block!(blocks,diagonals,X,grams,residual[1],
                    alpha[1],beta[1],delta[1],one(T),ones(T,11),0.65,ones(T,29),
                    independent)
            else
                pi=Dict([0.0,0.0]=>0.45,[1.0,0.0]=>0.20,
                        [0.0,1.0]=>0.15,[1.0,1.0]=>0.20)
                prior=JWAS.GlobalPiPrior(pi)
                R=T[1 0.1;0.1 1]
                G=[T[0.8 0.15;0.15 0.9] for _ in 1:11]
                if sampler==:I
                    fn=independent ? JWAS._MTBayesABC_block_samplerI_independent! :
                                     JWAS._MTBayesABC_block_samplerI!
                    fn(blocks,diagonals,X,grams,residual,beta,delta,alpha,R,G,prior,
                       ones(T,29))
                else
                    fn=independent ? JWAS._MTBayesABC_block_samplerII_independent! :
                                     JWAS._MTBayesABC_block_samplerII!
                    fn(blocks,diagonals,X,grams,residual,beta,delta,alpha,R,G,collect(keys(pi)),
                       prior,ones(T,29))
                end
            end
            for k in eachindex(alpha)
                @test isapprox(residual[k],before[k]+X*(oldalpha[k]-alpha[k]);atol=200eps(T),rtol=200eps(T))
                @test all(isfinite,residual[k]) && all(isfinite,alpha[k])
            end
            @test parent_X==original
        end
    end
end

@testset "Genotype-typed residual multipliers" begin
    for XT in (Float32,Float64), AT in (Float32,Float64), YT in (Float32,Float64)
        rng=MersenneTwister(73)
        X=randn(rng,XT,19,7);delta=randn(rng,AT,7);y=randn(rng,YT,19)
        expected=copy(y)
        mul!(expected,X,delta,one(XT),one(XT))
        JWAS.bayesc_block_residual!(y,X,delta)
        @test y==expected
        @test eltype(y) === YT
    end
end

@testset "Shared ABC wrapper preserves BayesB arithmetic" begin
    rng=MersenneTwister(73)
    X=randn(rng,Float32,29,11)
    blocks=[view(X,:,s:min(s+3,11)) for s in 1:4:11]
    grams=[b'b for b in blocks]
    diagonals=Float32[dot(view(X,:,j),view(X,:,j)) for j in 1:11]
    for independent in (false,true)
        alpha=[randn(rng,Float32,11)];beta=deepcopy(alpha);delta=[ones(Float32,11)]
        g=(;method="BayesB",MArray=blocks,mpRinvm=diagonals,genotypes=X,MpRinvM=grams,
            α=alpha,β=beta,δ=delta,π=0.65)
        expected_alpha=deepcopy(alpha);expected_beta=deepcopy(beta);expected_delta=deepcopy(delta)
        residual=randn(rng,Float32,29);expected=copy(residual)
        Random.seed!(123)
        JWAS.BayesABC_block!(g,residual,1f0,ones(Float32,11),ones(Float32,29),independent)
        Random.seed!(123)
        JWAS.BayesABC_block!(blocks,diagonals,X,grams,expected,
            expected_alpha[1],expected_beta[1],expected_delta[1],1f0,ones(Float32,11),0.65,
            ones(Float32,29),independent;residual_update! = JWAS.legacy_abc_block_residual!)
        @test residual==expected
        @test alpha==expected_alpha
    end
end
