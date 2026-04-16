# Unit tests for multi-trait MCMC with different Bayesian methods
using Test, JWAS, DataFrames, CSV, JWAS.Datasets
using LinearAlgebra
using Random

function exact_mt_bayesc_state_probs(x::AbstractVector,
                                     y_by_trait::AbstractVector,
                                     vare::AbstractMatrix,
                                     var_effect::AbstractMatrix,
                                     prior,
                                     state_labels::AbstractVector)
    xp = dot(x, x)
    Rinv = inv(vare)
    Ginv = inv(var_effect)
    w = [dot(x, y_trait) for y_trait in y_by_trait]
    log_delta = zeros(Float64, length(state_labels))

    for (idx, state) in enumerate(state_labels)
        D = diagm(state)
        lhs = D * Rinv * D * xp + Ginv
        rhs = (Rinv * D)' * w
        ghat = lhs \ rhs
        log_delta[idx] = -0.5 * (log(det(lhs)) - dot(rhs, ghat)) + JWAS.log_prior(prior, 1, state)
    end

    log_delta .-= maximum(log_delta)
    probs = exp.(log_delta)
    probs ./= sum(probs)
    return probs
end

function empirical_mt_bayesc_state_probs(mode::Symbol,
                                         x::AbstractVector,
                                         y_by_trait::AbstractVector,
                                         vare::AbstractMatrix,
                                         var_effect::AbstractMatrix,
                                         prior,
                                         state_labels::AbstractVector;
                                         niter::Int=20_000,
                                         burnin::Int=3_000)
    xArray = [collect(x)]
    xRinvArray = [collect(x)]
    xpRinvx = [dot(x, x)]
    beta = [zeros(1), zeros(1)]
    delta = [zeros(1), zeros(1)]
    alpha = [zeros(1), zeros(1)]
    # With one marker and alpha initialized at zero, ycorr starts at y - x * alpha = y.
    wArray = [collect(y_trait) for y_trait in y_by_trait]
    counts = zeros(Float64, length(state_labels))

    for iter in 1:niter
        if mode == :I
            JWAS._MTBayesABC_samplerI!(xArray, xRinvArray, xpRinvx, wArray, beta, delta, alpha, vare, [var_effect], prior)
        elseif mode == :II
            JWAS._MTBayesABC_samplerII!(xArray, xRinvArray, xpRinvx, wArray, beta, delta, alpha, vare, [var_effect], state_labels, prior)
        else
            error("mode must be :I or :II")
        end

        if iter > burnin
            state = [delta[1][1], delta[2][1]]
            counts[JWAS.annotated_bayesc_mt_state_index(state)] += 1
        end
    end

    counts ./= sum(counts)
    return counts
end

function empirical_mt_bayesc_block_state_probs(x::AbstractVector,
                                               y_by_trait::AbstractVector,
                                               vare::AbstractMatrix,
                                               var_effect::AbstractMatrix,
                                               prior,
                                               state_labels::AbstractVector;
                                               niter::Int=20_000,
                                               burnin::Int=3_000)
    X = reshape(collect(x), :, 1)
    XArray = [X]
    xpRinvx = [dot(x, x)]
    XpRinvX = [X' * X]
    beta = [zeros(1), zeros(1)]
    delta = [zeros(1), zeros(1)]
    alpha = [zeros(1), zeros(1)]
    wArray = [collect(y_trait) for y_trait in y_by_trait]
    counts = zeros(Float64, length(state_labels))
    Rinvw = ones(Float64, length(x))

    for iter in 1:niter
        JWAS._MTBayesABC_block_samplerII!(XArray, xpRinvx, X, XpRinvX,
                                          wArray, beta, delta, alpha,
                                          vare, [var_effect], state_labels,
                                          prior, Rinvw)
        if iter > burnin
            state = [delta[1][1], delta[2][1]]
            counts[JWAS.annotated_bayesc_mt_state_index(state)] += 1
        end
    end

    counts ./= sum(counts)
    return counts
end

phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
genofile = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
pedfile = Datasets.dataset("pedigree.txt", dataset_name="demo_7animals")
phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstring=["NA"])

R = [1.0 0.5; 0.5 1.0]

@testset "Multi-trait BayesC" begin
    G = [1.0 0.5; 0.5 1.0]
    global geno = get_genotypes(genofile, G, separator=',', method="BayesC")
    model = build_model("y1 = intercept + geno\ny2 = intercept + geno", R)

    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    output_folder="test_mt_bayesc",
                    seed=123)

    @test haskey(output, "location parameters")
    @test haskey(output, "residual variance")
    @test haskey(output, "marker effects geno")
    @test size(output["residual variance"], 1) == 4  # 2x2 covariance flattened
    isdir("test_mt_bayesc") && rm("test_mt_bayesc", recursive=true, force=true)
end

@testset "Multi-trait BayesC sampler default and auto selection" begin
    G = [1.0 0.5; 0.5 1.0]

    geno_default = get_genotypes(genofile, G, separator=',', method="BayesC")
    geno_auto = get_genotypes(genofile, G, separator=',', method="BayesC", multi_trait_sampler=:auto)

    @test geno_default.multi_trait_sampler == :I
    @test geno_auto.multi_trait_sampler == :auto
end

@testset "Multi-trait annotated BayesC" begin
    G = [1.0 0.5; 0.5 1.0]
    annotations = [
        0.0 1.0
        1.0 0.0
        1.0 1.0
        0.0 0.0
        0.5 0.5
    ]
    start_pi = Dict(
        [0.0, 0.0] => 0.45,
        [1.0, 0.0] => 0.20,
        [0.0, 1.0] => 0.15,
        [1.0, 1.0] => 0.20,
    )
    phenotypes_mt = DataFrame(
        ID=copy(phenotypes.ID),
        y1=copy(phenotypes.y1),
        y2=Float32.(coalesce.(phenotypes.y1, 0.0)),
    )
    global annotated_mt_geno = get_genotypes(
        genofile,
        G;
        separator=',',
        method="BayesC",
        quality_control=false,
        annotations=annotations,
        Pi=start_pi,
    )
    model = build_model("y1 = intercept + annotated_mt_geno\ny2 = intercept + annotated_mt_geno", R)

    output = runMCMC(
        model,
        phenotypes_mt,
        chain_length=40,
        burnin=10,
        output_samples_frequency=10,
        output_folder="test_mt_annotated_bayesc",
        seed=123,
    )

    @test model.M[1].annotations !== false
    @test model.M[1].annotations.nsteps == 3
    @test size(model.M[1].annotations.snp_pi, 2) == 4
    @test haskey(output, "annotation coefficients annotated_mt_geno")
    @test haskey(output, "pi_annotated_mt_geno")
    @test nrow(output["pi_annotated_mt_geno"]) == 4
    isdir("test_mt_annotated_bayesc") && rm("test_mt_annotated_bayesc", recursive=true, force=true)
end

@testset "Multi-trait annotated BayesC sampler II override" begin
    G = [1.0 0.5; 0.5 1.0]
    annotations = [
        0.0 1.0
        1.0 0.0
        1.0 1.0
        0.0 0.0
        0.5 0.5
    ]
    start_pi = Dict(
        [0.0, 0.0] => 0.45,
        [1.0, 0.0] => 0.20,
        [0.0, 1.0] => 0.15,
        [1.0, 1.0] => 0.20,
    )
    phenotypes_mt = DataFrame(
        ID=copy(phenotypes.ID),
        y1=copy(phenotypes.y1),
        y2=Float32.(coalesce.(phenotypes.y1, 0.0)),
    )
    global annotated_mt_geno_sampler2 = get_genotypes(
        genofile,
        G;
        separator=',',
        method="BayesC",
        quality_control=false,
        annotations=annotations,
        Pi=start_pi,
        multi_trait_sampler=:II,
    )
    model = build_model("y1 = intercept + annotated_mt_geno_sampler2\ny2 = intercept + annotated_mt_geno_sampler2", R)

    @test JWAS.mt_bayesc_sampler_mode(model.M[1], 2) == :II

    output = runMCMC(
        model,
        phenotypes_mt,
        chain_length=40,
        burnin=10,
        output_samples_frequency=10,
        output_folder="test_mt_annotated_bayesc_sampler2",
        seed=123,
    )

    @test haskey(output, "annotation coefficients annotated_mt_geno_sampler2")
    @test haskey(output, "pi_annotated_mt_geno_sampler2")
    isdir("test_mt_annotated_bayesc_sampler2") && rm("test_mt_annotated_bayesc_sampler2", recursive=true, force=true)
end

@testset "Multi-trait BayesC dense sampler II run" begin
    G = [1.0 0.5; 0.5 1.0]

    global plain_mt_dense_sampler2 = get_genotypes(
        genofile,
        G;
        separator=',',
        method="BayesC",
        quality_control=false,
        multi_trait_sampler=:II,
    )
    model = build_model("y1 = intercept + plain_mt_dense_sampler2\ny2 = intercept + plain_mt_dense_sampler2", R)

    output = runMCMC(
        model,
        phenotypes,
        chain_length=40,
        burnin=10,
        output_samples_frequency=10,
        output_folder="test_mt_bayesc_dense_sampler2",
        seed=123,
    )

    @test haskey(output, "marker effects plain_mt_dense_sampler2")
    @test haskey(output, "location parameters")
    isdir("test_mt_bayesc_dense_sampler2") && rm("test_mt_bayesc_dense_sampler2", recursive=true, force=true)
end

@testset "Multi-trait BayesC fast_blocks run" begin
    G = [1.0 0.5; 0.5 1.0]

    global plain_mt_fastblocks = get_genotypes(
        genofile,
        G;
        separator=',',
        method="BayesC",
        quality_control=false,
    )
    model = build_model("y1 = intercept + plain_mt_fastblocks\ny2 = intercept + plain_mt_fastblocks", R)

    output = runMCMC(
        model,
        phenotypes,
        chain_length=40,
        burnin=10,
        output_samples_frequency=10,
        output_folder="test_mt_bayesc_fastblocks",
        seed=123,
        fast_blocks=true,
    )

    @test haskey(output, "marker effects plain_mt_fastblocks")
    @test haskey(output, "location parameters")
    isdir("test_mt_bayesc_fastblocks") && rm("test_mt_bayesc_fastblocks", recursive=true, force=true)
end

@testset "Multi-trait BayesC fast_blocks sampler II run" begin
    G = [1.0 0.5; 0.5 1.0]

    global plain_mt_fastblocks_sampler2 = get_genotypes(
        genofile,
        G;
        separator=',',
        method="BayesC",
        quality_control=false,
        multi_trait_sampler=:II,
    )
    model = build_model("y1 = intercept + plain_mt_fastblocks_sampler2\ny2 = intercept + plain_mt_fastblocks_sampler2", R)

    output = runMCMC(
        model,
        phenotypes,
        chain_length=40,
        burnin=10,
        output_samples_frequency=10,
        output_folder="test_mt_bayesc_fastblocks_sampler2",
        seed=123,
        fast_blocks=true,
    )

    @test haskey(output, "marker effects plain_mt_fastblocks_sampler2")
    @test haskey(output, "location parameters")
    isdir("test_mt_bayesc_fastblocks_sampler2") && rm("test_mt_bayesc_fastblocks_sampler2", recursive=true, force=true)
end

@testset "Multi-trait annotated BayesC fast_blocks sampler II run" begin
    G = [1.0 0.5; 0.5 1.0]
    annotations = [
        0.0 1.0
        1.0 0.0
        1.0 1.0
        0.0 0.0
        0.5 0.5
    ]
    start_pi = Dict(
        [0.0, 0.0] => 0.45,
        [1.0, 0.0] => 0.20,
        [0.0, 1.0] => 0.15,
        [1.0, 1.0] => 0.20,
    )
    phenotypes_mt = DataFrame(
        ID=copy(phenotypes.ID),
        y1=copy(phenotypes.y1),
        y2=Float32.(coalesce.(phenotypes.y1, 0.0)),
    )

    global annotated_mt_fastblocks_sampler2 = get_genotypes(
        genofile,
        G;
        separator=',',
        method="BayesC",
        quality_control=false,
        annotations=annotations,
        Pi=start_pi,
        multi_trait_sampler=:II,
    )
    model = build_model("y1 = intercept + annotated_mt_fastblocks_sampler2\ny2 = intercept + annotated_mt_fastblocks_sampler2", R)

    output = runMCMC(
        model,
        phenotypes_mt,
        chain_length=40,
        burnin=10,
        output_samples_frequency=10,
        output_folder="test_mt_annotated_bayesc_fastblocks_sampler2",
        seed=123,
        fast_blocks=true,
    )

    @test haskey(output, "annotation coefficients annotated_mt_fastblocks_sampler2")
    @test haskey(output, "pi_annotated_mt_fastblocks_sampler2")
    isdir("test_mt_annotated_bayesc_fastblocks_sampler2") && rm("test_mt_annotated_bayesc_fastblocks_sampler2", recursive=true, force=true)
end

@testset "Multi-trait BayesC independent fast_blocks sampler I run" begin
    G = [1.0 0.5; 0.5 1.0]

    global plain_mt_independent_sampler1 = get_genotypes(
        genofile,
        G;
        separator=',',
        method="BayesC",
        quality_control=false,
        multi_trait_sampler=:I,
    )
    model = build_model("y1 = intercept + plain_mt_independent_sampler1\ny2 = intercept + plain_mt_independent_sampler1", R)
    @test JWAS.mt_bayesc_sampler_mode(model.M[1], 2) == :I

    output = runMCMC(
        model,
        phenotypes,
        chain_length=12,
        burnin=2,
        output_samples_frequency=10,
        output_folder="test_mt_bayesc_independent_sampler1",
        printout_frequency=13,
        seed=123,
        fast_blocks=[1, 3, 5],
        independent_blocks=true,
        outputEBV=false,
        output_heritability=false,
    )

    @test haskey(output, "marker effects plain_mt_independent_sampler1")
    @test haskey(output, "location parameters")
    @test model.MCMCinfo.independent_blocks == true
    @test model.MCMCinfo.fast_blocks == [1, 3, 5]
    isdir("test_mt_bayesc_independent_sampler1") && rm("test_mt_bayesc_independent_sampler1", recursive=true, force=true)
end

@testset "Multi-trait BayesC independent fast_blocks sampler II run" begin
    G = [1.0 0.5; 0.5 1.0]

    global plain_mt_independent_sampler2 = get_genotypes(
        genofile,
        G;
        separator=',',
        method="BayesC",
        quality_control=false,
        multi_trait_sampler=:II,
    )
    model = build_model("y1 = intercept + plain_mt_independent_sampler2\ny2 = intercept + plain_mt_independent_sampler2", R)
    @test JWAS.mt_bayesc_sampler_mode(model.M[1], 2) == :II

    output = runMCMC(
        model,
        phenotypes,
        chain_length=12,
        burnin=2,
        output_samples_frequency=10,
        output_folder="test_mt_bayesc_independent_sampler2",
        printout_frequency=13,
        seed=123,
        fast_blocks=[1, 3, 5],
        independent_blocks=true,
        outputEBV=false,
        output_heritability=false,
    )

    @test haskey(output, "marker effects plain_mt_independent_sampler2")
    @test haskey(output, "location parameters")
    @test model.MCMCinfo.independent_blocks == true
    @test model.MCMCinfo.fast_blocks == [1, 3, 5]
    isdir("test_mt_bayesc_independent_sampler2") && rm("test_mt_bayesc_independent_sampler2", recursive=true, force=true)
end

@testset "Multi-trait annotated BayesC independent fast_blocks sampler I run" begin
    G = [1.0 0.5; 0.5 1.0]
    annotations = [
        0.0 1.0
        1.0 0.0
        1.0 1.0
        0.0 0.0
        0.5 0.5
    ]
    start_pi = Dict(
        [0.0, 0.0] => 0.45,
        [1.0, 0.0] => 0.20,
        [0.0, 1.0] => 0.15,
        [1.0, 1.0] => 0.20,
    )
    phenotypes_mt = DataFrame(
        ID=copy(phenotypes.ID),
        y1=copy(phenotypes.y1),
        y2=Float32.(coalesce.(phenotypes.y1, 0.0)),
    )

    global annotated_mt_independent_sampler1 = get_genotypes(
        genofile,
        G;
        separator=',',
        method="BayesC",
        quality_control=false,
        annotations=annotations,
        Pi=start_pi,
        multi_trait_sampler=:I,
    )
    model = build_model("y1 = intercept + annotated_mt_independent_sampler1\ny2 = intercept + annotated_mt_independent_sampler1", R)
    @test JWAS.mt_bayesc_sampler_mode(model.M[1], 2) == :I

    output = runMCMC(
        model,
        phenotypes_mt,
        chain_length=12,
        burnin=2,
        output_samples_frequency=10,
        output_folder="test_mt_annotated_bayesc_independent_sampler1",
        printout_frequency=13,
        seed=123,
        fast_blocks=[1, 3, 5],
        independent_blocks=true,
        outputEBV=false,
        output_heritability=false,
    )

    @test haskey(output, "annotation coefficients annotated_mt_independent_sampler1")
    @test haskey(output, "pi_annotated_mt_independent_sampler1")
    @test model.MCMCinfo.independent_blocks == true
    @test model.MCMCinfo.fast_blocks == [1, 3, 5]
    isdir("test_mt_annotated_bayesc_independent_sampler1") && rm("test_mt_annotated_bayesc_independent_sampler1", recursive=true, force=true)
end

@testset "Multi-trait annotated BayesC independent fast_blocks sampler II run" begin
    G = [1.0 0.5; 0.5 1.0]
    annotations = [
        0.0 1.0
        1.0 0.0
        1.0 1.0
        0.0 0.0
        0.5 0.5
    ]
    start_pi = Dict(
        [0.0, 0.0] => 0.45,
        [1.0, 0.0] => 0.20,
        [0.0, 1.0] => 0.15,
        [1.0, 1.0] => 0.20,
    )
    phenotypes_mt = DataFrame(
        ID=copy(phenotypes.ID),
        y1=copy(phenotypes.y1),
        y2=Float32.(coalesce.(phenotypes.y1, 0.0)),
    )

    global annotated_mt_independent_sampler2 = get_genotypes(
        genofile,
        G;
        separator=',',
        method="BayesC",
        quality_control=false,
        annotations=annotations,
        Pi=start_pi,
        multi_trait_sampler=:II,
    )
    model = build_model("y1 = intercept + annotated_mt_independent_sampler2\ny2 = intercept + annotated_mt_independent_sampler2", R)
    @test JWAS.mt_bayesc_sampler_mode(model.M[1], 2) == :II

    output = runMCMC(
        model,
        phenotypes_mt,
        chain_length=12,
        burnin=2,
        output_samples_frequency=10,
        output_folder="test_mt_annotated_bayesc_independent_sampler2",
        printout_frequency=13,
        seed=123,
        fast_blocks=[1, 3, 5],
        independent_blocks=true,
        outputEBV=false,
        output_heritability=false,
    )

    @test haskey(output, "annotation coefficients annotated_mt_independent_sampler2")
    @test haskey(output, "pi_annotated_mt_independent_sampler2")
    @test model.MCMCinfo.independent_blocks == true
    @test model.MCMCinfo.fast_blocks == [1, 3, 5]
    isdir("test_mt_annotated_bayesc_independent_sampler2") && rm("test_mt_annotated_bayesc_independent_sampler2", recursive=true, force=true)
end

@testset "Multi-trait annotated BayesC samplers share the same target posterior" begin
    x = [1.0, -0.5, 0.75]
    y_by_trait = [
        [0.8, -0.1, 0.3],
        [0.2, 0.6, -0.4],
    ]
    vare = [1.0 0.25; 0.25 0.9]
    var_effect = [0.7 0.15; 0.15 0.8]
    state_labels = JWAS.annotated_bayesc_mt_state_keys()
    prior = JWAS.MarkerSpecificPiPrior(reshape([0.35, 0.20, 0.15, 0.30], 1, 4))

    exact = exact_mt_bayesc_state_probs(x, y_by_trait, vare, var_effect, prior, state_labels)

    Random.seed!(20260410)
    sampler_i = empirical_mt_bayesc_state_probs(:I, x, y_by_trait, vare, var_effect, prior, state_labels)
    Random.seed!(20260410)
    sampler_ii = empirical_mt_bayesc_state_probs(:II, x, y_by_trait, vare, var_effect, prior, state_labels)

    @test maximum(abs.(sampler_i .- exact)) < 0.02
    @test maximum(abs.(sampler_ii .- exact)) < 0.02
    @test maximum(abs.(sampler_i .- sampler_ii)) < 0.02
end

@testset "Multi-trait BayesC samplers share the same target posterior" begin
    x = [1.0, -0.5, 0.75]
    y_by_trait = [
        [0.8, -0.1, 0.3],
        [0.2, 0.6, -0.4],
    ]
    vare = [1.0 0.25; 0.25 0.9]
    var_effect = [0.7 0.15; 0.15 0.8]
    state_labels = JWAS.annotated_bayesc_mt_state_keys()
    prior = JWAS.GlobalPiPrior(Dict(
        [0.0, 0.0] => 0.35,
        [1.0, 0.0] => 0.20,
        [0.0, 1.0] => 0.15,
        [1.0, 1.0] => 0.30,
    ))

    exact = exact_mt_bayesc_state_probs(x, y_by_trait, vare, var_effect, prior, state_labels)

    Random.seed!(20260411)
    sampler_i = empirical_mt_bayesc_state_probs(:I, x, y_by_trait, vare, var_effect, prior, state_labels)
    Random.seed!(20260411)
    sampler_ii = empirical_mt_bayesc_state_probs(:II, x, y_by_trait, vare, var_effect, prior, state_labels)

    @test maximum(abs.(sampler_i .- exact)) < 0.02
    @test maximum(abs.(sampler_ii .- exact)) < 0.02
    @test maximum(abs.(sampler_i .- sampler_ii)) < 0.02
end

@testset "Multi-trait BayesC block sampler II matches dense sampler II on a one-marker case" begin
    x = [1.0, -0.5, 0.75]
    y_by_trait = [
        [0.8, -0.1, 0.3],
        [0.2, 0.6, -0.4],
    ]
    vare = [1.0 0.25; 0.25 0.9]
    var_effect = [0.7 0.15; 0.15 0.8]
    state_labels = JWAS.annotated_bayesc_mt_state_keys()
    plain_prior = JWAS.GlobalPiPrior(Dict(
        [0.0, 0.0] => 0.35,
        [1.0, 0.0] => 0.20,
        [0.0, 1.0] => 0.15,
        [1.0, 1.0] => 0.30,
    ))
    annotated_prior = JWAS.MarkerSpecificPiPrior(reshape([0.35, 0.20, 0.15, 0.30], 1, 4))

    exact_plain = exact_mt_bayesc_state_probs(x, y_by_trait, vare, var_effect, plain_prior, state_labels)
    exact_annotated = exact_mt_bayesc_state_probs(x, y_by_trait, vare, var_effect, annotated_prior, state_labels)

    Random.seed!(20260412)
    dense_plain = empirical_mt_bayesc_state_probs(:II, x, y_by_trait, vare, var_effect, plain_prior, state_labels)
    Random.seed!(20260412)
    block_plain = empirical_mt_bayesc_block_state_probs(x, y_by_trait, vare, var_effect, plain_prior, state_labels)

    Random.seed!(20260413)
    dense_annotated = empirical_mt_bayesc_state_probs(:II, x, y_by_trait, vare, var_effect, annotated_prior, state_labels)
    Random.seed!(20260413)
    block_annotated = empirical_mt_bayesc_block_state_probs(x, y_by_trait, vare, var_effect, annotated_prior, state_labels)

    @test maximum(abs.(block_plain .- exact_plain)) < 0.02
    @test maximum(abs.(block_plain .- dense_plain)) < 0.02
    @test maximum(abs.(block_annotated .- exact_annotated)) < 0.02
    @test maximum(abs.(block_annotated .- dense_annotated)) < 0.02
end

@testset "Multi-trait BayesC block wrapper dispatch" begin
    G = [1.0 0.5; 0.5 1.0]
    annotations = [
        0.0 1.0
        1.0 0.0
        1.0 1.0
        0.0 0.0
        0.5 0.5
    ]
    start_pi = Dict(
        [0.0, 0.0] => 0.45,
        [1.0, 0.0] => 0.20,
        [0.0, 1.0] => 0.15,
        [1.0, 1.0] => 0.20,
    )

    plain_geno = get_genotypes(
        genofile,
        G;
        separator=',',
        method="BayesC",
        quality_control=false,
        storage=:dense,
    )
    plain_auto_geno = get_genotypes(
        genofile,
        G;
        separator=',',
        method="BayesC",
        quality_control=false,
        storage=:dense,
        multi_trait_sampler=:auto,
    )
    annotated_geno = get_genotypes(
        genofile,
        G;
        separator=',',
        method="BayesC",
        quality_control=false,
        storage=:dense,
        annotations=annotations,
        Pi=start_pi,
        multi_trait_sampler=:II,
    )

    plain_prior, plain_mode = JWAS.mt_bayesc_block_context(plain_geno, 2)
    plain_auto_prior, plain_auto_mode = JWAS.mt_bayesc_block_context(plain_auto_geno, 2)
    annotated_prior, annotated_mode = JWAS.mt_bayesc_block_context(annotated_geno, 2)

    @test plain_prior isa JWAS.GlobalPiPrior
    @test plain_mode == :I
    @test plain_auto_prior isa JWAS.GlobalPiPrior
    @test plain_auto_mode == :II
    @test annotated_prior isa JWAS.MarkerSpecificPiPrior
    @test annotated_mode == :II
end

@testset "Multi-trait RR-BLUP" begin
    G = [1.0 0.5; 0.5 1.0]
    global geno = get_genotypes(genofile, G, separator=',', method="RR-BLUP")
    model = build_model("y1 = intercept + geno\ny2 = intercept + geno", R)

    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    output_folder="test_mt_rrblup",
                    seed=123)

    @test haskey(output, "location parameters")
    @test haskey(output, "residual variance")
    isdir("test_mt_rrblup") && rm("test_mt_rrblup", recursive=true, force=true)
end

@testset "Multi-trait with pedigree" begin
    G = [1.0 0.5; 0.5 1.0]
    ped = get_pedigree(pedfile, separator=",", header=true)
    model = build_model("y1 = intercept + ID\ny2 = intercept + ID", R)
    set_random(model, "ID", ped, G)

    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    output_folder="test_mt_ped",
                    seed=123)

    @test haskey(output, "polygenic effects covariance matrix")
    @test haskey(output, "location parameters")
    @test haskey(output, "residual variance")
    isdir("test_mt_ped") && rm("test_mt_ped", recursive=true, force=true)
end

@testset "Multi-trait with i.i.d. random effect" begin
    G_geno = [1.0 0.5; 0.5 1.0]
    global geno = get_genotypes(genofile, G_geno, separator=',', method="RR-BLUP")
    model = build_model("y1 = intercept + x2 + geno\ny2 = intercept + x2 + geno", R)
    G = [0.5 0.1; 0.1 0.5]
    set_random(model, "x2", G)

    output = runMCMC(model, phenotypes,
                    chain_length=100,
                    burnin=20,
                    output_samples_frequency=10,
                    output_folder="test_mt_iid",
                    seed=123)

    @test haskey(output, "location parameters")
    @test haskey(output, "residual variance")
    isdir("test_mt_iid") && rm("test_mt_iid", recursive=true, force=true)
end
