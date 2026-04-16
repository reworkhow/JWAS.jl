using CSV
using DataFrames
using JWAS
using JWAS.Datasets
using LinearAlgebra

BLAS.set_num_threads(1)

const REPORT_DATE = "2026-04-15"
const DEFAULT_OUTFILE = joinpath(
    @__DIR__,
    "reports",
    "$(REPORT_DATE)-independent-blocks-fast-block-parallelism-smoke-results-threads$(Threads.nthreads()).csv",
)

function common_annotations()
    return [
        0.0 1.0
        1.0 0.0
        1.0 1.0
        0.0 0.0
        0.5 0.5
    ]
end

function multitrait_pi()
    return Dict(
        [0.0, 0.0] => 0.45,
        [1.0, 0.0] => 0.20,
        [0.0, 1.0] => 0.15,
        [1.0, 1.0] => 0.20,
    )
end

function bind_genotype!(name::AbstractString, geno)
    @eval $(Symbol(name)) = $geno
    return nothing
end

function build_benchmark_case(case, genofile, phenotypes)
    annotations = common_annotations()
    if case.multitrait
        G = [1.0 0.5; 0.5 1.0]
        kwargs = (
            separator=',',
            method="BayesC",
            quality_control=false,
            multi_trait_sampler=case.sampler,
        )
        geno = if case.annotated
            get_genotypes(genofile, G; kwargs..., annotations=annotations, Pi=multitrait_pi())
        else
            get_genotypes(genofile, G; kwargs...)
        end
        bind_genotype!(case.name, geno)
        phenotypes_mt = DataFrame(
            ID=copy(phenotypes.ID),
            y1=copy(phenotypes.y1),
            y2=Float32.(coalesce.(phenotypes.y1, 0.0)),
        )
        model = build_model("y1 = intercept + $(case.name)\ny2 = intercept + $(case.name)", [1.0 0.5; 0.5 1.0])
        return model, phenotypes_mt
    end

    kwargs = (
        separator=',',
        method=case.method,
        quality_control=false,
    )
    geno = if case.annotated
        get_genotypes(genofile, 1.0; kwargs..., annotations=annotations)
    else
        get_genotypes(genofile, 1.0; kwargs...)
    end
    bind_genotype!(case.name, geno)
    model = build_model("y1 = intercept + $(case.name)", 1.0)
    return model, phenotypes
end

function run_case!(rows, case, run_mode, genofile, phenotypes, tmpdir)
    model, pheno = build_benchmark_case(case, genofile, phenotypes)
    outdir = joinpath(tmpdir, "$(case.name)_$(run_mode.label)")
    elapsed = @elapsed output = runMCMC(
        model,
        pheno;
        chain_length=20,
        burnin=5,
        output_samples_frequency=10,
        output_folder=outdir,
        printout_model_info=false,
        printout_frequency=21,
        seed=20260415,
        fast_blocks=[1, 3, 5],
        independent_blocks=run_mode.independent_blocks,
        outputEBV=false,
        output_heritability=false,
    )

    push!(rows, (
        threads=Threads.nthreads(),
        blas_threads=BLAS.get_num_threads(),
        label=case.label,
        method=case.method,
        annotated=case.annotated,
        multitrait=case.multitrait,
        sampler=String(case.sampler),
        run_mode=run_mode.label,
        independent_blocks=run_mode.independent_blocks,
        chain_length=model.MCMCinfo.chain_length,
        block_starts=join(model.MCMCinfo.fast_blocks, ";"),
        nmarkers=model.M[1].nMarkers,
        elapsed_seconds=elapsed,
        output_key_count=length(keys(output)),
    ))
    return nothing
end

function main(args=ARGS)
    outfile = isempty(args) ? DEFAULT_OUTFILE : abspath(args[1])
    mkpath(dirname(outfile))

    genofile = Datasets.dataset("genotypes.txt", dataset_name="demo_7animals")
    phenofile = Datasets.dataset("phenotypes.txt", dataset_name="demo_7animals")
    phenotypes = CSV.read(phenofile, DataFrame, delim=',', missingstring=["NA"])

    cases = [
        (label="BayesC", name="ib_bayesc", method="BayesC", annotated=false, multitrait=false, sampler=:auto),
        (label="BayesR", name="ib_bayesr", method="BayesR", annotated=false, multitrait=false, sampler=:auto),
        (label="Annotated_BayesC", name="ib_annotated_bayesc", method="BayesC", annotated=true, multitrait=false, sampler=:auto),
        (label="Annotated_BayesR", name="ib_annotated_bayesr", method="BayesR", annotated=true, multitrait=false, sampler=:auto),
        (label="MT_BayesC_I", name="ib_mt_bayesc_i", method="BayesC", annotated=false, multitrait=true, sampler=:I),
        (label="MT_BayesC_II", name="ib_mt_bayesc_ii", method="BayesC", annotated=false, multitrait=true, sampler=:II),
        (label="MT_Annotated_BayesC_I", name="ib_mt_annotated_bayesc_i", method="BayesC", annotated=true, multitrait=true, sampler=:I),
        (label="MT_Annotated_BayesC_II", name="ib_mt_annotated_bayesc_ii", method="BayesC", annotated=true, multitrait=true, sampler=:II),
    ]
    modes = [
        (label="exact_fast_blocks", independent_blocks=false),
        (label="independent_blocks", independent_blocks=true),
    ]

    rows = NamedTuple[]
    mktempdir() do tmpdir
        for case in cases
            for mode in modes
                run_case!(rows, case, mode, genofile, phenotypes, tmpdir)
            end
        end
    end

    summary = DataFrame(rows)
    CSV.write(outfile, summary)
    println("WROTE ", outfile)
    println(summary)
    return summary
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
