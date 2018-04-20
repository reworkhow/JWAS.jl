using Documenter, JWAS

makedocs(
#    modules = [Documenter],
    doctest=false,
    clean  =true,
    format = :html,
#    assets = ["assets/favicon.ico"],
    sitename = "JWAS.jl",
    authors = "Hao Cheng, Rohan Fernando, Dorian Garrick and contributors.",
    pages = [
        "Home" => "index.md",
        "Some Theory" => "theory/theory.md",
        "Manual" => Any[
            "Guide" => "man/guide.md",
            "man/examples.md",
            "Contributing" => "examples/genomicBLMM.md",
            ],

        ],
        "Examples" => Any[
            "Linear Mixed Model (conventional)" => "examples/conventionalBLMM.md",
            "Linear Additive Genetic Model" => "examples/LinearAdditiveGeneticModel.md",
            "Linear Mixed Model (Genomic data)" => "examples/genomicBLMM.md",
        ],
        "Library" => Any[
            "Public" => "lib/public.md",
            "Internals" => "lib/internals.md",
            # hide("Internals" => "lib/internals.md", Any[
            #     "lib/internals/anchors.md",
            #     "lib/internals/builder.md",
            #     "lib/internals/cross-references.md",
            #     "lib/internals/docchecks.md",
            #     "lib/internals/docsystem.md",
            #     "lib/internals/documenter.md",
            #     "lib/internals/documents.md",
            #     "lib/internals/dom.md",
            #     "lib/internals/expanders.md",
            #     "lib/internals/formats.md",
            #     "lib/internals/generator.md",
            #     "lib/internals/mdflatten.md",
            #     "lib/internals/selectors.md",
            #     "lib/internals/textdiff.md",
            #     "lib/internals/utilities.md",
            #     "lib/internals/walkers.md",
            #     "lib/internals/writers.md",
            #])
        ]
    ],
    # Use clean URLs, unless built as a "local" build
    #html_prettyurls = !("local" in ARGS),
    #html_canonical = "https://juliadocs.github.io/Documenter.jl/stable/",
)

deploydocs(
    repo = "github.com/reworkhow/JWAS.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
    julia  = "0.6",
    osname = "osx"
)
