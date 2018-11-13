using Documenter, JWAS

makedocs(
    modules = [JWAS,JWAS.misc,JWAS.Datasets,JWAS.PedModule],
    doctest=false,
    clean  =true,
    format = :html,
#    assets = ["assets/favicon.ico"],
    sitename = "JWAS.jl",
    authors = "Hao Cheng, Rohan Fernando, Dorian Garrick and contributors.",
    pages = [
        "Home" => "index.md",
        "Some Theory" => "theory/theory.md",
        "Citing" => "citing/citing.md",
        "Manual" => Any[
            "Get Started" => "manual/getstarted.md",
            "Workflow" => "manual/workflow.md",
            "Public" => "manual/public.md",
            "Internals" => "manual/internals.md",
            ],
        "Examples" => Any[
            "Examples" => "examples/examples.md",
        ],
    ],
)

deploydocs(
    repo = "github.com/reworkhow/JWAS.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
    julia  = "1.0",
    osname = "osx"
)
