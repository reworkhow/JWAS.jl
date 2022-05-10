using Documenter, JWAS

makedocs(
    modules = [JWAS,JWAS.Datasets,JWAS.PedModule],
    doctest=false,
    clean  =true,
#    format = Documenter.HTML(analytics = "UA-90474609-3",
#                             assets = ["assets/favicon.ico"],
#                             canonical="https://docs.juliadiffeq.org/stable/"),
    sitename = "JWAS.jl",
    authors = "Hao Cheng, Rohan Fernando, Dorian Garrick and contributors.",
    pages = Any[
        "Home" => "index.md",
        "Some Theory" => "theory/theory.md",
        "Citing" => "citing/citing.md",
        "Manual" => Any[
            "Get Started" => "manual/getstarted.md",
            "Workflow" => "manual/workflow.md",
            "Public" => "manual/public.md",
            "Internals" => "manual/internals.md"
            ],
        "Examples" => Any[
            "Examples" => "examples/examples.md"
        ],
        "Frequently Asked Questions" => "FrequentlyAskedQuestions/FrequentlyAskedQuestions.md",

    ],
)

deploydocs(
    repo = "github.com/reworkhow/JWAS.jl.git",
    target = "build",
)
