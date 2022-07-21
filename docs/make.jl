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
        "Mixed Effects Neural Networks (NNMM)" => Any[
                    "Part1.introduction" => "Part1_introduction.md"
                    "Part2.NNMM" => "Part2_NNMM.md"
                    "Part3.NNMM with intermediate omics" => "Part3_NNMMwithIntermediateOmics.md"
                    "Part4.partial connected neural network" => "Part4_PartialConnectedNeuralNetwork.md"
                    "Part5.user-defined onlinear function" => "Part5_UserDefinedNonlinearFunction.md"
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
