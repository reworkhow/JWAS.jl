# Long test: Single-step Bayesian Regression with block updates
# Uses simulated_omics data (3,534 genotyped / 6,473 pedigree animals)

using JWAS, DataFrames, CSV, Statistics, JWAS.Datasets, Random, Test

# Load data
phenofile = Datasets.dataset("phenotypes_sim.txt", dataset_name="simulated_omics")
pedfile   = Datasets.dataset("pedigree.txt", dataset_name="simulated_omics")
genofile  = Datasets.dataset("genotypes_1000snps.txt", dataset_name="simulated_omics")

phenotypes = CSV.read(phenofile, DataFrame, delim=',', header=true, missingstring=["NA"])
pedigree   = get_pedigree(pedfile, separator=",", header=true)

# Build single-step phenotype data (add non-genotyped pedigree animals)
pedfile_df = CSV.read(pedfile, DataFrame, delim=',', header=true)
geno_ids = Set(string.(phenotypes[!, :ID]))
non_geno_ids = setdiff(Set(string.(pedfile_df[!, :ID])), geno_ids)

extra_ids = collect(non_geno_ids)[1:min(500, length(non_geno_ids))]
Random.seed!(42)
extra_pheno = DataFrame(
    ID = extra_ids,
    trait1 = randn(length(extra_ids)),
    genetic_total = zeros(length(extra_ids))
)
phenotypes_ss = vcat(
    select(phenotypes, [:ID, :trait1, :genetic_total]),
    extra_pheno
)
phenotypes_ss[!, :ID] = string.(phenotypes_ss[!, :ID])

# Get genotypes
geno = get_genotypes(genofile, 1.0, separator=',', method="BayesC")

# Build model
model = build_model("trait1 = intercept + geno", 1.0)

# Run single-step with block updates
@time out = runMCMC(model, phenotypes_ss,
                    single_step_analysis=true,
                    pedigree=pedigree,
                    chain_length=1000,
                    fast_blocks=10,
                    seed=314)

# Check accuracy
results = innerjoin(out["EBV_trait1"], phenotypes, on=:ID)
accuracy = cor(Float64.(results[!, :EBV]), Float64.(results[!, :genetic_total]))
println("Single-step block BayesC prediction accuracy: $(round(accuracy, digits=3))")
@test accuracy > 0.0
