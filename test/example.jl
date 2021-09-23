using Revise
using JWAS
using DataFrames, CSV, LinearAlgebra, Kronecker, Statistics;
using GSL;

data_path = "/Users/apple/Downloads/"
geno_file = data_path * "genoRRM";

# generate dataframe having Day columns
RR_pheno = data_path * "phenoRRM.csv"
try
    mkdir(data_path * "RR")
catch
    "a folder has been existed."
end

pheno = CSV.read(RR_pheno, DataFrame)

ncoef = 2
tp = length(unique(pheno[!, :time]))

### Generate Lengendre Polynomial Covariates
function generatefullPhi(timevec, ncoeff = 3)
    # timevec: a vector of time points for all observations
    # ncoeff : the number of polynomial coefficients to be estimated
    times = sort(unique(timevec)) # vector of whole timepoints sorted
    tmin = minimum(times)
    tmax = maximum(times)
    # standardized time points
    qi = 2 * (times .- tmin) ./ (tmax .- tmin) .- 1
    Φ = Matrix{Float64}(undef, length(times), ncoeff) # Phi of size ntimes x ncoeff
    for i = 1:ncoeff
        n = i - 1
        ## Given the normalized function from L. R. Schaeffer
        Φ[:, i] = sqrt((2 * n + 1) / 2) * sf_legendre_Pl.(n, qi) # construct the matrix
    end
    return Φ
end

myPhi = generatefullPhi(pheno[!, :time], ncoef);


cd(data_path * "RR")


genotypes = get_genotypes(
    geno_file,
    separator = ',',
    header = true,
    method = "BayesC",
    estimatePi = true,
);

model_equation = "yobs = intercept + genotypes"

model = build_model(model_equation);

@time output = runMCMC(
    model,
    pheno,
    RRM = myPhi,
    chain_length = 5000,
    outputEBV = true,
    burnin = 2000,
    double_precision = true,
)

output["marker effects genotypes"]

#### Validation (sample codes, users need to do the training/testing split)

u1 = output["EBV_1"] # estimated first random regression coefficients for individuals
u2 = output["EBV_2"] # estimated second random regression coefficients for individuals
EBVs = DataFrame(
    ID = u1[!, :ID],
    u1 = u1[!, :EBV],
    u2 = u2[!, :EBV],
    #u3 = u3[!, :EBV],
)

# Calculate the BV for individuals
n = size(EBVs, 1)
Idt = Matrix{Float64}(I, n, n)
Phi = generatefullPhi(collect(1:tp), ncoef)
Z = Idt ⊗ Phi
us = [[] for i = 1:n]
for i = 1:n
    us[i] = Matrix(EBVs)[i, 2:end]
end
uhat = vcat(us...)
gBLUP = Z * uhat
gBLUP = convert(Array{Float64,1}, gBLUP)

correctEBVs = DataFrame(
    ID = repeat(EBVs[!, :ID], inner = tp),
    time = repeat(collect(1:tp), outer = n),
    EBV = gBLUP,
)

first(correctEBVs, 10)

accuracy = []
for timei = 1:tp
    println("Day$timei.")
    phenoi = filter(row -> row.time in [timei], copy(pheno))
    select!(phenoi, Not(:time))
    EBVi = filter(row -> row.time in [timei], copy(correctEBVs))
    select!(EBVi, Not(:time))
    final_dfi = innerjoin(phenoi, EBVi, on = :ID)
    res = cor(final_dfi[!, :EBV], final_dfi[!, :yobs])
    push!(accuracy, res)
end

accuracy
