using Revise
using JWAS
using DataFrames,CSV,LinearAlgebra,Kronecker,Statistics;
using GSL;

data_path = "/Users/apple/Desktop/JiayiQu/UCD_PhD/RRM/WU/simulation/wholegenome/example_data/data/"
geno_file = data_path * "Mqtl.csv";

# generate dataframe having Day columns
RR_pheno = data_path*"RR_PSA_sim.csv"
try
    mkdir(data_path*"RR")
catch
    "a folder has been existed."
end

pheno = CSV.read(RR_pheno, DataFrame)[:,2:end]

ncoef = 3
tp = length(unique(pheno[!,:time]))

function generatefullPhi(timevec, ncoeff=3)
    # timevec: a vector of time points for all observations
    # ncoeff : the number of polynomial coefficients to be estimated
    times = sort(unique(timevec)) # vector of whole timepoints sorted
    tmin  = minimum(times)
    tmax  = maximum(times)
    # standardized time points
    qi = 2*(times .- tmin)./(tmax .- tmin) .- 1
    Φ  = Matrix{Float64}(undef, length(times), ncoeff) # Phi of size ntimes x ncoeff
    for i in 1:ncoeff
        n = i - 1
        ## Given the normalized function from L. R. Schaeffer
        Φ[:,i] = sqrt((2*n+1)/2)*sf_legendre_Pl.(n,qi) # construct the matrix
    end
    return Φ
end

myPhi=generatefullPhi(pheno[!,:time], ncoef);


cd(data_path*"RR")


genotypes  = get_genotypes(geno_file,separator=',', header = true, method="BayesC", estimatePi =true);

model_equation = "PSA = intercept + genotypes"#* fixed_effect * PE_effect * hetero_residual



model = build_model(model_equation);
#=
for i in 1:ncoef
    set_covariate(model, "Phi$i")
end

for i in 1:(ncoef-1)
    set_random(model, "Phi$i*ID")
end

PE_cov = join("Phi" .* string.(collect(1:(ncoef-1))) .* "*ID ")

Kpe = [1 0.5;0.5 1]

set_random(model, PE_cov, Kpe)

vare = collect(1:20).^2 # usually get from ST analysis

for d in 1:20
    set_random(model, "Day$d", vare[d])
end
=#
@time output=runMCMC(model,pheno, RRM=myPhi, chain_length = 5000, outputEBV=true, burnin = 3000, double_precision = true )

output["marker effects genotypes"]
#### Cross validation
path = data_path * "RR/results1/"

u1 = CSV.read(path*"EBV_1.txt", DataFrame) # estimated first random regression coefficients for individuals
u2 = CSV.read(path*"EBV_2.txt", DataFrame) # estimated second random regression coefficients for individuals
u3 = CSV.read(path*"EBV_3.txt", DataFrame) # estimated third random regression coefficients for individuals
EBVs = DataFrame(ID = u1[!,:ID], u1 = u1[!,:EBV], u2 =u2[!,:EBV], u3 =u3[!,:EBV])

# Calculate the BV for individuals across 20 days
n = size(EBVs,1)
Idt = Matrix{Float64}(I, n, n)
Phi=generatefullPhi(collect(1:tp),ncoef)
Z = Idt ⊗ Phi
us = [[] for i in 1:n]
for i in 1:n
    us[i] = Matrix(EBVs)[i,2:end]
end
uhat = vcat(us...)
gBLUP = Z*uhat
gBLUP= convert(Array{Float64,1}, gBLUP)

correctEBVs = DataFrame(ID = repeat(EBVs[!,:ID], inner = tp),
    time = repeat(collect(1:tp),outer=n),EBV = gBLUP)

first(correctEBVs, 10)

accuracy = []
for timei in 1:20
    println("Day$timei.")
    phenoi = filter(row -> row.time in [timei], copy(pheno))
    select!(phenoi,Not(:time))
    EBVi = filter(row -> row.time in [timei], copy(correctEBVs))
    select!(EBVi,Not(:time))
    final_dfi = innerjoin(phenoi, EBVi, on = :ID)
    println(size(final_dfi))

    res = cor(final_dfi[!,:EBV], final_dfi[!,:PSA])
    push!(accuracy, res)
end

accuracy
