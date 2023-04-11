using Revise
using JWAS
using DataFrames, CSV, Statistics, HTTP;
using LinearAlgebra, Kronecker;
using Random
using GSL;
using Distributions;
using DelimitedFiles;
using Plots;

function getdata(file_name)  # function to load data from github folder
    http_obj = HTTP.get("https://raw.githubusercontent.com/Jiayi-Qu/bio_protocol/main/data/$file_name")
    data = CSV.read(http_obj.body, DataFrame, header=true, missingstrings=[""])
    return data
end
############################################################################################################
# whole model with real phenotypes (example)
############################################################################################################

geno_file = getdata("Genotype2_first5k.csv")
pheno = getdata("RR_BLUE_missing2_PSA.csv")
pheno[!, :PE] = pheno[!, :ID]
pheno

ncoef = 3
tp = length(unique(pheno[!, :time]));

### Generate Lengendre Polynomial Covariates
function generatefullPhi(timevec, ncoeff=3)
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
# the Ledengre polynomial matrix 
myPhi = generatefullPhi(pheno[!, :time], ncoef)

# split data into training and testing data set 
# number of testing individual = 119
IDs = unique(pheno[!, :ID])
nind = size(IDs, 1)
nt = 119
testingInd = sample(IDs, nt, replace=false)
trainingInd = filter!(x -> x ∉ testingInd, IDs);
pheno_training = filter(row -> row.ID in trainingInd, copy(pheno))
pheno_testing = filter(row -> row.ID in testingInd, copy(pheno));


# add corresponding Phi1 and Phi2 columns to the pheno_training data frame
nrec = size(pheno_training, 1)
wholePhi = Array{Float64}(undef, nrec, ncoef);
r = 1
for timei in pheno_training[!, :time]
    wholePhi[r, :] = myPhi[timei, :]
    global r += 1
end

for i in 1:ncoef
    pheno_training[!, Symbol("Phi$i")] = wholePhi[:, i]
end

# get_genotypes
geno = get_genotypes(geno_file, header=true, separator=',', method="BayesC")

#model equation 
fixed_effect = join(" + Phi" .* string.(collect(2:ncoef)))
pe_effect = join(" + Phi" .* string.(collect(1:(ncoef-1))) .* "*PE")
residual_effect = join(" + Day" .* string.(collect(1:20)))
model_equation = "PSA = Phi1" * fixed_effect * pe_effect * residual_effect * " + geno"

model = build_model(model_equation);

# set covariates
for i in 1:ncoef
    set_covariate(model, "Phi$i")
end
# set random effects
for i in 1:(ncoef-1)
    set_random(model, "Phi$i*PE")
end

for d in 1:tp
    set_random(model, "Day$d")
end


# run MCMC sampling 
@time output = runMCMC(outputEBV=true, model, pheno_training, burnin=100,
    chain_length=1000, RRM=myPhi, seed=123);


# estimated EBVs (coefficients) for individuals
u1 = CSV.read("results/EBV_1.txt", DataFrame)
u2 = CSV.read("results/EBV_2.txt", DataFrame)
u3 = CSV.read("results/EBV_3.txt", DataFrame)
EBVs = DataFrame(
    ID=u1[!, :ID],
    u1=u1[!, :EBV],
    u2=u2[!, :EBV],
    u3=u3[!, :EBV])
# Calculate the BV for individuals by time 
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
    ID=repeat(EBVs[!, :ID], inner=tp),
    time=repeat(collect(1:tp), outer=n),
    EBV=gBLUP)

# choose EBV for testing Individuals
testingEBVs = filter(row -> row.ID in testingInd, correctEBVs)

# calculate prediction accuracy across time for testing individuals
accuracy = []
for timei = 1:tp
    println("Day$timei.")
    pheno_testingi = filter(row -> row.time in [timei], pheno_testing)
    select!(pheno_testingi, Not(:time))
    EBV_testingi = filter(row -> row.time in [timei], testingEBVs)
    select!(EBV_testingi, Not(:time))
    final_dfi = innerjoin(pheno_testingi, EBV_testingi, on=:ID)
    res = cor(final_dfi[!, :EBV], final_dfi[!, :PSA])
    push!(accuracy, res)
end
accuracy

############################################################################################################
# simple model with simulated data 
############################################################################################################
pheno = getdata("RR_PSA_rep1.csv")[!, 2:end]
ncoef = 2
tp = length(unique(pheno[!, :time]));
myPhi = generatefullPhi(pheno[!, :time], ncoef)

# split data into training and testing data set 
# number of testing individual = 119
IDs = unique(pheno[!, :ID])
nind = size(IDs, 1)
nt = 119
testingInd = sample(IDs, nt, replace=false)
trainingInd = filter!(x -> x ∉ testingInd, IDs);
pheno_training = filter(row -> row.ID in trainingInd, copy(pheno))
pheno_testing = filter(row -> row.ID in testingInd, copy(pheno));

# add corresponding Phi1 and Phi2 columns to the pheno_training data frame
nrec = size(pheno_training, 1)
wholePhi = Array{Float64}(undef, nrec, ncoef);
r = 1
for timei in pheno_training[!, :time]
    wholePhi[r, :] = myPhi[timei, :]
    global r += 1
end

for i in 1:ncoef
    pheno_training[!, Symbol("Phi$i")] = wholePhi[:, i]
end

# get_genotypes
geno = get_genotypes(geno_file, header=true, separator=',', method="BayesC")

#model equation 
fixed_effect = join(" + Phi" .* string.(collect(2:ncoef)))
model_equation = "PSA = Phi1" * fixed_effect * "+ geno"

model = build_model(model_equation);

# set covariates
set_covariate(model, "Phi1")
set_covariate(model, "Phi2");

# run MCMC sampling 
@time output = runMCMC(outputEBV=true, model, pheno_training, burnin=100,
    chain_length=1000, RRM=myPhi, seed=123);

#### Validation (sample codes)
# estimated EBVs (coefficients) for individuals
u1 = CSV.read("results1/EBV_1.txt", DataFrame)
u2 = CSV.read("results1/EBV_2.txt", DataFrame)
EBVs = DataFrame(
    ID=u1[!, :ID],
    u1=u1[!, :EBV],
    u2=u2[!, :EBV])
# Calculate the BV for individuals by time 
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
    ID=repeat(EBVs[!, :ID], inner=tp),
    time=repeat(collect(1:tp), outer=n),
    EBV=gBLUP)

# choose EBV for testing Individuals
testingEBVs = filter(row -> row.ID in testingInd, correctEBVs)

# calculate prediction accuracy across time for testing individuals
accuracy = []
for timei = 1:tp
    println("Day$timei.")
    pheno_testingi = filter(row -> row.time in [timei], pheno_testing)
    select!(pheno_testingi, Not(:time))
    EBV_testingi = filter(row -> row.time in [timei], testingEBVs)
    select!(EBV_testingi, Not(:time))
    final_dfi = innerjoin(pheno_testingi, EBV_testingi, on=:ID)
    res = cor(final_dfi[!, :EBV], final_dfi[!, :PSA])
    push!(accuracy, res)
end
accuracy


# simulated true QTL
trueQTLs = getdata("QTLpos_rep1.csv")[:, 2:end]

# RR MF & WPPF
# calculate the model frequency for coef i (1: intercept, 2: slope)
coef = 2
println("Coef $coef")
mrk_eff_sample_c = "results1/MCMC_samples_marker_effects_geno_$coef.txt"
mrk_samples, markerID = readdlm(mrk_eff_sample_c, header=true, ',')
markerID = vec(strip.(markerID, ['\"']))
model_freq = vec(mean(mrk_samples .!= 0, dims=1))
# markerOrder is used to order SNPs in the plots
RR_MFc = DataFrame(markerID=markerID, modelfrequency=model_freq, markerOrder=collect(1:length(markerID)))

# generate plots for model frequency
trueQTL_union = vec(trueQTLs[!, Symbol("coef$coef")])
QTL_MF = filter(row -> row.markerID in unique(trueQTL_union), copy(RR_MFc))
noneQTLc = filter!(x -> x ∉ unique(trueQTL_union), copy(RR_MFc[!, :markerID]))
SNP_MF = filter(row -> row.markerID in noneQTLc, copy(RR_MFc));

plot(SNP_MF[!, :markerOrder], SNP_MF[!, :modelfrequency], seriestype=:scatter, legend=:none)
p = plot!(QTL_MF[!, :markerOrder], QTL_MF[!, :modelfrequency], seriestype=:scatter, legend=:none,
    xlabel="SNP", ylabel="MF")
    #savefig(p,RRM_path*"RR_MF_coef$coef.png")