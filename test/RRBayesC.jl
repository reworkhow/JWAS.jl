using Revise
using DataFrames,CSV,JWAS,Statistics,DelimitedFiles,Random,Distributions,LinearAlgebra;

data_path = "/Users/apple/Desktop/JiayiQu/UCD_PhD/RRM/WU/data/"
analysis_path = "/Users/apple/Box/Jiayi_Computer/UCD_PhD/RRM/peer_review/DIC/"
ncoef = 3;

using GSL
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

RR_pheno = data_path*"RR_BLUE_missing2_PSA.csv"
pheno = CSV.read(RR_pheno, DataFrame);

geno_file = data_path*"Genotype2.csv";
Random.seed!(123);

myPhi=generatefullPhi(pheno[!,:time],ncoef);

cd(analysis_path*"AllSamples")

nrec = size(pheno,1)
wholePhi = Array{Float64}(undef, nrec, ncoef);

r = 1
for timei in pheno[!, :time]
    wholePhi[r,:] = myPhi[timei,:]
    global r += 1
end

for i in 1:ncoef
    pheno[!,Symbol("Phi$i")] = wholePhi[:,i]
end

residual_var = readdlm("/Users/apple/Box/Jiayi_Computer/UCD_PhD/RRM/peer_review/DominanceEffect/data/ST_RRBLUP_rvar.txt");

geno = get_genotypes(geno_file, separator = ',', method = "BayesC", estimatePi = true);

# construct the model
hetero_residual = join(" + Day" .* string.(collect(1:20)))
fixed_effect = join(" + Phi" .* string.(collect(2:ncoef)))
pe = join(" + Phi" .* string.(collect(1:ncoef)) .* "*ID") # linear regression for PE (dominance effects )
model_equation = "PSA = Phi1" * fixed_effect * pe * hetero_residual * " + geno"
model = build_model(model_equation);

for i in 1:ncoef
    set_covariate(model, "Phi$i")
end

for i in 1:ncoef
    set_random(model, "Phi"*string(i)*"*ID")
end

for d in 1:20
    set_random(model, "Day$d", residual_var[d])
end

@time output=runMCMC(model,pheno,
    chain_length=1000,burnin=0,RRM=myPhi,double_precision=true,
    output_samples_frequency = 50, output_samples_for_all_parameters = true);
