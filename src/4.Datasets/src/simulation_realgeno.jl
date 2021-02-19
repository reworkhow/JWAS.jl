# A dataset of phenotype, genotype, and pedigree is simulated below. This dataset is used in the following examples. Genotype and pedigree data are simulated using XSim. Phenotypic data is simulated below based on
#
# ```
# y1 = intercept + x1 + x2 + x2*x3 + ID + dam + geno
# y2 = intercept + x1 + x2 + ID + geno
# y3 = intercept + x1 + ID + geno
#
# covariate:x1
# random: x2 (i.i.d), "ID","dam" (pedigree)
# ```

using DelimitedFiles,JWAS, Statistics, Distributions, DataFrames, LinearAlgebra,CSV
using Random
Random.seed!(314);

#load genotype and pedigree
M            = readdlm("genotypes.csv",',');
pedigree     = CSV.read("pedigree.csv",DataFrame)
ped          = get_pedigree("pedigree.csv",header=true);
IDs,Ai,inbred= get_info(ped,Ai=true);
A            = Symmetric(inv(Matrix(Ai))); #numerator relationship matrix

#data simulation
data      = DataFrame(ID=M[2:end,1])
data[!,:x1] = round.(randn(size(data,1)),digits=2)
data[!,:x2] = sample(["g1","g2","g3","g4","g5"],size(data,1))
data[!,:x3] = sample(["m","f"],size(data,1));
data      = leftjoin(data, pedigree[!,[:ID,:dam]], on = :ID)
data

#Genetic/Phenotypic Value Simulation
##simulation parameters
nQTL = 10
h2   = [0.7,0.5,0.2];

##genomic info
QTL_covariates = M[2:end,rand(2:size(M,2),nQTL)];
QTL_effects    = rand(MvNormal(zeros(3), [1.0 0.5 0.3
                                          0.5 1.0 0.5
                                          0.3 0.5 1.0]),nQTL);

BV = round.(QTL_covariates*QTL_effects',digits=2);
bv1 = BV[:,1]./std(BV[:,1])*sqrt(h2[1])
bv2 = BV[:,2]./std(BV[:,2])*sqrt(h2[2])
bv3 = BV[:,3]./std(BV[:,3])*sqrt(h2[3])
println("variance for bv1, bv2, and bv3 are: $(var(bv1)), $(var(bv2)), $(var(bv3))")
##maternal effects
dam_effect = rand(MvNormal(zeros(length(IDs)), A*sqrt(0.1)))
dam        = JWAS.mkmat_incidence_factor(data[!,:dam],[IDs;missing])*[dam_effect;0.0];

##non-genetic effects
x1_effect,x2_effect,x3_effect  = 1.0, randn(5), randn(2)
x1 = data[!,:x1]*x1_effect
x2 = JWAS.mkmat_incidence_factor(data[!,:x2],["g1","g2","g3","g4","g5"])*x2_effect
x3 = JWAS.mkmat_incidence_factor(data[!,:x3],["m","f"])*x3_effect;
##genetic/phenotypic value
data[!,:bv1] = bv1 + dam
data[!,:bv2] = bv2
data[!,:bv3] = bv3;

data[!,:y1] = data[!,:bv1] + x1 + x2 + x2.*x3 + randn(size(data,1))*(1-h2[1]-0.1)
data[!,:y2] = data[!,:bv2] + x1 + x2 + randn(size(data,1))*(1-h2[2])
data[!,:y3] = data[!,:bv3] + x1 + randn(size(data,1))*(1-h2[3]);

data=data[!,[:ID,:y1,:y2,:y3,:x1,:x2,:x3,:dam,:bv1,:bv2,:bv3]]

CSV.write("phenotypes.csv",data,missingstring="NA")
