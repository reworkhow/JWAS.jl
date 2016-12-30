include("../src/JWAS.jl")
cd("/Users/haocheng/Desktop/hailin-8.reByCECG.123only.aod")

using DataFrames
using JWAS

ped = get_pedigree("finalPED");

dat = readtable("finalDAT", separator = ' ');

model_equations = "bw    = intercept + bwcg + HH + XX + aod + animal + dam;
                   ce    = intercept + cecg + HH + XX + aod + animal + dam";

R = readdlm("startingValues.R");
G = readdlm("startingValues.G");
starting = readdlm("startingValues.sol");

model = build_model(model_equations,R);
set_covariate(model,"HH XX")
set_random(model,"animal dam", ped, G)

#@time out=solve(model, dat,solver="GaussSeidel", printout_frequency=100, tolerance = 0.0001);
#starting=map(Float64, out[:, 2]);
#writedlm("startingValues.sol", starting)

@time outMCMC=runMCMC(model, dat, chain_length=2, printout_frequency=1) #, starting_value=starting);

@profile outMCMC=runMCMC(model, dat, chain_length=2, printout_frequency=1)
