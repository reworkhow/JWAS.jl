# Bayesian Linear Mixed Models

## Univariate Linear Mixed Model (conventional)

```julia
using DataFrames,CSV,JWAS,JWAS.Datasets;
```


```julia
phenofile = Datasets.dataset("testMME","data.txt");
data      = CSV.read(phenofile,delim = ',',header=true);
head(data)

# 10×9 DataFrames.DataFrame
# │ Row │ sow       │ site │ yr   │ age │ geneticCode │ parity │ nwn │ SYS           │ bw   │
# ├─────┼───────────┼──────┼──────┼─────┼─────────────┼────────┼─────┼───────────────┼──────┤
# │ 1   │ 100-113   │ 113  │ 2005 │ 18  │ P1          │ 1      │ 8   │ 113_2005_WNTR │ 9.0  │
# │ 2   │ 100-113   │ 113  │ 2006 │ 18  │ P1          │ 2      │ 12  │ 113_2006_SPNG │ 8.0  │
# │ 3   │ 100-5     │ 5    │ 2008 │ 15  │ P2          │ 1      │ 10  │ 5_2008_ATMN   │ 7.5  │
# │ 4   │ 1000-5    │ 5    │ 2009 │ 17  │ P2          │ 1      │ 10  │ 5_2009_SPNG   │ 8.3  │
# │ 5   │ 10000-131 │ 13   │ 2004 │ 16  │ Commercial  │ 1      │ 9   │ 13_2004_WNTR  │ 4.3  │
# │ 6   │ 10000-131 │ 13   │ 2004 │ 18  │ Commercial  │ 2      │ 10  │ 13_2004_SMMR  │ 2.8  │
```


```julia
model_equation    = "nwn = intercept +parity + parity*site + yr + geneticCode + age"

residual_variance = 2.97;
model             = build_model(model_equation,residual_variance)

set_covariate(model,"age");

geneticCode_variance = 0.26;
set_random(model,"geneticCode",geneticCode_variance);
```

```julia
outputMCMCsamples(model,"parity","age");
```


```julia
out=runMCMC(data,model,chain_length=50000,output_samples_frequency=100);
```

## Multivariate Linear Mixed Model (conventional)
