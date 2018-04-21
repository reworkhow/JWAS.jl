# Bayesian Linear Mixed Models (Genomic Data)

## Univariate Linear Mixed Model (Genomic data)


## Multivariate Linear Mixed Model (Genomic data)


```@example 1
using DataFrames,CSV,JWAS,JWAS.Datasets;
```


```@example 1
phenofile = Datasets.dataset("testMME","data.txt");
data      = CSV.read(phenofile,delim = ',',header=true);
head(data);
```


```@example 1
model_equation    = "nwn = intercept +parity + parity*site + yr + geneticCode + age"

residual_variance = 2.97;
model             = build_model(model_equation,residual_variance)

set_covariate(model,"age");

geneticCode_variance = 0.26;
set_random(model,"geneticCode",geneticCode_variance);
```

```@example 1
outputMCMCsamples(model,"parity","age");
```


```@example 1
out=runMCMC(data,model,chain_length=50000,output_samples_frequency=100)
```
