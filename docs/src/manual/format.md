# Public functions

Documentation for JWAS.jl's public interface. Below are functions available to general users.


## Input

```@index
Pages = ["public.md"]
Modules = [JWAS]
```

```@index
Pages = ["public.md"]
Modules = [JWAS.misc]
```


## Output file

MCMC samples for marker effects

each row is one sample for all marker effects
the number of rows is equal to the number of samples


MCMC samples for genetic or residual covariance matrix

each row is one sample for  vec(covariances)
the number of rows is equal to the number of samples
