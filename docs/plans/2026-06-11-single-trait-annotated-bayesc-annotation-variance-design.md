# Single-Trait Annotated BayesC Annotation Variance Design

## Goal

Sample the annotation slope prior variance in single-trait annotated BayesC,
matching the Jian SBayesRC update already used by JWAS nested annotation
samplers.

## Background

Single-trait annotated BayesC uses one binary probit model for marker inclusion:

```text
delta_j = 1(l_j > 0)
l_j = W_j * alpha + e_j
e_j ~ N(0, 1)
```

The latent probit variance remains fixed at `1` for identifiability. The
annotation-effect variance is a different component: it is the prior variance
for annotation slopes,

```text
alpha_k ~ N(0, sigma_alpha^2), k > 1
```

Currently, single-trait BayesC uses fixed `ann.variance` for slope shrinkage.
Jian's `sbayesrc.R` samples this variance after coefficient updates:

```r
sigmaSqAlpha[j] = (sum(alpha[-1,j]^2) + 2) / rchisq(1, nAnno-1+2)
```

The intercept is excluded.

## Proposed Change

After sampling the single-trait annotated BayesC probit coefficients:

1. If there are annotation slopes, sample scalar `ann.variance` as:

   ```text
   ann.variance = (sum(alpha_slopes^2) + 2) / chisq(number_of_slopes + 2)
   ```

2. If the model is intercept-only, leave `ann.variance` unchanged.

3. Continue rebuilding marker-specific `Pi_j` from the updated probit mean:

   ```text
   Pi_j = 1 - Phi(W_j * alpha)
   ```

## Tests

Add unit coverage that:

- constructs a small single-trait annotated BayesC state with one slope
- samples liabilities and coefficients from a fixed seed
- then expects `ann.variance` to be sampled from the same random stream using
  the Jian formula
- confirms intercept-only annotated BayesC does not change scalar
  `ann.variance`

## Validation

Run:

```bash
julia --project=. --startup-file=no test/unit/test_annotated_bayesc.jl
```
