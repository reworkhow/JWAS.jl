# Annotated BayesC Standard Probit Implementation

## Summary

Aligned annotated BayesC with standard probit identification and the existing
annotated BayesR annotation semantics.

## Code Changes

Updated [annotation_updates.jl](/Users/haocheng/Github/JWAS.jl/src/1.JWAS/src/MCMC/annotation_updates.jl):

- `sample_annotation_liabilities_single_step!` now samples latent liabilities
  with `latent_sd = 1.0`
- `update_annotation_priors_single_step!` now refreshes `π_j` with the
  standard-normal CDF
- added an inline comment clarifying that `ann.variance` is now reserved for
  coefficient shrinkage in annotated BayesC

Updated [test_annotated_bayesc.jl](/Users/haocheng/Github/JWAS.jl/test/unit/test_annotated_bayesc.jl):

- rewrote the BayesC annotation-semantics test to expect standard-probit latent
  variance
- verified the old implementation failed this test before the code change

Updated [annotated_bayesc.md](/Users/haocheng/Github/JWAS.jl/docs/src/manual/annotated_bayesc.md):

- documented that the annotation latent error variance is fixed to `1`
- documented that the annotation variance parameter controls coefficient
  shrinkage rather than latent noise scale

## Verification

Run:

- `julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'`
- `julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'`
- `julia --project=docs --startup-file=no docs/make.jl`

Attempt:

- `julia --project=. --startup-file=no test/runtests.jl`

## Result

Annotated BayesC and annotated BayesR now share the same standard-probit
annotation identification convention:

- latent annotation error variance fixed to `1`
- annotation variance used for coefficient shrinkage only
