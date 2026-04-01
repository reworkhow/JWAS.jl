# Annotated BayesC Standard Probit Design

## Goal

Correct the annotated BayesC annotation submodel so its latent probit error
variance is fixed to `1`, matching standard probit identification and the
existing annotated BayesR implementation.

## Problem

Annotated BayesC previously used one scalar `ann.variance` for two different
roles:

- latent annotation error variance in the liability draw
- prior variance for annotation coefficients in the Gibbs update

That is not a clean model parameterization. In a probit model, latent error
variance should be fixed for identification, while coefficient shrinkage should
remain a separate quantity.

## Design

Keep the current one-step BayesC annotation architecture, but change its
semantics:

- latent liabilities are sampled from `N(mu, 1)` truncated by the current
  inclusion indicator
- the BayesC exclusion prior refresh uses the standard-normal CDF:
  `π_j = 1 - Φ(mu_j)`
- `ann.variance` remains in use only for the annotation-coefficient prior /
  shrinkage path

## Scope

In scope:

- annotated BayesC one-step annotation updater
- BayesC annotation unit tests
- user documentation and implementation notes

Out of scope:

- BayesR annotation logic
- phenotype residual variance updates
- redesign of the `MarkerAnnotations` container
