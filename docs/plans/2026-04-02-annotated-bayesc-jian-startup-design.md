# Annotated BayesC Jian-Style Startup Design

## Goal

Align dense annotated BayesC startup with Jian's `sbayesrc.R` startup structure:

- start from the supplied common prior
- start annotation coefficients at zero
- do not fit the annotation model before the first phenotype-informed marker sweep

The immediate reason is seed stability. The earlier BayesC startup path fit the annotation model from synthetic startup labels before the sampler had seen phenotype information, which created seed-specific annotation states.

## Design

### 1. Startup prior

For annotated BayesC, JWAS should convert the supplied scalar `Pi` into a
marker-level exclusion vector:

`π_j = Pi` for all markers `j`.

If the caller already supplies a marker-level `Pi` vector, JWAS should preserve it
after validating the length.

### 2. Startup annotation state

The annotation state should start at:

- coefficients `α = 0`
- linear predictor `μ = 0`

This intentionally does **not** force the annotation model to reproduce the
starting `π_j` yet. The first marker sweep should use `Mi.π` directly, and only
after that should the annotation model be fit from phenotype-informed marker
indicators.

This matches Jian's structure more closely than the earlier JWAS intercept-based
initialization.

### 3. MCMC startup flow

Dense annotated BayesC should skip the old pre-MCMC annotation initialization and
pre-fit. The first annotation update should happen only inside the normal
`estimatePi` section after the first marker sweep.

Annotated BayesR startup is unchanged in this pass.

### 4. Marker-level priors

Annotated BayesC now starts with marker-level `π_j`, so the BayesC variance-prior
setup must accept a length-`nMarkers` `Pi` vector and interpret it as marker-level
exclusion probabilities instead of BayesR class weights.

## Expected Outcome

If the unstable startup pre-fit was the main source of seed sensitivity, dense
annotated BayesC should become much more stable across seeds on the large
reproducer.

## Scope

This pass is intentionally narrow:

- dense annotated BayesC only
- no fast-block changes
- no BayesR algorithm changes
- no annotation-regularization redesign
