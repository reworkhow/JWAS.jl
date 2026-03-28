# Annotated BayesR Implementation

## Goal

Add dense single-trait annotated BayesR to JWAS as an individual-level analogue of `temp/BayesR_Jian/sbayesrc.R`, with annotation-driven four-class mixture priors and step-specific annotation-coefficient output.

## What Changed

### 1. Step-aware annotation container

`MarkerAnnotations` now supports both:

- `nsteps = 1` for Annotated BayesC
- `nsteps = 3` for annotated BayesR

For annotated BayesR, the container stores:

- `coefficients`, `mean_coefficients`, `mean_coefficients2` as `p x 3`
- `variance` as a length-3 vector
- `liability`, `mu`, `lower_bound`, and `upper_bound` as `m x 3`
- `snp_pi` as `m x 4`

BayesC behavior is preserved by keeping the single-step shapes and single-step update helpers intact.

That is a deliberate v1 stability choice. Annotated BayesC is still routed through dedicated one-step helpers rather than the exact same generalized helper path used for the new three-step BayesR annotation update.

### 2. Annotated BayesR initialization

`get_genotypes(...; method="BayesR", annotations=...)` now:

- accepts annotations for dense single-trait BayesR
- validates the annotation design after JWAS applies genotype QC/filtering
- initializes the three probit intercepts from the default BayesR prior
  - `pi = [0.95, 0.03, 0.015, 0.005]`
  - `p1 = 0.05`
  - `p2 = 0.40`
  - `p3 = 0.25`
- sets all non-intercept annotation coefficients to zero
- initializes `snp_pi` so every marker starts with the same BayesR default prior
- rejects custom annotated BayesR `Pi` vectors that would create degenerate conditional splits and infinite probit intercepts

This gives a common starting prior across markers, while still parameterizing the model through annotation coefficients.

### 3. Annotation Gibbs updates

Annotated BayesR uses three conditional probit models:

- step 1: `delta > 1`
- step 2: `delta > 2`, conditional on `delta > 1`
- step 3: `delta > 3`, conditional on `delta > 2`

Each step has:

- its own liabilities
- its own coefficients
- its own annotation-effect variance

After sampling the three conditional models, JWAS rebuilds the four-class per-marker probabilities:

- `pi_j1 = 1 - p1_j`
- `pi_j2 = p1_j * (1 - p2_j)`
- `pi_j3 = p1_j * p2_j * (1 - p3_j)`
- `pi_j4 = p1_j * p2_j * p3_j`

This is the same annotation-model structure as `sbayesrc.R`. The difference is that JWAS still uses the production individual-level BayesR marker update.

### 4. Dense BayesR prior hookup

The dense BayesR sweep now accepts either:

- one shared four-class prior vector
- one `nMarkers x 4` per-marker prior matrix

When annotations are present, dense BayesR uses row `j` of `ann.snp_pi` as the class prior for marker `j`.

The dense BayesR hot path only performs cheap prior-shape checks at sweep time. The expensive row-by-row validation is kept out of the production marker loop to avoid an O(m) penalty on every annotated BayesR sweep.

### 5. Output and documentation

Annotated BayesR adds:

- `output["annotation coefficients <genotype name>"]`

with columns:

- `Annotation`
- `Step`
- `Estimate`
- `SD`

The step labels are:

- `step1_zero_vs_nonzero`
- `step2_small_vs_larger`
- `step3_medium_vs_large`

The new manual page is:

- `docs/src/manual/annotated_bayesr.md`

It explains the conditional-probability parameterization, default initialization, Gibbs order, current restrictions, and the output schema.

## Mapping To `sbayesrc.R`

The annotated BayesR implementation follows the same statistical structure as `sbayesrc.R`:

- three conditional annotation models for four classes
- step-up indicators derived from `delta`
- reconstruction of joint class probabilities from conditional probabilities
- separate annotation-effect variances by step

The intended difference is computational:

- `sbayesrc.R` uses summary-statistics marker updates
- JWAS uses the existing individual-level BayesR marker equations

## V1 Scope

Included:

- single-trait only
- dense only
- `method="BayesR"` with `annotations=...`
- standard BayesR outputs
- step-specific annotation coefficient output
- documentation page for the sampler

Deferred:

- multi-trait annotated BayesR
- `fast_blocks`
- `storage=:stream`
- random regression models
- annotation-specific phenotype or genetic-value contribution summaries
- extra SNP-level prior-probability output tables

## Verification

Targeted verification completed:

- `julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesc.jl")'`
- `julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'`
- `julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'`
- `julia --project=docs --startup-file=no docs/make.jl`

Full-suite verification:

- `julia --project=. --startup-file=no test/runtests.jl`
- result: `601/601` tests passed in `2m10.3s`
