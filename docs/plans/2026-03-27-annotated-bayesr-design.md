# Annotated BayesR Design

## Goal

Add an annotation-aware single-trait dense BayesR to JWAS as an individual-level
analogue of `temp/BayesR_Jian/sbayesrc.R`.

## Scope

Included in v1:

- single-trait only
- dense storage only
- `method="BayesR"` with `annotations=...`
- annotation-controlled full 4-class BayesR mixture prior for each SNP
- standard BayesR outputs
- step-specific annotation coefficient output
- manual documentation describing the sampler

Explicitly excluded in v1:

- multi-trait annotated BayesR
- `fast_blocks`
- `storage=:stream`
- RRM
- annotation-specific phenotype or genetic-value contribution summaries
- extra SNP-level annotation-induced class-probability output tables

## Reference Alignment

The target is to stay close to `sbayesrc.R` statistically while remaining native
to JWAS computationally.

Shared with `sbayesrc.R`:

- 4 mixture classes
- 3 conditional annotation models
- one annotation-effect variance per conditional step
- step-up indicator construction from sampled class labels
- conversion from conditional probabilities to joint class probabilities
- one-shot 4-way class sampling once the class probabilities are available

JWAS-specific differences:

- marker likelihood is individual-level and uses `yCorr`, `xRinv`, and
  `xpRinvx` rather than summary-statistics inputs
- variance updates follow the existing JWAS BayesR production path
- output and documentation follow JWAS conventions

## Model

Let the BayesR class labels be:

- class 1: zero effect
- class 2: small effect
- class 3: medium effect
- class 4: large effect

For SNP `j`, let the annotation design row including the intercept be `a_j`.
Annotated BayesR uses three conditional probit models:

- `p1_j = Pr(class > 1 | a_j) = Phi(a_j' * alpha_1)`
- `p2_j = Pr(class > 2 | class > 1, a_j) = Phi(a_j' * alpha_2)`
- `p3_j = Pr(class > 3 | class > 2, a_j) = Phi(a_j' * alpha_3)`

These imply the joint per-SNP class probabilities:

- `pi_j1 = 1 - p1_j`
- `pi_j2 = p1_j * (1 - p2_j)`
- `pi_j3 = p1_j * p2_j * (1 - p3_j)`
- `pi_j4 = p1_j * p2_j * p3_j`

The marker sampler uses `log(pi_jk)` in the BayesR class update and then samples
the class in one 4-way draw. This keeps the annotation model hierarchical while
keeping the BayesR class update simple and aligned with `sbayesrc.R`.

## Initialization

Annotated BayesR should start from the existing JWAS BayesR default:

- `pi = [0.95, 0.03, 0.015, 0.005]`

Convert that joint prior into the three conditional probabilities:

- `p1 = pi2 + pi3 + pi4 = 0.05`
- `p2 = (pi3 + pi4) / (pi2 + pi3 + pi4) = 0.4`
- `p3 = pi4 / (pi3 + pi4) = 0.25`

Then initialize the three probit intercepts as:

- `b01 = Phi^{-1}(0.05) = -1.64485362695`
- `b02 = Phi^{-1}(0.40) = -0.25334710314`
- `b03 = Phi^{-1}(0.25) = -0.67448975020`

Initialization rule:

- intercepts use the three values above
- all non-intercept annotation coefficients start at zero

This means:

- at initialization, all SNPs share the same BayesR default prior
- after annotation coefficients update, the class probabilities become SNP-specific

## Sampler

The annotated BayesR Gibbs order in JWAS should be:

1. update location parameters and refresh `yCorr`
2. sample BayesR marker classes and effects using current SNP-specific `pi_j`
3. build step-up indicators from `delta`
4. update the 3 conditional annotation probit models
5. rebuild all SNP-specific class probabilities `pi_j`
6. sample BayesR shared marker variance `sigmaSq`
7. sample residual variance

### Step-up indicators

After sampling `delta_j` for all SNPs, build:

- `z1_j = 1(delta_j > 1)`
- `z2_j = 1(delta_j > 2)`
- `z3_j = 1(delta_j > 3)`

The three conditional annotation models are then fit as:

- step 1 uses all SNPs with response `z1`
- step 2 uses only SNPs with `z1 = 1` and response `z2`
- step 3 uses only SNPs with `z2 = 1` and response `z3`

Each step:

- samples truncated-normal liabilities
- samples annotation coefficients
- samples its own annotation-effect variance
- refreshes `mu`

## Annotation-Effect Variance

Annotated BayesR uses one annotation-effect variance per conditional step, in the
same spirit as `sbayesrc.R`.

For 4 classes:

- `sigmaSqAlpha[1]` for step 1
- `sigmaSqAlpha[2]` for step 2
- `sigmaSqAlpha[3]` for step 3

The intercept is not shrunk; shrinkage applies to the non-intercept annotation
coefficients.

## Data Structure

Do not create a separate BayesR-only annotation subsystem. Generalize the
existing annotation container so Annotated BayesC remains the 1-step case and
annotated BayesR becomes the 3-step case.

For annotated BayesR, `MarkerAnnotations` should become step-aware and store:

- `design_matrix`: `m x p`
- `coefficients`: `p x 3`
- `mean_coefficients`: `p x 3`
- `mean_coefficients2`: `p x 3`
- `variance`: length `3`
- `liability`: `m x 3`
- `mu`: `m x 3`
- `lower_bound`: `m x 3`
- `upper_bound`: `m x 3`
- `lhs`: annotation-design Gram matrix
- `snp_pi`: `m x 4`

For Annotated BayesC:

- use the same structure with `nsteps = 1`

## Integration Points

Expected code changes:

- `src/1.JWAS/src/types.jl`
  - generalize `MarkerAnnotations`
- `src/1.JWAS/src/markers/readgenotypes.jl`
  - allow `annotations=` for `method="BayesR"`
  - initialize the step-aware annotation container
- `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`
  - generalize annotation initialization and update logic from 1 step to `nsteps`
  - drive annotated BayesR through the dense BayesR sweep with per-SNP class probabilities
- `src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`
  - accept SNP-specific class priors in the dense BayesR sweep
- `src/1.JWAS/src/input_data_validation.jl`
  - remove the current BayesR annotation rejection while keeping all other v1 exclusions
- `src/1.JWAS/src/output.jl`
  - output step-specific annotation coefficients for annotated BayesR
- `docs/src/manual/annotated_bayesr.md`
  - explain the model, initialization, sampler, restrictions, and output

## Output

Keep all existing BayesR outputs unchanged:

- marker effects
- `Model_Frequency = Pr(delta > 1 | data)`
- `pi`
- `sigmaSq`
- residual variance

Add one annotation-parameter output table:

- `output["annotation coefficients "*Mi.name]`

For annotated BayesR, this should be a long table with columns:

- `Annotation`
- `Step`
- `Estimate`
- `SD`

Recommended step labels:

- `step1_zero_vs_nonzero`
- `step2_small_vs_larger`
- `step3_medium_vs_large`

Deferred from v1:

- SNP-level annotation-induced class-probability output
- annotation-level phenotype or genetic-value contribution summaries

## Documentation

Add a manual page:

- `docs/src/manual/annotated_bayesr.md`

It should cover:

1. what annotated BayesR is
2. how the 3 conditional probabilities map to the 4 BayesR classes
3. default initialization from the JWAS BayesR prior
4. the Gibbs sampler order
5. current support and exclusions
6. a usage example
7. explanation of the annotation-coefficient output table

## Testing

Required validation for v1:

- API and validation tests for annotated BayesR input handling
- initialization tests for the default intercept conversion
- small dense sampler tests that confirm annotated BayesR runs end-to-end
- regression tests that Annotated BayesC still works after the annotation-container generalization
- full `julia --project=. --startup-file=no test/runtests.jl`
- docs build with `julia --project=docs --startup-file=no docs/make.jl`

## Deferred Work

Keep these out of v1:

- annotation-specific contribution to phenotype
- annotation-specific contribution to total genetic value
- multi-trait annotated BayesR
- block or streaming annotated BayesR
- extra derived annotation summaries beyond the sampled coefficients
