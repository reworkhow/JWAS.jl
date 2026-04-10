# Multi-Trait Annotated BayesC Design

## Goal

Add dense 2-trait annotated BayesC to the production JWAS MCMC path.

The method should:

- keep the existing dense multi-trait BayesC likelihood and marker-effect updates
- reuse the current annotation Gibbs machinery pattern based on binary probit updates
- avoid introducing a multinomial or logistic annotation sampler
- preserve existing single-trait annotated BayesC behavior

This design is intentionally limited to:

- `method="BayesC"`
- `storage=:dense`
- `ntraits == 2`
- `constraint = false`
- `fast_blocks = false`

Unsupported paths should continue to fail explicitly.

## Existing Constraints

JWAS already has:

- dense multi-trait BayesC samplers in `MTBayesABC.jl`
- single-trait annotated BayesC with a 1-step binary probit annotation update
- single-trait annotated BayesR with a 3-step sequential annotation parameterization that rebuilds a joint per-marker class prior

The missing piece is a dense 2-trait annotated BayesC prior that:

- produces a coherent joint prior on the four 2-trait inclusion states
- stays compatible with Gibbs updates
- can be consumed by both current multi-trait BayesC sampler styles

## Rejected Designs

### Direct 4-state multinomial prior

Learning `Pr(00)`, `Pr(10)`, `Pr(01)`, and `Pr(11)` directly from annotations is statistically reasonable, but it does not fit the current simple Gibbs/probit annotation framework. It would require a multinomial logit/probit-style sampler or additional augmentation machinery.

### Hurdle-only design

The hurdle design only lets annotations act on `00` versus `active`, with a shared annotation-free split among `10`, `01`, and `11`. This is too restrictive for the intended model.

### Four independent edge-conditionals

Modeling the four one-bit edge probabilities separately over-parameterizes the 4-state joint prior. A valid 4-state prior has three free degrees of freedom, so four independent binary models do not automatically reconstruct one coherent joint distribution.

## Chosen Design: 3-Step Tree Joint Prior

For each marker `j`, define three annotation-driven binary probabilities:

- `p1_j = Pr(delta_j != 00)`
- `p2_j = Pr(delta_j = 11 | delta_j != 00)`
- `p3_j = Pr(delta_j = 10 | delta_j in {10, 01})`

These are fit with three binary probit models:

- `p1_j = Phi(a_j' gamma_1)`
- `p2_j = Phi(a_j' gamma_2)`
- `p3_j = Phi(a_j' gamma_3)`

where `a_j` is the annotation row for marker `j`.

The full joint per-marker prior is reconstructed as:

- `Pi_j(00) = 1 - p1_j`
- `Pi_j(11) = p1_j * p2_j`
- `Pi_j(10) = p1_j * (1 - p2_j) * p3_j`
- `Pi_j(01) = p1_j * (1 - p2_j) * (1 - p3_j)`

This gives a valid 4-state prior for every marker and keeps the annotation updates in the current Gibbs-friendly binary probit form.

## Why This Tree Is Acceptable

This parameterization is not an ordered-scale model in the same sense as single-trait annotated BayesR. However, it is still a meaningful binary tree:

1. null versus active
2. pleiotropic versus singleton
3. singleton trait 1 versus singleton trait 2

The first two splits are naturally interpretable. The final `10` versus `01` split is not an ordered effect-size split, but it is still a coherent binary partition of the singleton states. That is sufficient for a valid Gibbs-friendly joint prior.

## Interaction With Multi-Trait BayesC Samplers

For two traits, the joint states are:

- `00`
- `10`
- `01`
- `11`

The dense multi-trait BayesC code already supports two sampler styles.

### Sampler I: one indicator changes at a time

Sampler I only compares neighboring states. For marker `j`, the relevant prior odds become:

- `00` versus `10`: `Pi_j(10) / Pi_j(00)`
- `01` versus `11`: `Pi_j(11) / Pi_j(01)`
- `00` versus `01`: `Pi_j(01) / Pi_j(00)`
- `10` versus `11`: `Pi_j(11) / Pi_j(10)`

So sampler I does not need a direct `00` versus `11` jump. It only needs marker-specific joint prior values for the compared neighboring states.

### Sampler II: full-state draw

Sampler II can consume the same `Pi_j(00:11)` values directly in a 4-state draw.

### Consistency

The tree prior is path-consistent under sampler I:

- `(Pi_j(10) / Pi_j(00)) * (Pi_j(11) / Pi_j(10)) = Pi_j(11) / Pi_j(00)`
- `(Pi_j(01) / Pi_j(00)) * (Pi_j(11) / Pi_j(01)) = Pi_j(11) / Pi_j(00)`

So the prior does not depend on whether the path to `11` goes through `10` or `01`.

## Annotation Update

After each multi-trait BayesC marker sweep, derive three binary responses:

- `z1_j = 1(delta_j != 00)` on all markers
- `z2_j = 1(delta_j == 11)` on active markers `{10, 01, 11}`
- `z3_j = 1(delta_j == 10)` on singleton markers `{10, 01}`

Then run the usual probit Gibbs update separately for each step:

- sample liabilities
- sample annotation coefficients
- optionally sample step-specific annotation variances using the existing framework
- rebuild `Pi_j(00:11)` for all markers

The first annotation update should happen only after the first phenotype-informed marker sweep.

## Startup

Startup should preserve the current annotated BayesC principle:

- use the supplied starting prior first
- start annotation coefficients at zero
- do not pre-fit the annotation model before the first marker sweep

For 2-trait BayesC, if the caller supplies the usual joint starting prior:

- `Pi(00)`
- `Pi(10)`
- `Pi(01)`
- `Pi(11)`

then the initial per-marker prior matrix should simply repeat that same row for every marker. Annotation coefficients and annotation linear predictors start at zero and are fit only after the first sweep.

## Data Representation

The current single-trait annotated BayesC path specializes too early to a scalar/vector `Pi`. For 2-trait annotated BayesC, the annotation container should carry the current per-marker 4-state prior matrix directly, analogous to annotated BayesR:

- `nsteps = 3`
- `nclasses = 4`
- `snp_pi[m, 4]`

Single-trait annotated BayesC remains the 1-step case.

## Output

The marker sampler should consume per-marker joint priors. User-facing output should still remain interpretable:

- annotation coefficients should be written step-by-step
- the reported `pi_<geno>` summary for the chain should be the average marker-level joint prior over the four states

## Acceptance Criteria

Dense 2-trait `method="BayesC"` with `annotations=...` should:

- initialize successfully on the production path
- run end-to-end through the dense MCMC loop
- use marker-specific annotation priors in the multi-trait marker sampler
- update the three annotation steps after each marker sweep
- preserve current single-trait annotated BayesC behavior
- reject unsupported combinations explicitly
