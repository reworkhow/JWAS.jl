# Multi-Trait Annotated BayesR Design

## Goal

Add dense 2-trait annotated BayesR to JWAS by extending the current dense
2-trait annotated BayesC design with BayesR effect-magnitude classes.

The method should let marker annotations affect both:

- which traits a marker affects
- the marker-effect variance class for each affected trait

The design preserves the core multi-trait BayesR idea:

- one global marker-effect covariance matrix `G` for the genotype term
- fixed BayesR variance multipliers
- marker-specific class probabilities learned from annotations
- correlated marker effects for pleiotropic markers

This design is intentionally scoped to a first production version:

- `method="BayesR"`
- two traits
- dense genotype storage
- no RRM
- no structural equation model special handling
- no diagonal-only covariance shortcut

Cross-trait marker-effect covariance is essential and is retained.

## Existing Pieces To Reuse

JWAS already has:

- single-trait BayesR with fixed class multipliers
- single-trait annotated BayesR with nested binary probit annotation updates
- dense multi-trait BayesC samplers
- dense 2-trait annotated BayesC with marker-specific 4-state priors
- inverse-Wishart sampling for multi-trait covariance matrices

The new method should reuse those pieces where possible instead of introducing
a multinomial annotation model.

## Minimum Update Principle

The implementation should be a small extension of the current dense
multi-trait annotated BayesC path. The important existing separation is:

```text
beta_j  = full latent base marker effect
alpha_j = realized marker effect used in residual updates
```

Multi-trait BayesR should keep that separation and only add BayesR class
scales:

```text
alpha_j = S_j * beta_j
```

The current single-trait `BayesR!` function should not be called literally for
each trait inside a multi-trait sampler, because it assumes scalar residual
variance and scalar marker-effect variance. That would ignore the cross-trait
terms in `R` and `G`. Instead, the implementation should reuse the single-trait
BayesR four-class logic by factoring out a small helper for "given four class
log weights, sample the class and Gaussian effect", then call that helper inside
an adapted multi-trait BayesC sampler.

If the implementation requires a large rewrite of the marker sampler, that is a
warning sign. The preferred v1 path is:

1. reuse the multi-trait BayesC Gibbs-sampler-I structure
2. replace each binary active/inactive update with a four-class BayesR update
3. keep `beta` for all markers and traits
4. compute realized `alpha` by multiplying `beta` by the BayesR scale

## Plain-Text Math Notation

The equations in this document are written in plain text so they are readable in
any Markdown viewer. The notation is:

```text
Pr(x | y)       probability of x conditional on y
Phi(x)          standard normal cumulative distribution function
N(mean, var)    normal distribution
N_2(mean, cov)  two-dimensional normal distribution
IW(df, scale)   inverse-Wishart distribution
A'              transpose of A
inv(A)          matrix inverse
sqrt(x)         square root of x
sum_j f(j)      sum over markers j
```

Subscripts are written with underscores. For example, `c_j1` means the class
for marker `j` and trait 1.

## BayesR Class Convention

For the statistical description, use class labels:

- `0`: zero effect
- `1`: small nonzero effect
- `2`: medium nonzero effect
- `3`: large nonzero effect

In JWAS implementation terms, this maps naturally to existing BayesR labels:

- class `1`: zero effect, multiplier `0.0`
- class `2`: small effect, multiplier `0.01`
- class `3`: medium effect, multiplier `0.1`
- class `4`: large effect, multiplier `1.0`

Let

```text
gamma = (0, 0.01, 0.1, 1.0)
```

in JWAS class-label order.

## State Space

For marker `j`, define a two-trait BayesR class vector:

```text
c_j = (c_j1, c_j2)
```

where each trait-specific class is one of:

```text
c_jt in {0, 1, 2, 3}
```

The complete state space has 16 states:

```text
(0,0)

(1,0), (2,0), (3,0)
(0,1), (0,2), (0,3)

(1,1), (1,2), (1,3)
(2,1), (2,2), (2,3)
(3,1), (3,2), (3,3)
```

The active-trait pattern is derived from the class vector:

```text
d_j = (1(c_j1 > 0), 1(c_j2 > 0))
d_j is one of: 00, 10, 01, 11
```

This means the method has two conceptual layers:

1. choose the active-trait pattern `00`, `10`, `01`, or `11`
2. for each active trait, choose its BayesR magnitude class

Unlike a 10-state design, this 16-state design does not force shared markers to
have the same magnitude class in both traits. A brain annotation can increase
large-effect probability for a brain-related trait without forcing the
non-brain trait to have the same magnitude.

## Marker-Effect Model

Let `G` be the global base marker-effect covariance matrix:

```text
G = [ G11  G12
      G21  G22 ]
```

For every marker, introduce an unscaled base effect:

```text
beta_j = [beta_j1, beta_j2]'
beta_j ~ N_2(mean = [0, 0]', covariance = G)
```

The realized marker effect used in the phenotype model is:

```text
alpha_j = S_j * beta_j
```

where

```text
S_j = [ sqrt(gamma[c_j1])  0
        0                  sqrt(gamma[c_j2]) ]
```

If a trait class is zero, its multiplier is zero and the realized effect for
that trait is zero.

For example, if a marker has classes `(large, small)`, then:

```text
S_j = [ sqrt(1.0)   0
        0           sqrt(0.01) ]
```

The SNP-specific realized covariance is:

```text
Var(alpha_j | c_j, G) = S_j * G * S_j
```

Therefore:

```text
Var(alpha_j1) = gamma[c_j1] * G11
```

```text
Var(alpha_j2) = gamma[c_j2] * G22
```

and

```text
Cov(alpha_j1, alpha_j2) = sqrt(gamma[c_j1] * gamma[c_j2]) * G12
```

The correlation implied by `G` is preserved for shared nonzero states, while
the covariance is scaled by the two trait-specific magnitude classes.

## Phenotype Model

For each individual, the two-trait marker contribution is:

```text
marker contribution for individual i = sum_j x_ij * alpha_j
```

Equivalently:

```text
y_i = X_i * b + sum_j x_ij * S_j * beta_j + e_i
```

with

```text
e_i ~ N_2(mean = [0, 0]', covariance = R)
```

The new BayesR machinery changes only the marker prior and marker sampler. The
fixed effects, residual covariance, and broader MCMC structure remain the same.

## Annotation Prior Factorization

For each marker, annotations determine:

```text
Pr(c_j | a_j)
```

where `a_j` is the marker annotation row with an intercept.

The proposed prior factorization is:

```text
Pr(c_j | a_j)
  = Pr(d_j | a_j)
    * Pr(c_j1 | c_j1 > 0, a_j) if trait 1 is active
    * Pr(c_j2 | c_j2 > 0, a_j) if trait 2 is active
```

This says:

- annotations first determine which traits are active
- each active trait then has its own annotation-driven BayesR magnitude
  probabilities
- trait 1 and trait 2 do not have to share magnitude probabilities

This is important for annotations that are trait-relevant, such as a brain-cell
annotation when only one trait is brain-related.

## Stage 1 Annotation Model: Active-Trait Pattern

Reuse the current multi-trait annotated BayesC tree.

Define:

```text
p1_j = Pr(d_j != 00 | a_j)
     = Phi(a_j' * theta_1)
```

```text
p2_j = Pr(d_j = 11 | d_j != 00, a_j)
     = Phi(a_j' * theta_2)
```

```text
p3_j = Pr(d_j = 10 | d_j is 10 or 01, a_j)
     = Phi(a_j' * theta_3)
```

Then:

```text
pi_j00 = 1 - p1_j
```

```text
pi_j11 = p1_j * p2_j
```

```text
pi_j10 = p1_j * (1 - p2_j) * p3_j
```

```text
pi_j01 = p1_j * (1 - p2_j) * (1 - p3_j)
```

These are the active-trait pattern probabilities.

## Stage 2 Annotation Model: Trait-Specific Magnitude

For each trait `t`, define two nested BayesR magnitude probabilities:

```text
q_t1j = Pr(c_jt >= 2 | c_jt > 0, a_j)
      = Phi(a_j' * eta_t1)
```

```text
q_t2j = Pr(c_jt = 3 | c_jt >= 2, a_j)
      = Phi(a_j' * eta_t2)
```

Then the three nonzero class probabilities for trait `t` are:

```text
m_tj1 = Pr(c_jt = 1 | c_jt > 0, a_j)
      = 1 - q_t1j
```

```text
m_tj2 = Pr(c_jt = 2 | c_jt > 0, a_j)
      = q_t1j * (1 - q_t2j)
```

```text
m_tj3 = Pr(c_jt = 3 | c_jt > 0, a_j)
      = q_t1j * q_t2j
```

Trait 1 and trait 2 use separate coefficients:

```text
eta_11, eta_12
```

for trait 1 and

```text
eta_21, eta_22
```

for trait 2.

Therefore, a marker can have high large-effect probability for trait 1 and low
large-effect probability for trait 2, even when it is active for both.

## Full 16-State Prior

Using the active-pattern probabilities and trait-specific magnitude
probabilities:

```text
Pr(c_j = (0, 0) | a_j) = pi_j00
```

For trait-1-only states:

```text
Pr(c_j = (k, 0) | a_j) = pi_j10 * m_1jk
where k is 1, 2, or 3
```

For trait-2-only states:

```text
Pr(c_j = (0, l) | a_j) = pi_j01 * m_2jl
where l is 1, 2, or 3
```

For shared states:

```text
Pr(c_j = (k, l) | a_j) = pi_j11 * m_1jk * m_2jl
where k and l are each 1, 2, or 3
```

The probabilities sum to one:

```text
pi_j00
  + pi_j10 * sum_k m_1jk
  + pi_j01 * sum_l m_2jl
  + pi_j11 * sum_k sum_l m_1jk * m_2jl
  = 1
```

## Annotation Coefficient Priors

Each binary probit annotation step uses the same prior style as current
annotated BayesC/BayesR.

For a generic annotation coefficient vector `b_h`:

- intercept has a flat or very weak prior
- annotation slopes have normal shrinkage priors

For slope `k > 1`:

```text
b_hk | sigma_h^2 ~ N(mean = 0, variance = sigma_h^2)
```

The slope variance uses the existing scaled inverse-chi-square form:

```text
sigma_h^2 = (sum over slopes k>1 of b_hk^2 + 2) / chi_square(df = p + 1)
```

where `p` is the number of annotation coefficients including the intercept.

The latent probit residual variance is fixed to one for identifiability.

## Annotation Latent Liability Updates

Each binary probit step uses a latent liability:

```text
l_jh = a_j' * b_h + epsilon_jh
epsilon_jh ~ N(mean = 0, variance = 1)
```

For binary response `z_jh`:

```text
z_jh = 1(l_jh > 0)
```

The full conditional for the liability is:

```text
l_jh | z_jh, b_h ~ N(mean = a_j' * b_h, variance = 1)
```

truncated to:

- `(0, Inf)` when `z_jh = 1`
- `(-Inf, 0]` when `z_jh = 0`

Given liabilities, each annotation submodel is a Gaussian regression and the
coefficient update remains conditionally conjugate.

## Annotation Responses

After each marker sweep, derive binary responses from the sampled class states.

### Active-Pattern Responses

Let:

```text
d_j = (1(c_j1 > 0), 1(c_j2 > 0))
```

Use the existing annotated BayesC responses:

```text
z1_j = 1(d_j != 00)
```

on all markers,

```text
z2_j = 1(d_j = 11)
```

on active markers only, and

```text
z3_j = 1(d_j = 10)
```

on singleton markers only.

### Magnitude Responses

For trait 1:

```text
h11_j = 1(c_j1 >= 2)
```

on markers with `c_j1 > 0`, and

```text
h12_j = 1(c_j1 = 3)
```

on markers with `c_j1 >= 2`.

For trait 2:

```text
h21_j = 1(c_j2 >= 2)
```

on markers with `c_j2 > 0`, and

```text
h22_j = 1(c_j2 = 3)
```

on markers with `c_j2 >= 2`.

The complete annotation model has seven binary probit steps:

1. active versus `00`
2. shared `11` versus singleton
3. `10` versus `01`
4. trait 1 medium-or-large versus small
5. trait 1 large versus medium
6. trait 2 medium-or-large versus small
7. trait 2 large versus medium

## Marker State And Effect Sampling

The marker sampler should be a direct extension of the existing dense
multi-trait BayesC Gibbs sampler I. For each marker, the sampler already loops
over traits and samples a latent `beta` value even when that trait is inactive.
BayesR should keep the same loop structure, but each trait update considers four
BayesR classes instead of two BayesC states.

Let `w_j` be the current marker right-hand-side vector, as in the dense
multi-trait BayesC sampler:

```text
w_j = [
    x_j' * (y_corr_1 + x_j * alpha_j1_old),
    x_j' * (y_corr_2 + x_j * alpha_j2_old)
]'
```

Let:

```text
q_j = x_j'x_j
```

or the weighted equivalent already used by the dense sampler.

For trait `k`, hold the other trait's current class and `beta` fixed. Evaluate
the four candidate classes for trait `k`:

```text
r in {0, 1, 2, 3}
```

with candidate scale:

```text
s_r = sqrt(gamma[r])
```

Let `l` be the other trait and let:

```text
s_l = sqrt(gamma[c_jl])
```

The scalar Gaussian conditional for candidate class `r` has precision:

```text
C_r = inv(G)[k,k] + q_j * s_r^2 * inv(R)[k,k]
```

and linear term:

```text
b_r =
    s_r * (w_j' * inv(R)[:, k])
    - ( inv(G)[k,l] + q_j * s_r * s_l * inv(R)[k,l] ) * beta_jl
```

Then:

```text
beta_jk | r, rest ~ N(mean = b_r / C_r, var = 1 / C_r)
```

The class log weight, up to constants common across the four candidate classes,
is:

```text
log_weight_r =
    log Pr(c_jk = r, c_jl = current class | a_j)
    - 0.5 * log(C_r)
    + 0.5 * b_r^2 / C_r
```

Sampling class `0` is not a special case mathematically. It has `s_r = 0`, so
the realized effect is zero, but `beta_jk` is still sampled from its conditional
latent distribution. This is the same idea as the current multi-trait BayesC
sampler, where inactive dimensions still have a sampled `beta`.

## Marker Sampler Pseudocode

```text
for marker j in 1:nMarkers:

    old_alpha = [alpha[1][j], alpha[2][j]]
    w = marker_rhs_with_old_alpha(j)

    for trait k in {1, 2}:

        l = the other trait
        current_other_class = delta[l][j]
        current_other_beta = beta[l][j]
        current_other_scale = sqrt(gamma[current_other_class])

        for candidate class r in {zero, small, medium, large}:

            candidate_scale = sqrt(gamma[r])

            C = inv(G)[k,k] + xpRinvx[j] * candidate_scale^2 * inv(R)[k,k]

            b =
                candidate_scale * dot(w, inv(R)[:, k])
                - (
                    inv(G)[k,l]
                    + xpRinvx[j] * candidate_scale * current_other_scale * inv(R)[k,l]
                  ) * current_other_beta

            mean = b / C
            var = 1 / C

            log_weight[r] =
                log_joint_annotation_prior(j, candidate class r, current_other_class)
                - 0.5 * log(C)
                + 0.5 * b^2 / C

            store mean[r], var[r]

        sampled_class = categorical_softmax(log_weight)
        delta[k][j] = sampled_class
        beta[k][j] = Normal(mean[sampled_class], var[sampled_class])
        alpha[k][j] = sqrt(gamma[sampled_class]) * beta[k][j]

        update residual array for trait k by old_alpha[k] - alpha[k][j]
```

## Sampling `G` With Full Latent `beta`

The marker sampler above produces a full two-dimensional `beta_j` for every
marker, regardless of whether the realized marker-effect pattern is `00`, `10`,
`01`, or `11`. This mirrors the current multi-trait BayesC implementation:
`beta_j` is the latent base marker effect, and `alpha_j` is the realized marker
effect after gating or scaling.

Inactive components are not realized marker effects and do not enter the
phenotype model. They are used only to keep the `G` update conditionally
conjugate and to keep cross-trait covariance on the same full latent scale as
multi-trait BayesC.

### Derivation

For every marker:

```text
beta_j | G ~ N_2(mean = [0, 0]', covariance = G)
```

Use an inverse-Wishart prior:

```text
G ~ IW(df = nu_0, scale = S_0)
```

Given full latent base effects for all markers:

```text
base effects = {beta_j for all markers j}
```

the posterior is:

```text
G | base effects
  ~ IW(
        df    = nu_0 + nMarkers,
        scale = S_0 + sum over all markers of beta_j * beta_j'
    )
```

This is the direct multi-trait analogue of the existing multi-trait BayesC
implementation, where `G` is sampled from latent `beta` rather than realized
`alpha`. It differs from the current single-trait BayesR update, which uses
only nonzero marker classes:

```text
single-trait BayesR uses:

sum over nonzero markers of alpha_j^2 / gamma[c_j]
```

That single-trait convention remains unchanged. The dense multi-trait BayesR
extension follows the multi-trait BayesC convention because cross-trait
covariance is modeled through the full latent base effect vector.

### Future Alternative: Active-Only `G` Update

The current v1 design samples and stores `beta_j` for all SNPs and uses all SNPs
in the `G` update. This is the closest match to the existing multi-trait BayesC
implementation.

A later experiment should consider using only SNPs with at least one nonzero
BayesR class in the `G` update:

```text
active for G update means c_j1 > 0 or c_j2 > 0
```

That alternative would make the multi-trait BayesR variance update closer to the
current single-trait BayesR convention and would reduce the influence of `00`
latent prior draws on `G`. The cost is that singleton SNPs would still need
their inactive trait's latent `beta` component if we want a conjugate full
cross-trait `G` update. It would also make the multi-trait BayesR covariance
interpretation less aligned with multi-trait BayesC.

This is a useful strategy to test later, especially if diagnostics show that
the all-SNP update makes `G12` too prior-driven. It is intentionally not part of
the first implementation.

## `G` Update Pseudocode

```text
SSE = prior_scale

for marker j in 1:nMarkers:

    beta_j = [beta[1][j], beta[2][j]]

    SSE += beta_j * beta_j'

G = rand(InverseWishart(G.df + nMarkers, G.scale + SSE_without_prior_if_needed))
```

Implementation should follow JWAS's existing `sample_variance` convention for
whether the prior scale is passed as `G.scale + SSE` or initialized inside the
helper. The statistical posterior is:

```text
G | base effects
  ~ IW(df = nu_0 + nMarkers,
       scale = S_0 + sum over all markers of beta_j * beta_j')
```

## Initial Marker Variance From Genetic Variance

JWAS often initializes marker-effect covariance from a target genetic covariance
and marker allele frequencies.

For multi-trait annotated BayesR, the expected genetic covariance is:

```text
V_g approx sum_j 2 * p_j * q_j * E(S_j * G * S_j | a_j)
```

Elementwise:

```text
V_g_ab approx G_ab * sum_j 2 * p_j * q_j * E(s_ja * s_jb | a_j)
```

where:

```text
s_jt = sqrt(gamma[c_jt])
```

Therefore, initialize:

```text
G_ab = V_g_ab / denominator_ab

denominator_ab = sum_j 2 * p_j * q_j * E(s_ja * s_jb | a_j)
```

For diagonal entries:

```text
E(s_jt^2 | a_j) = E(gamma[c_jt] | a_j)
```

For the off-diagonal entry:

```text
E(s_j1 * s_j2 | a_j)
  = sum over k=1..3 and l=1..3 of
      Pr(c_j = (k, l) | a_j) * sqrt(gamma[k] * gamma[l])
```

Only shared states contribute to the off-diagonal denominator because one scale
is zero for singleton states.

## Annotation Prior Refresh Pseudocode

```text
for marker j:

    p1 = Phi(A[j, :] * theta[:, 1])
    p2 = Phi(A[j, :] * theta[:, 2])
    p3 = Phi(A[j, :] * theta[:, 3])

    pattern[00] = 1 - p1
    pattern[11] = p1 * p2
    pattern[10] = p1 * (1 - p2) * p3
    pattern[01] = p1 * (1 - p2) * (1 - p3)

    for trait t in {1,2}:
        q1 = Phi(A[j, :] * eta[t, 1])
        q2 = Phi(A[j, :] * eta[t, 2])

        mag[t, small]  = 1 - q1
        mag[t, medium] = q1 * (1 - q2)
        mag[t, large]  = q1 * q2

    prior[(0,0)] = pattern[00]

    for k in small:large:
        prior[(k,0)] = pattern[10] * mag[1,k]
        prior[(0,k)] = pattern[01] * mag[2,k]

    for k in small:large:
        for l in small:large:
            prior[(k,l)] = pattern[11] * mag[1,k] * mag[2,l]

    normalize prior row defensively
```

## Full MCMC Order

For each MCMC iteration:

1. sample fixed and non-marker effects as usual
2. sample residual covariance `R` as usual
3. sample multi-trait BayesR marker states and marker effects
4. update annotation probit steps from sampled class states
5. rebuild marker-specific 16-state priors
6. sample global marker-effect covariance `G` using full latent `beta`
7. accumulate posterior means and write MCMC samples

The order of steps 4 and 6 can be aligned with the existing JWAS MCMC order as
long as the marker sampler uses the current iteration's priors and `G`.

## Data Structures

The implementation should keep these concepts separate:

- `delta[trait][marker]`: BayesR class label for that trait
- `alpha[trait][marker]`: realized scaled marker effect used in residual updates
- `beta[trait][marker]`: unscaled latent base effect used in the `G` update
- `annotations.snp_pi`: marker-specific full 16-state prior matrix
- `annotations.coefficients`: annotation coefficients for seven probit steps
- `annotations.mu`, `liability`, `lower_bound`, `upper_bound`: `nMarkers x 7`
- `G.val`: global 2-trait base marker-effect covariance matrix

The 16-state order should be explicit and centralized. A recommended order is:

```text
1  (0,0)
2  (1,0)
3  (2,0)
4  (3,0)
5  (0,1)
6  (0,2)
7  (0,3)
8  (1,1)
9  (1,2)
10 (1,3)
11 (2,1)
12 (2,2)
13 (2,3)
14 (3,1)
15 (3,2)
16 (3,3)
```

In implementation labels, add one to each class if reusing JWAS's existing
BayesR class labels `1:4`.

## Output

The method should report:

- posterior mean marker effects by trait
- posterior class probabilities by marker and trait
- posterior active-pattern probabilities `00`, `10`, `01`, `11`
- posterior shared probability `Pr(c1>0,c2>0)`
- posterior large-class probability per trait
- annotation coefficients for all seven steps
- `G` samples and posterior summary

For GWAS-style summaries, useful marker-level quantities include:

```text
Pr(c_j1 > 0 | y)
```

```text
Pr(c_j2 > 0 | y)
```

```text
Pr(c_j1 > 0 and c_j2 > 0 | y)
```

```text
Pr(c_j1 = 3 | y)
```

```text
Pr(c_j2 = 3 | y)
```

## Tests

Focused tests should cover:

- state-index mapping for all 16 states
- annotation prior rows sum to one
- trait-specific magnitude priors differ when annotation coefficients differ
- every marker, including `00`, stores a full latent `beta` vector
- `G` update uses the outer-product sum over all latent `beta` vectors
- shared markers allow unequal class pairs such as `(large, small)`
- realized `alpha` remains zero for inactive trait dimensions
- marker sampler probabilities match exact one-marker calculations
- dense 2-trait annotated BayesR runs end to end on a tiny fixture

Benchmark tests should compare:

- plain multi-trait BayesC
- dense multi-trait annotated BayesC
- single-trait annotated BayesR per trait
- dense multi-trait annotated BayesR

Primary benchmark summaries should include shared-marker recovery, trait-specific
large-effect recovery, prediction accuracy, runtime, and `G12` mixing.

## Risks

The largest statistical risk is mixing and interpretation of the cross-trait
covariance `G12`. Full latent `beta` sampling is conjugate and simple, but if
most markers are `00` or singleton states, `G12` is informed mainly through
data augmentation rather than directly observed shared realized effects.

Diagnostics should track:

- posterior number of `11` markers
- trace and ESS for `G12`
- posterior state occupancy across all 16 states
- annotation coefficient trace behavior

If `G12` mixes poorly or appears too prior-driven, a later implementation can
replace the conjugate full-latent update with an observed-block likelihood
update for `G`. That would avoid using inactive latent components in the
covariance update, but it would require Metropolis-Hastings or slice sampling
and is intentionally out of scope for v1.

## Acceptance Criteria

Dense 2-trait annotated BayesR should:

- run end to end on a small two-trait fixture
- preserve cross-trait marker-effect covariance through full `G`
- allow different magnitude classes for the two traits at shared markers
- allow trait-specific annotation effects on magnitude classes
- use all latent `beta` vectors in the `G` update, matching multi-trait BayesC
- reuse the multi-trait BayesC sampler structure and BayesR four-class logic
  instead of introducing a large independent sampler
- retain conditionally conjugate annotation coefficient updates
- reject unsupported combinations explicitly
- leave existing single-trait BayesR, annotated BayesR, and annotated BayesC
  behavior unchanged
