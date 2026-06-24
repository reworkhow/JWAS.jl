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

```math
\gamma = (0,\;0.01,\;0.1,\;1.0)
```

in JWAS class-label order.

## State Space

For marker `j`, define a two-trait BayesR class vector:

```math
c_j = (c_{j1}, c_{j2}),
```

where each trait-specific class is one of:

```math
c_{jt} \in \{0,1,2,3\}.
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

```math
d_j =
\left(
1(c_{j1}>0),
1(c_{j2}>0)
\right)
\in \{00,10,01,11\}.
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

```math
G =
\begin{bmatrix}
G_{11} & G_{12}\\
G_{21} & G_{22}
\end{bmatrix}.
```

For every active marker, introduce an unscaled base effect:

```math
u_j =
\begin{bmatrix}
u_{j1}\\
u_{j2}
\end{bmatrix}
\sim N_2(0,G).
```

The realized marker effect used in the phenotype model is:

```math
\alpha_j = S_j u_j,
```

where

```math
S_j =
\begin{bmatrix}
\sqrt{\gamma_{c_{j1}}} & 0\\
0 & \sqrt{\gamma_{c_{j2}}}
\end{bmatrix}.
```

If a trait class is zero, its multiplier is zero and the realized effect for
that trait is zero.

For example, if a marker has classes `(large, small)`, then:

```math
S_j =
\begin{bmatrix}
\sqrt{1.0} & 0\\
0 & \sqrt{0.01}
\end{bmatrix}.
```

The SNP-specific realized covariance is:

```math
\operatorname{Var}(\alpha_j \mid c_j,G)
=
S_j G S_j.
```

Therefore:

```math
\operatorname{Var}(\alpha_{j1}) = \gamma_{c_{j1}}G_{11},
```

```math
\operatorname{Var}(\alpha_{j2}) = \gamma_{c_{j2}}G_{22},
```

and

```math
\operatorname{Cov}(\alpha_{j1},\alpha_{j2})
=
\sqrt{\gamma_{c_{j1}}\gamma_{c_{j2}}}G_{12}.
```

The correlation implied by `G` is preserved for shared nonzero states, while
the covariance is scaled by the two trait-specific magnitude classes.

## Phenotype Model

For each individual, the two-trait marker contribution is:

```math
\sum_j x_{ij}\alpha_j.
```

Equivalently:

```math
y_i = X_i b + \sum_{j=1}^m x_{ij}S_j u_j + e_i,
```

with

```math
e_i \sim N_2(0,R).
```

The new BayesR machinery changes only the marker prior and marker sampler. The
fixed effects, residual covariance, and broader MCMC structure remain the same.

## Annotation Prior Factorization

For each marker, annotations determine:

```math
Pr(c_j \mid a_j),
```

where `a_j` is the marker annotation row with an intercept.

The proposed prior factorization is:

```math
Pr(c_j \mid a_j)
=
Pr(d_j \mid a_j)
Pr(c_{j1}\mid c_{j1}>0,a_j)^{1(d_{j1}=1)}
Pr(c_{j2}\mid c_{j2}>0,a_j)^{1(d_{j2}=1)}.
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

```math
p_{1j} = Pr(d_j \ne 00 \mid a_j)
       = \Phi(a_j'\theta_1),
```

```math
p_{2j} = Pr(d_j = 11 \mid d_j \ne 00, a_j)
       = \Phi(a_j'\theta_2),
```

```math
p_{3j} = Pr(d_j = 10 \mid d_j \in \{10,01\}, a_j)
       = \Phi(a_j'\theta_3).
```

Then:

```math
\pi_{j,00} = 1 - p_{1j},
```

```math
\pi_{j,11} = p_{1j}p_{2j},
```

```math
\pi_{j,10} = p_{1j}(1-p_{2j})p_{3j},
```

```math
\pi_{j,01} = p_{1j}(1-p_{2j})(1-p_{3j}).
```

These are the active-trait pattern probabilities.

## Stage 2 Annotation Model: Trait-Specific Magnitude

For each trait `t`, define two nested BayesR magnitude probabilities:

```math
q_{t1,j} =
Pr(c_{jt} \ge 2 \mid c_{jt}>0,a_j)
=
\Phi(a_j'\eta_{t1}),
```

```math
q_{t2,j} =
Pr(c_{jt} = 3 \mid c_{jt}\ge 2,a_j)
=
\Phi(a_j'\eta_{t2}).
```

Then the three nonzero class probabilities for trait `t` are:

```math
m_{t,j,1} = Pr(c_{jt}=1 \mid c_{jt}>0,a_j)
          = 1 - q_{t1,j},
```

```math
m_{t,j,2} = Pr(c_{jt}=2 \mid c_{jt}>0,a_j)
          = q_{t1,j}(1-q_{t2,j}),
```

```math
m_{t,j,3} = Pr(c_{jt}=3 \mid c_{jt}>0,a_j)
          = q_{t1,j}q_{t2,j}.
```

Trait 1 and trait 2 use separate coefficients:

```math
\eta_{11}, \eta_{12}
```

for trait 1 and

```math
\eta_{21}, \eta_{22}
```

for trait 2.

Therefore, a marker can have high large-effect probability for trait 1 and low
large-effect probability for trait 2, even when it is active for both.

## Full 16-State Prior

Using the active-pattern probabilities and trait-specific magnitude
probabilities:

```math
Pr(c_j=(0,0)\mid a_j) = \pi_{j,00}.
```

For trait-1-only states:

```math
Pr(c_j=(k,0)\mid a_j)
=
\pi_{j,10}m_{1,j,k},
\quad k \in \{1,2,3\}.
```

For trait-2-only states:

```math
Pr(c_j=(0,l)\mid a_j)
=
\pi_{j,01}m_{2,j,l},
\quad l \in \{1,2,3\}.
```

For shared states:

```math
Pr(c_j=(k,l)\mid a_j)
=
\pi_{j,11}m_{1,j,k}m_{2,j,l},
\quad k,l \in \{1,2,3\}.
```

The probabilities sum to one:

```math
\pi_{j,00}
+\pi_{j,10}\sum_k m_{1,j,k}
+\pi_{j,01}\sum_l m_{2,j,l}
+\pi_{j,11}\sum_k\sum_l m_{1,j,k}m_{2,j,l}
=1.
```

## Annotation Coefficient Priors

Each binary probit annotation step uses the same prior style as current
annotated BayesC/BayesR.

For a generic annotation coefficient vector `b_h`:

- intercept has a flat or very weak prior
- annotation slopes have normal shrinkage priors

For slope `k > 1`:

```math
b_{hk}\mid \sigma_h^2 \sim N(0,\sigma_h^2).
```

The slope variance uses the existing scaled inverse-chi-square form:

```math
\sigma_h^2 =
\frac{\sum_{k>1}b_{hk}^2 + 2}{\chi^2_{p+1}},
```

where `p` is the number of annotation coefficients including the intercept.

The latent probit residual variance is fixed to one for identifiability.

## Annotation Latent Liability Updates

Each binary probit step uses a latent liability:

```math
\ell_{jh} = a_j'b_h + \epsilon_{jh},
\quad
\epsilon_{jh}\sim N(0,1).
```

For binary response `z_jh`:

```math
z_{jh}=1(\ell_{jh}>0).
```

The full conditional for the liability is:

```math
\ell_{jh}\mid z_{jh},b_h
\sim
N(a_j'b_h,1)
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

```math
d_j=(1(c_{j1}>0),1(c_{j2}>0)).
```

Use the existing annotated BayesC responses:

```math
z_{1j}=1(d_j\ne 00)
```

on all markers,

```math
z_{2j}=1(d_j=11)
```

on active markers only, and

```math
z_{3j}=1(d_j=10)
```

on singleton markers only.

### Magnitude Responses

For trait 1:

```math
h_{11,j}=1(c_{j1}\ge 2)
```

on markers with `c_j1 > 0`, and

```math
h_{12,j}=1(c_{j1}=3)
```

on markers with `c_j1 >= 2`.

For trait 2:

```math
h_{21,j}=1(c_{j2}\ge 2)
```

on markers with `c_j2 > 0`, and

```math
h_{22,j}=1(c_{j2}=3)
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

The marker sampler should use a two-stage implementation while preserving the
correct BayesR likelihood contribution.

For marker `j`, define all 16 possible class states `c`.

For each state `c`, build:

```math
S_c =
\begin{bmatrix}
\sqrt{\gamma_{c_1}} & 0\\
0 & \sqrt{\gamma_{c_2}}
\end{bmatrix}.
```

Let `w_j` be the current marker right-hand-side vector, analogous to the dense
multi-trait BayesC sampler's `w`:

```math
w_j =
\begin{bmatrix}
x_j'(y_{\text{corr},1}+x_j\alpha_{j1}^{old})\\
x_j'(y_{\text{corr},2}+x_j\alpha_{j2}^{old})
\end{bmatrix}.
```

Let:

```math
q_j = x_j'x_j
```

or the weighted equivalent already used by the dense sampler.

Conditional on state `c`, the unscaled base effect has Gaussian full
conditional:

```math
u_j\mid c,\text{rest}
\sim
N_2(\mu_{jc}, C_{jc}^{-1}),
```

where:

```math
C_{jc}
=
G^{-1} + q_j S_c R^{-1}S_c,
```

and:

```math
b_{jc}
=
S_c R^{-1}w_j,
```

```math
\mu_{jc}=C_{jc}^{-1}b_{jc}.
```

The state log weight, up to constants common across states, is:

```math
\log W_{jc}
=
\log Pr(c_j=c\mid a_j)
-\frac{1}{2}\log|C_{jc}|
+\frac{1}{2}b_{jc}'C_{jc}^{-1}b_{jc}.
```

This is the same Gaussian integration idea used in multi-trait BayesC, but with
state-specific scaling matrix `S_c`.

### Two-Stage Sampling

Do not ignore magnitude classes when sampling the active-trait pattern. Instead,
aggregate the full 16-state weights.

For each active pattern `d`, compute:

```math
W_{jd} = \sum_{c:\;d(c)=d} W_{jc}.
```

Then sample:

```math
d_j \sim \text{Categorical}(W_{j,00},W_{j,10},W_{j,01},W_{j,11}).
```

After `d_j` is chosen, sample the class state within that pattern:

```math
c_j \mid d_j
\sim
\text{Categorical}\{W_{jc}:d(c)=d_j\}.
```

Finally sample:

```math
u_j \sim N_2(\mu_{jc_j},C_{jc_j}^{-1})
```

and set:

```math
\alpha_j = S_{c_j}u_j.
```

The realized effect `alpha_j` updates the phenotype residual. The unscaled
effect `u_j` is stored for the `G` update.

For `00`, `S_c=0`, `alpha_j=0`, and the sampled `u_j` is not needed for the
`G` update.

## Marker Sampler Pseudocode

```text
for marker j in 1:nMarkers:

    old_alpha = alpha[j, :]
    w = marker_rhs_with_old_alpha(j)

    for each state c in 16 states:

        prior[c] = annotation_prior[j, c]

        S = diag(sqrt(gamma[c1]), sqrt(gamma[c2]))

        C = inv(G) + xpRinvx[j] * S * inv(R) * S
        b = S * inv(R) * w
        mu = inv(C) * b

        log_weight[c] =
            log(prior[c])
            - 0.5 * logdet(C)
            + 0.5 * b' * inv(C) * b

        store mu[c], inv(C)

    for d in {00, 10, 01, 11}:
        log_pattern_weight[d] =
            logsumexp(log_weight[c] for states c with pattern d)

    sampled_pattern = categorical_softmax(log_pattern_weight)

    sampled_state =
        categorical_softmax(log_weight[c] for states c with sampled_pattern)

    u[j, :] = Normal(mu[sampled_state], invC[sampled_state])

    S = scale_matrix(sampled_state)
    alpha[j, :] = S * u[j, :]

    update residual arrays by old_alpha - alpha[j, :]

    delta[1][j] = sampled class for trait 1
    delta[2][j] = sampled class for trait 2
```

## Sampling `G` With Latent `u` Completion

The marker sampler above already produces a full two-dimensional `u_j` for each
nonzero marker. For singleton markers, the inactive trait's component of `u_j`
is a latent completion sampled from the conditional Gaussian posterior implied
by the current `G` and the observed active-trait effect.

The inactive component is not a realized marker effect and does not enter the
phenotype model. It is only used to keep the `G` update conjugate.

### Derivation

For active markers:

```math
u_j\mid G \sim N_2(0,G).
```

Use an inverse-Wishart prior:

```math
G \sim IW(\nu_0,S_0).
```

Given completed base effects for currently active markers:

```math
\{u_j:c_j\ne(0,0)\},
```

the posterior is:

```math
G\mid \{u_j\}
\sim
IW\left(
\nu_0+n_{\text{active}},
S_0+\sum_{j:c_j\ne(0,0)}u_ju_j'
\right).
```

This is the direct multi-trait analogue of single-trait BayesR, where the
global variance update uses:

```math
\sum_{j:c_j>0}\alpha_j^2/\gamma_{c_j}.
```

Here:

```math
u_j = S_j^{-1}\alpha_j
```

for active dimensions, with missing dimensions completed by the Gaussian
posterior.

### Explicit Singleton Completion

If the marker sampler stores only active dimensions, the missing dimension can
be sampled explicitly.

For a `10` marker:

```math
u_{j1} = \alpha_{j1}/\sqrt{\gamma_{c_{j1}}}.
```

Then:

```math
u_{j2}\mid u_{j1},G
\sim
N\left(
\frac{G_{21}}{G_{11}}u_{j1},
G_{22}-\frac{G_{21}^2}{G_{11}}
\right).
```

For a `01` marker:

```math
u_{j2} = \alpha_{j2}/\sqrt{\gamma_{c_{j2}}},
```

and:

```math
u_{j1}\mid u_{j2},G
\sim
N\left(
\frac{G_{12}}{G_{22}}u_{j2},
G_{11}-\frac{G_{12}^2}{G_{22}}
\right).
```

For a `11` marker:

```math
u_j =
\begin{bmatrix}
\alpha_{j1}/\sqrt{\gamma_{c_{j1}}}\\
\alpha_{j2}/\sqrt{\gamma_{c_{j2}}}
\end{bmatrix}.
```

For a `00` marker, skip the marker in the `G` update.

## `G` Update Pseudocode

```text
SSE = prior_scale
n_active = 0

for marker j in 1:nMarkers:

    c1 = delta[1][j]
    c2 = delta[2][j]

    if c1 == zero_class and c2 == zero_class:
        continue

    n_active += 1

    if u[j, :] was stored by marker sampler:
        uj = u[j, :]

    else:
        if c1 nonzero and c2 nonzero:
            u1 = alpha[1][j] / sqrt(gamma[c1])
            u2 = alpha[2][j] / sqrt(gamma[c2])

        else if c1 nonzero and c2 zero:
            u1 = alpha[1][j] / sqrt(gamma[c1])
            u2 = draw_conditional_normal(u2 | u1, current_G)

        else if c1 zero and c2 nonzero:
            u2 = alpha[2][j] / sqrt(gamma[c2])
            u1 = draw_conditional_normal(u1 | u2, current_G)

        uj = [u1, u2]

    SSE += uj * uj'

G = rand(InverseWishart(G.df + n_active, G.scale + SSE_without_prior_if_needed))
```

Implementation should follow JWAS's existing `sample_variance` convention for
whether the prior scale is passed as `G.scale + SSE` or initialized inside the
helper. The statistical posterior is:

```math
IW(\nu_0+n_{\text{active}}, S_0+\sum u_ju_j').
```

## Initial Marker Variance From Genetic Variance

JWAS often initializes marker-effect covariance from a target genetic covariance
and marker allele frequencies.

For multi-trait annotated BayesR, the expected genetic covariance is:

```math
V_g
\approx
\sum_j 2p_jq_j \; E(S_jGS_j\mid a_j).
```

Elementwise:

```math
V_{g,ab}
\approx
G_{ab}
\sum_j 2p_jq_j \; E(s_{ja}s_{jb}\mid a_j),
```

where:

```math
s_{jt} = \sqrt{\gamma_{c_{jt}}}.
```

Therefore, initialize:

```math
G_{ab}
=
\frac{V_{g,ab}}
{\sum_j 2p_jq_j E(s_{ja}s_{jb}\mid a_j)}.
```

For diagonal entries:

```math
E(s_{jt}^2\mid a_j)
=
E(\gamma_{c_{jt}}\mid a_j).
```

For the off-diagonal entry:

```math
E(s_{j1}s_{j2}\mid a_j)
=
\sum_{k=1}^3\sum_{l=1}^3
Pr(c_j=(k,l)\mid a_j)
\sqrt{\gamma_k\gamma_l}.
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
6. sample global marker-effect covariance `G` using completed base effects `u`
7. accumulate posterior means and write MCMC samples

The order of steps 4 and 6 can be aligned with the existing JWAS MCMC order as
long as the marker sampler uses the current iteration's priors and `G`.

## Data Structures

The implementation should keep these concepts separate:

- `delta[trait][marker]`: BayesR class label for that trait
- `alpha[trait][marker]`: realized scaled marker effect used in residual updates
- `beta[trait][marker]` or a new equivalent: unscaled base effect `u`
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

```math
Pr(c_{j1}>0\mid y),
```

```math
Pr(c_{j2}>0\mid y),
```

```math
Pr(c_{j1}>0,c_{j2}>0\mid y),
```

```math
Pr(c_{j1}=3\mid y),
```

```math
Pr(c_{j2}=3\mid y).
```

## Tests

Focused tests should cover:

- state-index mapping for all 16 states
- annotation prior rows sum to one
- trait-specific magnitude priors differ when annotation coefficients differ
- `00` markers do not contribute to `G`
- singleton markers produce completed `u` vectors for the `G` update
- shared markers allow unequal class pairs such as `(large, small)`
- `G` update reduces to single-trait BayesR sufficient statistics on diagonal
  when only one trait is active
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

The largest statistical risk is mixing of the cross-trait covariance `G12`.
Latent completion is conjugate and simple, but if most active markers are
singleton states, `G12` is weakly informed by observed shared marker effects.

Diagnostics should track:

- posterior number of `11` markers
- trace and ESS for `G12`
- posterior state occupancy across all 16 states
- annotation coefficient trace behavior

If `G12` mixes poorly, a later implementation can replace the conjugate latent
completion update with an observed-block likelihood update for `G`. That would
avoid latent completions, but it would require Metropolis-Hastings or slice
sampling and is intentionally out of scope for v1.

## Acceptance Criteria

Dense 2-trait annotated BayesR should:

- run end to end on a small two-trait fixture
- preserve cross-trait marker-effect covariance through full `G`
- allow different magnitude classes for the two traits at shared markers
- allow trait-specific annotation effects on magnitude classes
- use only nonzero markers in the `G` update
- retain conditionally conjugate annotation coefficient updates
- reject unsupported combinations explicitly
- leave existing single-trait BayesR, annotated BayesR, and annotated BayesC
  behavior unchanged
