# Annotated BayesR

Annotated BayesR extends single-trait dense BayesR by letting marker annotations change the full four-class mixture prior for each marker.
It is the individual-level JWAS analogue of the `sbayesrc.R` summary-statistics sampler.

## Method Overview

Standard BayesR uses one shared class-probability vector:

```math
\pi = (\pi_1, \pi_2, \pi_3, \pi_4)
```

where:

- `pi_1` is the zero-effect class
- `pi_2`, `pi_3`, and `pi_4` are the nonzero mixture classes

Annotated BayesR replaces that shared prior with marker-specific class probabilities `pi_j`.
JWAS does this with three conditional probit models:

```math
p_{1j} = \Pr(\delta_j > 1 \mid a_j, \alpha_1) = \Phi(a_j^\top \alpha_1)
```

```math
p_{2j} = \Pr(\delta_j > 2 \mid \delta_j > 1, a_j, \alpha_2) = \Phi(a_j^\top \alpha_2)
```

```math
p_{3j} = \Pr(\delta_j > 3 \mid \delta_j > 2, a_j, \alpha_3) = \Phi(a_j^\top \alpha_3)
```

with:

- `a_j` as the annotation row for marker `j`
- `alpha_1`, `alpha_2`, `alpha_3` as step-specific annotation coefficients
- `Phi` as the standard normal CDF

JWAS then reconstructs the four-class prior for each marker:

```math
\pi_{j1} = 1 - p_{1j}
```

```math
\pi_{j2} = p_{1j}(1 - p_{2j})
```

```math
\pi_{j3} = p_{1j}p_{2j}(1 - p_{3j})
```

```math
\pi_{j4} = p_{1j}p_{2j}p_{3j}
```

The BayesR marker update still samples `delta_j` in one four-way draw.
The sequential structure is in the annotation model that generates `pi_j`, not in a chained marker-state sampler.

## Initialization

Annotated BayesR starts from the same default prior as standard BayesR:

```math
\pi = (0.95, 0.03, 0.015, 0.005)
```

JWAS expands this to a constant marker-level prior matrix `snp_pi`, so every
marker starts from the same four-class BayesR prior.

Annotation coefficients start at zero and JWAS does not fit the annotation model
before the first phenotype-informed marker sweep. The first BayesR sweep uses the
supplied starting `snp_pi` directly, and only then does the annotation model
update the conditional class probabilities.

## Sampler Order

For each MCMC iteration, JWAS runs annotated BayesR in this order:

1. update location parameters and `yCorr`
2. sample BayesR marker classes and marker effects using the current marker-specific `pi_j`
3. build the step-up indicators:
   - `z1_j = 1(delta_j > 1)`
   - `z2_j = 1(delta_j > 2)`
   - `z3_j = 1(delta_j > 3)`
4. update the three conditional annotation models:
   - step 1 on all markers
   - step 2 on markers with `z1_j = 1`
   - step 3 on markers with `z2_j = 1`
5. rebuild all marker-specific class probabilities `pi_j`
6. sample the shared BayesR marker variance `sigmaSq`
7. sample the residual variance

This follows the same high-level Gibbs ordering as `sbayesrc.R`, but uses JWAS individual-level marker updates instead of summary-statistics equations.

## Annotation Update In Detail

The BayesR class label for marker `j`

```math
\delta_j \in \{1,2,3,4\}
```

is still sampled in one four-way draw during the marker sweep.
The annotation model does not sample `\delta_j` sequentially.
Instead, it re-parameterizes the four class probabilities with three nested binary indicators:

```math
z_{1j} = 1(\delta_j > 1)
```

```math
z_{2j} = 1(\delta_j > 2)
```

```math
z_{3j} = 1(\delta_j > 3)
```

These indicators define three binary probit regressions:

- `z1` is fit on all markers
- `z2` is fit only on markers with `\delta_j > 1`
- `z3` is fit only on markers with `\delta_j > 2`

If `A_1`, `A_2`, and `A_3` denote those three marker sets, then:

```math
A_1 = \{1,\ldots,m\}
```

```math
A_2 = \{j : \delta_j > 1\}
```

```math
A_3 = \{j : \delta_j > 2\}
```

For each conditional model `s \in \{1,2,3\}`, JWAS uses the annotation design matrix restricted to the active marker set:

```math
X_s = A[A_s, :]
```

where `A` is the full annotation design matrix after JWAS has added the intercept column.

## Latent Probit Form

For each conditional model `s`, JWAS introduces a latent liability:

```math
l_{sj} = x_j^\top \alpha_s + \varepsilon_{sj}
```

with:

- `x_j` as the annotation row for marker `j`, including the intercept term
- `\alpha_s` as the coefficient vector for conditional model `s`
- `\varepsilon_{sj} \sim N(0,1)`

The observed binary indicator is recovered by thresholding at zero:

```math
z_{sj} = 1(l_{sj} > 0)
```

So after the current marker classes `\delta_j` are sampled, JWAS updates each conditional model by:

1. sampling latent liabilities from a truncated normal distribution
2. sampling the intercept and annotation effects by Gibbs
3. sampling the annotation-effect variance for the slope terms

## Liability Update

For active marker `j` in conditional model `s`, let

```math
\mu_{sj} = x_j^\top \alpha_s
```

Then the liability update is:

- if `z_{sj} = 0`, sample `l_{sj} \sim N(\mu_{sj}, 1)` truncated to `(-\infty, 0]`
- if `z_{sj} = 1`, sample `l_{sj} \sim N(\mu_{sj}, 1)` truncated to `[0, \infty)`

The latent error variance is fixed to `1` for identifiability.

## Coefficient Update

For one conditional model, write:

```math
l = X \alpha + \varepsilon
```

After liabilities are sampled, JWAS works with the residual

```math
r = l - X\alpha
```

and updates the coefficient vector `\alpha` one element at a time.

The first coefficient is the intercept. Its prior is flat, so if `n = |A_s|` is the number of active markers in the current conditional model, then:

```math
\alpha_0 \mid \cdot \sim N(\hat{\alpha}_0, 1/n)
```

with:

```math
\hat{\alpha}_0 = \frac{\sum_i r_i + n \alpha_0^{old}}{n}
```

The remaining coefficients are annotation effects. Each annotation effect uses a normal prior:

```math
\alpha_k \sim N(0, \sigma^2_{\alpha,s})
```

For annotation column `x_k`, the Gibbs update is:

```math
\alpha_k \mid \cdot \sim N(\hat{\alpha}_k, v_k)
```

where:

```math
v_k = \frac{1}{x_k^\top x_k + 1/\sigma^2_{\alpha,s}}
```

```math
\hat{\alpha}_k = v_k \left(x_k^\top r + (x_k^\top x_k)\alpha_k^{old}\right)
```

So the intercept is updated from the active marker count, while the annotation effects are updated from the active-subset annotation crossproducts plus the step-specific shrinkage variance.

## Annotation-Effect Variance Update

Each conditional model has its own slope-variance parameter `\sigma^2_{\alpha,s}`.
After the slope coefficients are updated, JWAS samples that variance with the same scaled inverse-chi-square form used in the reference implementation:

```math
\sigma^2_{\alpha,s} = \frac{\sum_{k>1}\alpha_{k,s}^2 + 2}{\chi^2_{p+1}}
```

where `p` is the total number of coefficients in that conditional model, including the intercept.

The intercept is not shrunk by this variance update.

## Prior Reconstruction

After all three conditional models are updated, JWAS evaluates:

```math
p_{1j} = \Phi(x_j^\top \alpha_1)
```

```math
p_{2j} = \Phi(x_j^\top \alpha_2)
```

```math
p_{3j} = \Phi(x_j^\top \alpha_3)
```

for every marker `j`, and then reconstructs the per-marker four-class prior matrix:

```math
\pi_{j1} = 1 - p_{1j}
```

```math
\pi_{j2} = p_{1j}(1 - p_{2j})
```

```math
\pi_{j3} = p_{1j}p_{2j}(1 - p_{3j})
```

```math
\pi_{j4} = p_{1j}p_{2j}p_{3j}
```

These `\pi_{jk}` values are the class priors used in the next BayesR marker sweep.

## Pseudocode

The annotation part of one MCMC iteration is:

```text
given current marker classes delta_j:
    z1_j = 1(delta_j > 1)
    z2_j = 1(delta_j > 2)
    z3_j = 1(delta_j > 3)

    A1 = all markers
    A2 = markers with delta_j > 1
    A3 = markers with delta_j > 2

    for s in {1, 2, 3}:
        X_s = annotation design on active set A_s
        z_s = binary response on active set A_s

        sample latent liabilities l_s from truncated normals
        sample intercept with flat prior using n_s = number of active markers
        sample annotation slopes with N(0, sigmaSqAlpha_s) priors
        sample sigmaSqAlpha_s from the updated slopes

    for every marker j:
        p1_j = Phi(x_j' alpha_1)
        p2_j = Phi(x_j' alpha_2)
        p3_j = Phi(x_j' alpha_3)

        pi_j1 = 1 - p1_j
        pi_j2 = p1_j * (1 - p2_j)
        pi_j3 = p1_j * p2_j * (1 - p3_j)
        pi_j4 = p1_j * p2_j * p3_j
```

## Input Requirements

- Current support is **single-trait `method="BayesR"` only**.
- Current support is **dense storage only**.
- Pass annotations through `get_genotypes(...; annotations=...)`.
- `annotations` must be a numeric matrix with one row per marker in the raw genotype input.
- JWAS applies the same marker QC/filtering mask to `annotations` as it applies to genotypes.
- JWAS prepends an intercept column automatically after filtering. Users should not include an intercept column.

Current v1 exclusions:

- `storage=:stream`
- multi-trait BayesR
- random regression models (`RRM`)

`fast_blocks` is supported for dense annotated BayesR. The block sampler uses the
same annotation-induced marker-specific class probabilities `pi_j` as the dense
sampler. As with ordinary BayesR, block mode is an accelerated approximation to
the dense transition kernel rather than the exact same sampler.
Explicit block starts such as `fast_blocks=[1, 501, 975]` are also supported.
Set `independent_blocks=true` only when you intentionally want the approximate
independent-block mode for block-level thread parallelism. See
[Block BayesC](block_bayesc.md) for the shared block-sampler interpretation.

## Dense Example

```julia
using JWAS, CSV, DataFrames

phenotypes = CSV.read("phenotypes.txt", DataFrame, delim=',', missingstring=["NA"])

annotations = [
    0.0 1.0
    1.0 0.0
    1.0 1.0
    0.0 0.0
    0.5 0.5
]

genotypes = get_genotypes(
    "genotypes.txt", 1.0;
    method="BayesR",
    separator=',',
    quality_control=false,
    annotations=annotations,
)

model = build_model("y1 = intercept + genotypes", 1.0)

output = runMCMC(
    model,
    phenotypes;
    chain_length=2000,
    burnin=500,
    output_samples_frequency=10,
    outputEBV=false,
    output_heritability=false,
)
```

## Dense Block Example

```julia
output = runMCMC(
    model,
    phenotypes;
    chain_length=2000,
    burnin=500,
    output_samples_frequency=10,
    fast_blocks=true,
    outputEBV=false,
    output_heritability=false,
)
```

## Output

Annotated BayesR keeps the standard BayesR outputs, including:

- marker effects
- posterior `Model_Frequency`
- posterior shared marker variance
- posterior residual variance

It also keeps the `pi_<genotype name>` table, but the meaning is slightly different from ordinary BayesR.
For annotated BayesR, that table reports posterior means of the current annotation-induced class probabilities averaged across markers.

It also adds a step-specific annotation-coefficient table:

```julia
output["annotation coefficients genotypes"]
```

with columns:

- `Annotation`
- `Step`
- `Estimate`
- `SD`

The `Step` labels are:

- `step1_zero_vs_nonzero`
- `step2_small_vs_larger`
- `step3_medium_vs_large`

These are the sampled annotation-model parameters.
They describe how annotations change the current conditional prior class probabilities, not the final posterior class probabilities.

## Practical Notes

- Standard BayesR is unchanged when no `annotations` are provided.
- Annotation rows are defined on the raw marker order, not the post-QC marker order.
- If QC drops markers, JWAS drops the corresponding annotation rows before adding the intercept column.
- Posterior PIP is still read from the marker-effects output as `Model_Frequency = Pr(delta_j > 1 | data)`.
