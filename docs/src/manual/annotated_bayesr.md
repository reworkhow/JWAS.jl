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
