# Annotated BayesC

Annotated BayesC lets marker annotations change marker-state priors.

This page focuses on the original single-trait production paths:

- single-trait dense BayesC with marker-specific exclusion probabilities `π_j`
- single-trait dense BayesC with exact fast-block sweeps
- single-trait dense BayesC with optional approximate independent-block sweeps
- single-trait streaming BayesC

JWAS also supports dense 2-trait annotated BayesC with a marker-specific
4-state joint prior over `00`, `10`, `01`, and `11`. That method now has its
own guide: [Multi-Trait Annotated BayesC](multitrait_annotated_bayesc.md).

## Method Overview

Standard BayesC uses:

- `δ_j = 0` when marker `j` is excluded
- `δ_j = 1` when marker `j` is included
- a common exclusion probability `π`

Annotated BayesC replaces that common prior with annotation-driven probabilities:

```math
\Pr(\delta_j = 1 \mid a_j, \gamma) = \Phi(a_j^\top \gamma)
```

where:

- `a_j` is the annotation row for marker `j`
- `γ` is a vector of annotation coefficients
- `Φ` is the standard normal CDF

JWAS stores BayesC `π` as an exclusion probability, so the internal update is:

```math
\pi_j = 1 - \Phi(a_j^\top \gamma)
```

During MCMC, JWAS alternates between:

1. sampling marker inclusion indicators and marker effects using current `π_j`
2. sampling latent liabilities for the annotation model
3. sampling annotation coefficients `γ`
4. refreshing the per-marker exclusion probabilities `π_j`

JWAS uses the standard probit identification convention for this annotation
submodel:

- latent annotation error variance is fixed to `1`
- the annotation variance parameter controls coefficient shrinkage, not the
  latent probit noise scale

## Input Requirements

- Current support is:
  - single-trait dense `method="BayesC"`
  - single-trait dense `method="BayesC"` with `fast_blocks != false` and `independent_blocks=false`
  - single-trait dense `method="BayesC"` with `fast_blocks != false` and `independent_blocks=true`
  - single-trait streaming `method="BayesC"`
- Pass annotations through `get_genotypes(...; annotations=...)`.
- `annotations` must be a numeric matrix with one row per marker in the raw genotype input.
- JWAS applies the same marker QC/filtering mask to `annotations` as it applies to genotypes.
- JWAS prepends an intercept column automatically after filtering. Users should not include an intercept column.
- Annotation columns must be non-constant and not perfectly collinear after JWAS adds the intercept.

For the dense 2-trait method and its restrictions, see
[Multi-Trait Annotated BayesC](multitrait_annotated_bayesc.md).

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
    "genotypes.txt";
    method="BayesC",
    separator=',',
    quality_control=false,
    annotations=annotations,
)

model = build_model("y1 = intercept + genotypes")

output = runMCMC(
    model,
    phenotypes;
    chain_length=2000,
    burnin=500,
    output_samples_frequency=10,
)
```

The returned output includes an annotation-coefficient table:

```julia
output["annotation coefficients genotypes"]
```

The first row is the automatically added intercept, followed by `Annotation_1`,
`Annotation_2`, and so on.

## Multi-Trait Method

Dense 2-trait annotated BayesC is documented separately because it has a
different prior structure, startup `Pi` rules, sampler selection, and runtime
restrictions than the original single-trait method.

See [Multi-Trait Annotated BayesC](multitrait_annotated_bayesc.md) for:

- the 4-state prior over `00`, `10`, `01`, `11`
- the 3-step tree annotation model
- dense 2-trait examples
- `multi_trait_sampler=:auto|:I|:II`
- output interpretation for the joint `pi_<geno>` summary

## Dense Example with `fast_blocks`

Annotated BayesC also works with the dense block sampler:

```julia
output = runMCMC(
    model,
    phenotypes;
    chain_length=6000,
    burnin=500,
    output_samples_frequency=10,
    fast_blocks=true,
)
```

This path is still limited to the **single-trait** annotated BayesC model.
It uses the same annotation-driven `π_j` updates, but samples marker effects through the block BayesC path.
JWAS rescales the outer `chain_length` in block mode, so the nominal `chain_length` should be chosen large enough that the effective post-scaling chain is still longer than `burnin`.
You can also provide explicit block starts, for example `fast_blocks=[1, 501, 975]`.
With the default `independent_blocks=false`, this is the exact sequential block sweep.
Set `independent_blocks=true` only when you intentionally want the approximate independent-block mode for block-level thread parallelism.
For block-update details, see [Block BayesC](block_bayesc.md).

## Streaming Example

The streaming backend is also supported for **single-trait** Annotated BayesC:

```julia
using JWAS, CSV, DataFrames

prefix = prepare_streaming_genotypes(
    "genotypes.txt";
    separator=',',
    header=true,
    quality_control=false,
    center=true,
)

annotations = [
    0.0 1.0
    1.0 0.0
    1.0 1.0
    0.0 0.0
    0.5 0.5
]

genotypes = get_genotypes(
    prefix;
    method="BayesC",
    storage=:stream,
    annotations=annotations,
)

phenotypes = DataFrame(
    ID=copy(genotypes.obsID),
    y1=Float32[1.1, -0.3, 0.8, -0.9, 0.5, -0.1, 0.2],
)

model = build_model("y1 = intercept + genotypes")

output = runMCMC(
    model,
    phenotypes;
    chain_length=2000,
    burnin=500,
    output_samples_frequency=10,
    outputEBV=false,
    output_heritability=false,
    memory_guard=:off,
)
```

Streaming keeps the same single-trait annotation logic, but it still inherits the current streaming MVP constraints:

- single trait only
- BayesC only
- `fast_blocks=false`
- unit residual weights only
- `double_precision=false`
- phenotype IDs must exactly match genotype IDs in the same order
- annotated streaming requires a backend prepared by the current `prepare_streaming_genotypes` so raw-marker mapping metadata is available

For more on the streaming backend, see [Streaming Genotype Walkthrough](streaming_genotype_walkthrough.md).

## Practical Notes

- Standard BayesC is still the same model when no `annotations` are provided; JWAS simply fills the `π` vector with one repeated value.
- Annotated BayesC starts from a deterministic marker-level exclusion vector `π_j` that matches the supplied starting `Pi`.
- Annotated BayesC annotation coefficients start at zero, and JWAS does not fit the annotation model before the first phenotype-informed marker sweep.
- Annotation rows are defined on the raw marker order, not the post-QC marker order.
- If QC drops markers, JWAS drops the corresponding annotation rows before adding the intercept column.

## Equivalence With Multi-Class BayesC

For non-overlapping categorical annotations, Annotated BayesC assigns one
exclusion probability to each annotation class. This section states when that
single annotated model is equivalent to a multi-class representation that splits
the marker matrix by class.

Write the single annotated model as:

```math
y = Xb + M\alpha + e, \qquad e \sim N(0, \sigma_e^2 I).
```

For marker `j`:

```math
\alpha_j = \delta_j \beta_j,
```

where `δ_j = 0` excludes marker `j` and `δ_j = 1` includes marker `j`. The
nonzero marker effect has the BayesC prior:

```math
\beta_j \mid \sigma_\alpha^2 \sim N(0, \sigma_\alpha^2).
```

With one-hot categorical annotations, marker `j` belongs to exactly one class
`c(j) ∈ {1, ..., K}`:

```math
A_{jc} =
\begin{cases}
1, & c = c(j),\\
0, & c \ne c(j),
\end{cases}
\qquad
\sum_{c=1}^K A_{jc} = 1.
```

The annotated probit model then reduces to one class-level inclusion probability:

```math
\Pr(\delta_j = 1 \mid c(j)=c) = \Phi(\theta_c),
```

or, using the BayesC exclusion probability,

```math
\pi_j = \pi_{c(j)}, \qquad \pi_c = 1 - \Phi(\theta_c).
```

A multi-class representation splits the full marker matrix into class-specific
blocks:

```math
M = [M_1 \; M_2 \; \cdots \; M_K],
```

and writes the model as:

```math
y = Xb + M_1\alpha_1 + M_2\alpha_2 + \cdots + M_K\alpha_K + e.
```

The two likelihoods are the same when:

```math
\alpha = [\alpha_1^\top \; \alpha_2^\top \; \cdots \; \alpha_K^\top]^\top
```

and each marker appears in exactly one class block. The priors are equivalent
when both formulations use:

```math
\delta_j \mid \pi_{c(j)} \sim \mathrm{Bernoulli}(1 - \pi_{c(j)}),
```

```math
\beta_j \mid \sigma_\alpha^2 \sim N(0, \sigma_\alpha^2).
```

Thus, class membership changes the exclusion probability `π_c`, but does not
change the marker-effect variance. A multi-class model with separate
class-specific marker-effect variances,

```math
\beta_j \mid c(j)=c, \sigma_{\alpha,c}^2 \sim N(0, \sigma_{\alpha,c}^2),
```

is a different model unless the annotated model is also extended to use the same
class-specific variances.

For one-hot annotations, both formulations can have one `π_c` per class. They
are fully equivalent only if they also use the same parameterization and update
for `π_c`. Current Annotated BayesC uses the probit parameterization:

```math
\pi_c = 1 - \Phi(\theta_c).
```

A multi-class model that instead updates each class probability directly with a
BayesC-style Beta update,

```math
\pi_c \sim \mathrm{Beta}(a_c, b_c),
```

```math
\pi_c \mid \delta \sim \mathrm{Beta}(a_c + N_{0c}, b_c + N_{1c}),
```

where `N_{0c}` and `N_{1c}` are the excluded and included marker counts in class
`c`, is not the same prior model as the probit annotation model. It becomes
equivalent only if Annotated BayesC is changed to use the same direct
class-specific Beta update, or if the multi-class model is changed to use the
same probit-generated class probabilities as Annotated BayesC.
