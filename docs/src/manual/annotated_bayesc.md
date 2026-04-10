# Annotated BayesC

Annotated BayesC lets marker annotations change marker-state priors.

This page focuses on the original single-trait production paths:

- single-trait dense BayesC with marker-specific exclusion probabilities `π_j`
- single-trait dense BayesC with `fast_blocks=true`
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
  - single-trait dense `method="BayesC"` with `fast_blocks=true`
  - single-trait streaming `method="BayesC"`
- Pass annotations through `get_genotypes(...; annotations=...)`.
- `annotations` must be a numeric matrix with one row per marker in the raw genotype input.
- JWAS applies the same marker QC/filtering mask to `annotations` as it applies to genotypes.
- JWAS prepends an intercept column automatically after filtering. Users should not include an intercept column.

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
    "genotypes.txt", 1.0;
    method="BayesC",
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
    prefix, 1.0;
    method="BayesC",
    storage=:stream,
    annotations=annotations,
)

phenotypes = DataFrame(
    ID=copy(genotypes.obsID),
    y1=Float32[1.1, -0.3, 0.8, -0.9, 0.5, -0.1, 0.2],
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

For more on the streaming backend, see [Streaming Genotype Walkthrough](streaming_genotype_walkthrough.md).

## Practical Notes

- Standard BayesC is still the same model when no `annotations` are provided; JWAS simply fills the `π` vector with one repeated value.
- Annotated BayesC starts from a deterministic marker-level exclusion vector `π_j` that matches the supplied starting `Pi`.
- Annotated BayesC annotation coefficients start at zero, and JWAS does not fit the annotation model before the first phenotype-informed marker sweep.
- Annotation rows are defined on the raw marker order, not the post-QC marker order.
- If QC drops markers, JWAS drops the corresponding annotation rows before adding the intercept column.
