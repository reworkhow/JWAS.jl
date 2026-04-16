# Multi-Trait Annotated BayesC

Dense 2-trait annotated BayesC lets marker annotations change a marker's joint
inclusion prior over the four joint states:

- `00`: excluded from both traits
- `10`: included for trait 1 only
- `01`: included for trait 2 only
- `11`: included for both traits

JWAS currently supports this as a production path for **dense 2-trait BayesC**
with a marker-specific 4-state prior, under both the standard and fast-block
marker sweeps. The single-trait annotated BayesC path is documented separately
in [Annotated BayesC](annotated_bayesc.md).

## Method Overview

For each marker `j`, JWAS builds a marker-specific joint prior over the four
states above. The annotation submodel uses a 3-step binary probit tree:

- `p1_j = Pr(delta_j != 00)`
- `p2_j = Pr(delta_j = 11 | delta_j != 00)`
- `p3_j = Pr(delta_j = 10 | delta_j in {10, 01})`

The full per-marker joint prior is reconstructed as:

- `Pi_j(00) = 1 - p1_j`
- `Pi_j(11) = p1_j * p2_j`
- `Pi_j(10) = p1_j * (1 - p2_j) * p3_j`
- `Pi_j(01) = p1_j * (1 - p2_j) * (1 - p3_j)`

This keeps the annotation update Gibbs-friendly while still giving a coherent
joint prior over the four marker states.

During MCMC, JWAS alternates between:

1. sampling the multi-trait BayesC marker states and marker effects using the current marker-level joint prior
2. constructing the three binary tree responses from the current joint marker states
3. sampling latent liabilities and annotation coefficients for the three probit steps
4. rebuilding the full marker-level joint prior matrix for the next marker sweep

JWAS uses the standard probit identification convention for these annotation
submodels:

- latent annotation error variance is fixed to `1`
- the annotation variance parameters control coefficient shrinkage, not the
  latent probit noise scale

## Input Requirements

- Current support is:
  - dense 2-trait `method="BayesC"`
- Pass annotations through `get_genotypes(...; annotations=...)`.
- `annotations` must be a numeric matrix with one row per marker in the raw genotype input.
- JWAS applies the same marker QC/filtering mask to `annotations` as it applies to genotypes.
- JWAS prepends an intercept column automatically after filtering. Users should not include an intercept column.

Dense 2-trait annotated BayesC v1 requires:

- exactly 2 traits
- `storage=:dense`
- `constraint=false`
- either the standard sweep (`fast_blocks=false`) or the fast-block sweep (`fast_blocks=true`)

Dense 2-trait annotated BayesC v1 does **not** support:

- `storage=:stream`
- `constraint=true`
- more than 2 traits

If you need single-trait streaming or single-trait fast-block annotated BayesC,
those remain available on the original single-trait annotated BayesC path.

## Startup Behavior

Dense 2-trait annotated BayesC starts from a marker-level joint prior row.
JWAS accepts two startup forms:

- `Pi = 0.0`
- a joint `Pi` dictionary over the four 2-trait states

If `Pi=0.0`, JWAS keeps the legacy dense multi-trait BayesC default and starts
from:

- `Pi(00) = 0`
- `Pi(10) = 0`
- `Pi(01) = 0`
- `Pi(11) = 1`

So the startup row is `[0, 0, 0, 1]` in the fixed state order
`00, 10, 01, 11`.

If you provide a joint `Pi` dictionary, JWAS checks that:

- the probabilities sum to one
- trait 1 has positive startup mass in `{10, 11}`
- trait 2 has positive startup mass in `{01, 11}`
- the shared state `11` has positive startup mass

The dictionary is interpreted by its 2-trait state labels, not by insertion
order. Use keys such as:

- `[0.0, 0.0]`
- `[1.0, 0.0]`
- `[0.0, 1.0]`
- `[1.0, 1.0]`

Those checks happen during model setup, before marker variance scaling.

## Sampler Selection

Dense 2-trait annotated BayesC supports explicit sampler selection through
`multi_trait_sampler`:

- `:auto`
- `:I`
- `:II`

The default is `:I` for multi-trait BayesC. `:auto` still preserves the current
JWAS dispatch rule when you request it explicitly. For annotated 2-trait
BayesC, `:auto` typically uses sampler I because the model starts with all four
joint states available.

Use:

- `:I` to force the one-bit-flip sampler
- `:II` to force the full joint-state sampler

The same sampler choices apply to both:

- the standard dense sweep (`fast_blocks=false`)
- the fast-block sweep (`fast_blocks=true`)

The fast-block sweep also supports explicit block starts, for example
`fast_blocks=[1, 501, 975]`. `independent_blocks=false` is the default exact
sequential block sweep. Set `independent_blocks=true` only when you intentionally
want the approximate independent-block mode for block-level thread parallelism.
See [Block BayesC](block_bayesc.md) for the statistical assumption and server
threading guidance.

## Example

```julia
using JWAS, CSV, DataFrames

phenotypes = CSV.read("phenotypes.txt", DataFrame, delim=',', missingstring=["NA"])
phenotypes_mt = DataFrame(
    ID=copy(phenotypes.ID),
    y1=copy(phenotypes.y1),
    y2=Float32.(coalesce.(phenotypes.y1, 0.0)),
)

annotations = [
    0.0 1.0
    1.0 0.0
    1.0 1.0
    0.0 0.0
    0.5 0.5
]

Pi0 = Dict(
    [0.0, 0.0] => 0.45,
    [1.0, 0.0] => 0.20,
    [0.0, 1.0] => 0.15,
    [1.0, 1.0] => 0.20,
)

genotypes = get_genotypes(
    "genotypes.txt",
    [1.0 0.5; 0.5 1.0];
    method="BayesC",
    separator=',',
    quality_control=false,
    annotations=annotations,
    Pi=Pi0,
)

model = build_model(
    "y1 = intercept + genotypes\ny2 = intercept + genotypes",
    [1.0 0.5; 0.5 1.0],
)

output = runMCMC(
    model,
    phenotypes_mt;
    chain_length=2000,
    burnin=500,
    output_samples_frequency=10,
    fast_blocks=true,
)
```

This example uses the default `multi_trait_sampler=:I`. Add
`multi_trait_sampler=:auto` or `multi_trait_sampler=:II` in `get_genotypes(...)`
only when you want to override the default explicitly. Set
`fast_blocks=false` if you want the original non-block sweep instead.
Add `independent_blocks=true` only for the approximate independent-block block
sweep.

## Output Interpretation

The annotation output for dense 2-trait annotated BayesC uses three named steps:

- `step1_zero_vs_active`
- `step2_11_vs_singleton`
- `step3_10_vs_01`

The returned `pi_<geno>` summary is the average marker-level joint prior over the
four states `00`, `10`, `01`, and `11`.

If you inspect MCMC output files directly, the joint state labels follow the
same fixed order used throughout JWAS:

- `00`
- `10`
- `01`
- `11`

## Practical Notes

- The multi-trait annotation model is fit after the first marker sweep, so the
  annotation coefficients start at zero.
- The startup `Pi` is only the initial prior. After the first annotation update,
  JWAS learns marker-specific joint priors from the data and the annotations.
- If you want the old one-trait annotated BayesC behavior, use the single-trait
  page instead: [Annotated BayesC](annotated_bayesc.md).
