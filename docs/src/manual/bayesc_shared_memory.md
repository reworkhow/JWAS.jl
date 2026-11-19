# BayesC shared genotype memory

Dense BayesC supports shared genotype storage for single- and multi-trait models,
with ordinary marker updates or fast blocks. Both existing multi-trait samplers
I and II and sequential/independent blocks are supported.

For matching-precision matrix input, BayesC defaults to copy_genotypes=false.
Imputation and centering can modify the input. Use copy_genotypes=true to retain
the original input or permit precision conversion. With quality_control=true,
JWAS keeps the genotype matrix when every SNP passes QC. Removing any SNP still
allocates a filtered matrix. File/DataFrame imports retain their existing loading
and intermediate allocations.

Use quality_control=false only when missing genotypes and marker filtering have
already been handled; it skips JWAS's imputation and marker QC checks.

BayesC automatically enables shared_genotypes for one dense genotype term whose
initial Matrix precision matches the requested run precision, in a standard
non-single-step/non-RRM analysis with unit residual weights.
Training and prediction selections must each be contiguous and ordered in
genotype order, without duplicates. Other methods and unsupported model modes
retain upstream copying defaults. Use shared_genotypes=false for noncontiguous
or repeated records. Traits can have different missing phenotypes within the
common training rows; prediction-only animals can have all traits missing.

```julia
# X is Matrix{Float32}; rows are ordered training first, prediction-only last.
geno = get_genotypes(X, G; method="BayesC", quality_control=false,
    double_precision=false)
# Set geno.obsID and geno.markerID to the input IDs.
model = build_model("y1 = intercept + geno\ny2 = intercept + geno", R)
outputEBV(model, prediction_ids)
# Omitting double_precision inherits the genotype matrix's Float32 precision.
result = runMCMC(model, phenotypes;
    fast_blocks=100, output_folder="bayesc_shared")
report = JWAS.genotype_storage_report(only(model.M))
@assert report.backing_matrix_count == 1
```

Use fast_blocks=false for the ordinary sampler and multi_trait_sampler=:II
in get_genotypes to select multi-trait sampler II. Set single-/multi-trait
priors G, R and Pi as required by the existing upstream API.

## Precision and conversion

get_genotypes defaults to double_precision=false, producing Float32 genotype
storage. When double_precision is omitted from runMCMC, the run inherits the
stored genotype precision: all Float32 terms select false and all Float64 terms
select true. Thus, for a Float64 matrix, set double_precision=true at import and
omit it from runMCMC; supported shared views can retain the imported parent
without a precision-conversion copy.

Models without genotype terms retain false as the default. Streaming genotype
storage uses Float32 and also selects false. Models with a mixture of Float32
and Float64 genotype terms require an explicit double_precision choice rather
than an implicit conversion.

Explicit true or false keeps its existing behavior. In particular, explicitly
requesting double_precision=true for an imported Float32 matrix disables
automatic sharing and uses the guarded Float64 conversion and copying path.
Explicit shared_genotypes=true still rejects a precision mismatch before
conversion.

With a Float32 import and the inferred run setting, genotype storage, marker
effect arrays and the main residual buffer use Float32. The analysis uses mixed
arithmetic: inclusion probabilities, annotation state and some multi-trait
workspaces use Float64. Multi-trait mixed-model equations can also use Float64
even when marker effects remain Float32. Sampled BayesC marker and residual
variances are cast back to Float32 when the resolved run setting is false.
Inferring the run setting preserves the existing arithmetic of the selected
precision mode; it does not make every computation use one floating-point type.

## Residual updates

BayesC blocks always update residuals in place using BLAS multipliers of the
genotype block's element type. This applies to single- and multi-trait sequential
and independent blocks, with no separate residual-mode option. Multi-trait
sequential updates copy only the current alpha block and need no product workspace.
The residual buffer keeps its existing element type; the update does not convert
or replace it. Ordinary marker updates retain their existing arithmetic.

Different floating-point rounding from upstream can alter seeded stochastic
trajectories. Prediction formulas, priors and sampling schedules are unchanged.

## Memory accounting and limits

Views keep the entire genotype parent alive. Memory reports count unique dense
genotype parents and block Gram caches separately; they do not measure peak
process RSS or all other model state. Dense memory estimates use the actual
stored element type, and the existing guard after alignment includes the full
retained parent for shared genotypes.

Before a requested Float64 genotype conversion, an additional check counts the
known unique retained dense genotype parents plus the new Float64 destination.
It uses memory_guard and memory_guard_ratio: :error stops before allocating the
destination when the estimate exceeds the threshold, :warn reports the estimate
and continues, and :off skips enforcement. This check covers known genotype
buffers rather than total available memory or peak RSS.

For example, 50,000 individuals by one million SNPs require 200 GB for one
Float32 matrix or 400 GB for one Float64 matrix (decimal units). Converting the
Float32 matrix while retaining its source requires at least 600 GB for those
two buffers alone. Loading intermediates, actual SNP filtering, alignment
copies, block caches and other model state can increase the total.

Neither guard protects the earlier dense file import or its intermediate
allocations from running out of memory. Efficient CSV ingestion and avoiding
copies when SNPs are actually removed require separate changes. The shared
dense path does not provide out-of-core storage.
