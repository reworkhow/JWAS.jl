# Multi-Trait BayesC Plain vs Empty-Annotated Sampler Benchmark Implementation

## Summary

Extended the production benchmark harness so it can compare four explicit
multi-trait cases on the packaged `simulated_annotations` fixture:

- `MT_BayesC_I`
- `MT_BayesC_II`
- `MT_EmptyAnnotated_BayesC_I`
- `MT_EmptyAnnotated_BayesC_II`

The implementation also adds:

- a focused `:plain_empty_sampler` selector
- an empty-annotation mode that activates the annotated path with an
  intercept-only design
- low-level plain multi-trait sampler-I versus sampler-II exact coverage

This report/implementation pair is written to answer the concrete overnight
question: whether sampler I and sampler II still match on a realistic
multi-trait dataset, both for plain multi-trait BayesC and for an
empty-annotated annotated-BayesC path.

Final answer from the production reruns:

- plain MT BayesC sampler I and II match on the headline shared top-8 metric
  on this fixture
- empty-annotated sampler I and II do not fully match in finite-chain
  production runs, and the direction of the gap can flip as the chain gets
  longer
- the exact tiny-case regression still supports the same target posterior for
  both samplers, so the remaining production gap is interpreted as a practical
  ranking/mixing effect rather than a target mismatch

## Files Changed

- Modified `benchmarks/simulated_annotations_multitrait_comparison.jl`
- Modified `test/unit/test_misc_coverage.jl`
- Modified `test/unit/test_multitrait_mcmc.jl`
- Added report:
  - `benchmarks/reports/2026-04-11-multitrait-bayesc-plain-empty-sampler-report.md`
- Added this implementation note:
  - `docs/plans/2026-04-11-multitrait-bayesc-plain-empty-sampler-implementation.md`

## What Changed

### Benchmark harness

In `benchmarks/simulated_annotations_multitrait_comparison.jl`:

- added explicit method cases:
  - `MT_BayesC_I`
  - `MT_BayesC_II`
  - `MT_EmptyAnnotated_BayesC_I`
  - `MT_EmptyAnnotated_BayesC_II`
- added `annotation_mode` per case
- added `annotation_mode_matrix(...)`
- added focused selector `selected_method_cases(:plain_empty_sampler)`
- added `JWAS_SIMULATED_MT_FOCUS_MODE`

### Empty annotation behavior

The first attempt used zero-valued annotation columns and failed validation
because JWAS rejects constant columns after adding the annotation intercept.

The final implementation uses an annotation matrix with zero user columns:

- shape `nmarkers x 0`
- JWAS still prepends the intercept internally
- the annotated BayesC path remains active
- no marker-level annotation signal is supplied

This is the intended meaning of the `Empty annotated BayesC` cases in the
report.

### Test coverage

In `test/unit/test_multitrait_mcmc.jl`:

- added exact plain multi-trait BayesC sampler-I versus sampler-II regression

In `test/unit/test_misc_coverage.jl`:

- updated the benchmark contract for the four-case selector
- added expectations for the new variants and sampler metadata

## Benchmark Runs

### Smoke

```bash
JWAS_SIMULATED_MT_FOCUS_MODE=plain_empty_sampler \
JWAS_SIMULATED_MT_SEEDS=11 \
JWAS_SIMULATED_MT_CHAIN_LENGTH=200 \
JWAS_SIMULATED_MT_BURNIN=50 \
JWAS_SIMULATED_MT_OUTPUT_FREQ=50 \
JWAS_SIMULATED_MT_WARMUP=false \
julia --project=. --startup-file=no \
  benchmarks/simulated_annotations_multitrait_comparison.jl \
  /tmp/jwas_mt_plain_empty_sampler_smoke
```

### First production rerun

```bash
JWAS_SIMULATED_MT_FOCUS_MODE=plain_empty_sampler \
JWAS_SIMULATED_MT_SEEDS=101,202,303,404 \
JWAS_SIMULATED_MT_CHAIN_LENGTH=4000 \
JWAS_SIMULATED_MT_BURNIN=1000 \
JWAS_SIMULATED_MT_OUTPUT_FREQ=40 \
JWAS_SIMULATED_MT_WARMUP=false \
julia --project=. --startup-file=no \
  benchmarks/simulated_annotations_multitrait_comparison.jl \
  /tmp/jwas_mt_plain_empty_sampler_20260411
```

### Longer warmup-enabled sensitivity rerun

```bash
JWAS_SIMULATED_MT_FOCUS_MODE=plain_empty_sampler \
JWAS_SIMULATED_MT_SEEDS=101,202 \
JWAS_SIMULATED_MT_CHAIN_LENGTH=12000 \
JWAS_SIMULATED_MT_BURNIN=3000 \
JWAS_SIMULATED_MT_OUTPUT_FREQ=120 \
JWAS_SIMULATED_MT_WARMUP=true \
julia --project=. --startup-file=no \
  benchmarks/simulated_annotations_multitrait_comparison.jl \
  /tmp/jwas_mt_plain_empty_sampler_long_20260411
```

### Stronger four-seed warmup-enabled rerun

```bash
JWAS_SIMULATED_MT_FOCUS_MODE=plain_empty_sampler \
JWAS_SIMULATED_MT_SEEDS=101,202,303,404 \
JWAS_SIMULATED_MT_CHAIN_LENGTH=12000 \
JWAS_SIMULATED_MT_BURNIN=3000 \
JWAS_SIMULATED_MT_OUTPUT_FREQ=120 \
JWAS_SIMULATED_MT_WARMUP=true \
julia --project=. --startup-file=no \
  benchmarks/simulated_annotations_multitrait_comparison.jl \
  /tmp/jwas_mt_plain_empty_sampler_long4_20260411
```

### Empty-only 12k control

```bash
JWAS_SIMULATED_MT_FOCUS_MODE=empty_annotated_sampler \
JWAS_SIMULATED_MT_SEEDS=101,202,303,404 \
JWAS_SIMULATED_MT_CHAIN_LENGTH=12000 \
JWAS_SIMULATED_MT_BURNIN=3000 \
JWAS_SIMULATED_MT_OUTPUT_FREQ=120 \
JWAS_SIMULATED_MT_WARMUP=true \
julia --project=. --startup-file=no \
  benchmarks/simulated_annotations_multitrait_comparison.jl \
  /tmp/jwas_mt_empty_only_12k_20260411
```

### Empty-only 24k longer run

```bash
JWAS_SIMULATED_MT_FOCUS_MODE=empty_annotated_sampler \
JWAS_SIMULATED_MT_SEEDS=101,202,303,404 \
JWAS_SIMULATED_MT_CHAIN_LENGTH=24000 \
JWAS_SIMULATED_MT_BURNIN=6000 \
JWAS_SIMULATED_MT_OUTPUT_FREQ=120 \
JWAS_SIMULATED_MT_WARMUP=true \
julia --project=. --startup-file=no \
  benchmarks/simulated_annotations_multitrait_comparison.jl \
  /tmp/jwas_mt_empty_only_longer_20260411
```

## Final Findings

From the strongest rerun at `/tmp/jwas_mt_plain_empty_sampler_long4_20260411`:

- `MT_BayesC_I`
  - runtime mean `17.70s`
  - shared precision/recall/F1 `0.3750`
  - mean true shared recovered inside top-8 shared set `3.00`
- `MT_BayesC_II`
  - runtime mean `69.02s`
  - shared precision/recall/F1 `0.3750`
  - mean true shared recovered inside top-8 shared set `3.00`
- `MT_EmptyAnnotated_BayesC_I`
  - runtime mean `17.97s`
  - shared precision/recall/F1 `0.53125`
  - mean true shared recovered inside top-8 shared set `4.25`
- `MT_EmptyAnnotated_BayesC_II`
  - runtime mean `66.85s`
  - shared precision/recall/F1 `0.46875`
  - mean true shared recovered inside top-8 shared set `3.75`

The marker-level follow-up showed:

- plain sampler I and II still have the same headline shared count even though
  their exact top-8 shared sets differ
- the empty-annotated mismatch is caused by a few boundary markers near the
  top-8 `P11` cutoff, where sampler I more often keeps one extra weak true
  shared marker above the line

The empty-only controls showed:

- the `12k` empty-only rerun matched the earlier `12k` mixed-family result
  exactly, so the empty-family gap is not a benchmark-order artifact
- the `24k` empty-only rerun flipped the `12k` ordering:
  - `MT_EmptyAnnotated_BayesC_I`: shared F1 `0.46875`
  - `MT_EmptyAnnotated_BayesC_II`: shared F1 `0.53125`

So the empty-family sampler difference is chain-length sensitive and not yet
stabilized by the current production settings.

## Verification

Executed so far:

- `julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'`
- `julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'`
- smoke benchmark above
- first production rerun above
- longer warmup-enabled sensitivity rerun above
- stronger four-seed warmup-enabled rerun above

Final full-suite verification after the benchmark/report updates:

- `julia --project=. --startup-file=no -e 'include("test/unit/test_multitrait_mcmc.jl")'`
- `julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'`
- `julia --project=. --startup-file=no test/runtests.jl`
