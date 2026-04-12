# Multi-Trait BayesC Plain vs Empty-Annotated Sampler Report

## Format

This report uses the following presentation format for the packaged
`simulated_annotations` 2-trait fixture.

| Requested label | Benchmark variant |
| --- | --- |
| `Mod MTBase-C sampler 1` | `MT_BayesC_I` |
| `Mod MTBase-C sampler 2` | `MT_BayesC_II` |
| `Empty annotated BayesC sampler 1` | `MT_EmptyAnnotated_BayesC_I` |
| `Empty annotated BayesC sampler 2` | `MT_EmptyAnnotated_BayesC_II` |

`Empty annotated` does not mean ordinary multi-trait BayesC. It means the
annotated BayesC path is active with an intercept-only annotation design:

- user annotation matrix has zero columns
- JWAS still adds the annotation intercept internally
- prior updating still uses the annotated 3-step probit tree

So the plain and empty-annotated cases are different production models even
though the user-supplied marker annotations are empty.

## Benchmark

Focused production benchmark on the packaged `simulated_annotations` 2-trait
fixture using the production harness in
`benchmarks/simulated_annotations_multitrait_comparison.jl`.

Primary metric:

- pleiotropic recovery from the true joint posterior shared-state score
  `P11 = Pr(δ = (1,1) | data)`

Secondary metrics:

- genomic prediction accuracy `cor(y, EBV)` for each trait
- per-trait top-`k` causal recovery
- any-active top-`k` recovery

The packaged truth contains:

- `944` null markers
- `6` trait-1-only markers
- `6` trait-2-only markers
- `8` shared markers

So each trait has `14` active markers and there are `20` any-active markers.

## Exact-Posterior Backstop

The low-level regression in `test/unit/test_multitrait_mcmc.jl` now covers:

- plain multi-trait BayesC sampler I versus II
- annotated multi-trait BayesC sampler I versus II

In both tiny one-marker cases, sampler I and sampler II agree with the same
exact posterior target within Monte Carlo tolerance.

That backstop means any production differences below should be interpreted as:

- finite-chain mixing/efficiency differences, or
- plain versus empty-annotated model differences

not as evidence that sampler I and sampler II target different posteriors.

## Main Production Rerun

Cold production rerun:

- output directory: `/tmp/jwas_mt_plain_empty_sampler_20260411`
- seeds: `101`, `202`, `303`, `404`
- `chain_length = 4000`
- `burnin = 1000`
- `output_samples_frequency = 40`
- warmup: disabled

### Compact summary by method

| Requested label | Variant | Runtime mean (s) | Mean `cor(y, EBV)` | Any-active recall mean | Shared precision/recall/F1 mean |
| --- | --- | ---: | ---: | ---: | ---: |
| `Mod MTBase-C sampler 1` | `MT_BayesC_I` | `14.37` | `0.7518` | `0.3625` | `0.3750` |
| `Mod MTBase-C sampler 2` | `MT_BayesC_II` | `26.19` | `0.7486` | `0.3000` | `0.3750` |
| `Empty annotated BayesC sampler 1` | `MT_EmptyAnnotated_BayesC_I` | `6.91` | `0.7425` | `0.3500` | `0.4688` |
| `Empty annotated BayesC sampler 2` | `MT_EmptyAnnotated_BayesC_II` | `23.36` | `0.7494` | `0.3875` | `0.4063` |

### Pleiotropic recovery

Shared recovery uses `P11` top-`8`:

- rank markers by reconstructed `P11`
- declare the top `8` markers as shared
- compare to the true `8` shared loci

| Requested label | True shared inside top 8 mean | Mean `P11` on true shared | Mean `P11` on nonshared |
| --- | ---: | ---: | ---: |
| `Mod MTBase-C sampler 1` | `3.00` | `0.3696` | `0.0149` |
| `Mod MTBase-C sampler 2` | `3.00` | `0.3733` | `0.0124` |
| `Empty annotated BayesC sampler 1` | `3.75` | `0.4208` | `0.0163` |
| `Empty annotated BayesC sampler 2` | `3.25` | `0.4279` | `0.0218` |

Per-seed shared counts inside the top-8 shared set:

| Variant | Seed 101 | Seed 202 | Seed 303 | Seed 404 |
| --- | ---: | ---: | ---: | ---: |
| `MT_BayesC_I` | `3` | `3` | `3` | `3` |
| `MT_BayesC_II` | `3` | `3` | `3` | `3` |
| `MT_EmptyAnnotated_BayesC_I` | `4` | `3` | `4` | `4` |
| `MT_EmptyAnnotated_BayesC_II` | `3` | `3` | `3` | `4` |

### Main reading from the first production run

1. Plain multi-trait BayesC sampler I and II are effectively aligned on the
   primary shared metric.
   Both plain samplers recovered exactly `3` shared loci inside the top-8
   shared set for all four seeds, with shared precision/recall/F1 all equal to
   `0.3750`.

2. The remaining practical difference in this dataset shows up in the
   empty-annotated path, not the plain path.
   `MT_EmptyAnnotated_BayesC_I` outperformed `MT_EmptyAnnotated_BayesC_II` on
   shared recovery in the first production rerun:
   `0.4688` versus `0.4063`.

3. The runtime numbers from this cold run should not be treated as final
   sampler-performance results.
   The first case in the process pays JIT compilation overhead, and that makes
   the raw runtime comparison less clean than the recovery comparison.

## Longer-Chain Sensitivity Rerun

Longer warmup-enabled sensitivity rerun:

- output directory: `/tmp/jwas_mt_plain_empty_sampler_long_20260411`
- seeds: `101`, `202`
- `chain_length = 12000`
- `burnin = 3000`
- `output_samples_frequency = 120`
- warmup: enabled

Two-seed long-run result:

| Requested label | Variant | Shared precision/recall/F1 mean | Mean `P11` on true shared | Mean `P11` on nonshared |
| --- | --- | ---: | ---: | ---: |
| `Mod MTBase-C sampler 1` | `MT_BayesC_I` | `0.3750` | `0.3733` | `0.0166` |
| `Mod MTBase-C sampler 2` | `MT_BayesC_II` | `0.3750` | `0.3733` | `0.0137` |
| `Empty annotated BayesC sampler 1` | `MT_EmptyAnnotated_BayesC_I` | `0.5000` | `0.4600` | `0.0211` |
| `Empty annotated BayesC sampler 2` | `MT_EmptyAnnotated_BayesC_II` | `0.4375` | `0.4383` | `0.0221` |

The two-seed long run kept the same qualitative pattern as the first production
rerun:

- plain sampler I and II still matched on the primary shared metric
- empty-annotated sampler I still outperformed empty-annotated sampler II on
  the shared top-8 summary

## Stronger Four-Seed Warmup-Enabled Rerun

Final stronger production rerun:

- output directory: `/tmp/jwas_mt_plain_empty_sampler_long4_20260411`
- seeds: `101`, `202`, `303`, `404`
- `chain_length = 12000`
- `burnin = 3000`
- `output_samples_frequency = 120`
- warmup: enabled

### Final compact summary by method

| Requested label | Variant | Runtime mean (s) | Mean `cor(y, EBV)` | Any-active recall mean | Shared precision/recall/F1 mean |
| --- | --- | ---: | ---: | ---: | ---: |
| `Mod MTBase-C sampler 1` | `MT_BayesC_I` | `17.70` | `0.7554` | `0.3250` | `0.3750` |
| `Mod MTBase-C sampler 2` | `MT_BayesC_II` | `69.02` | `0.7489` | `0.3375` | `0.3750` |
| `Empty annotated BayesC sampler 1` | `MT_EmptyAnnotated_BayesC_I` | `17.97` | `0.7445` | `0.3000` | `0.5313` |
| `Empty annotated BayesC sampler 2` | `MT_EmptyAnnotated_BayesC_II` | `66.85` | `0.7467` | `0.3750` | `0.4688` |

### Final pleiotropic recovery summary

| Requested label | True shared inside top 8 mean | Mean `P11` on true shared | Mean `P11` on nonshared |
| --- | ---: | ---: | ---: |
| `Mod MTBase-C sampler 1` | `3.00` | `0.3692` | `0.0161` |
| `Mod MTBase-C sampler 2` | `3.00` | `0.3763` | `0.0129` |
| `Empty annotated BayesC sampler 1` | `4.25` | `0.4600` | `0.0241` |
| `Empty annotated BayesC sampler 2` | `3.75` | `0.4533` | `0.0239` |

Per-seed shared counts inside the top-8 shared set:

| Variant | Seed 101 | Seed 202 | Seed 303 | Seed 404 |
| --- | ---: | ---: | ---: | ---: |
| `MT_BayesC_I` | `3` | `3` | `3` | `3` |
| `MT_BayesC_II` | `3` | `3` | `3` | `3` |
| `MT_EmptyAnnotated_BayesC_I` | `4` | `4` | `5` | `4` |
| `MT_EmptyAnnotated_BayesC_II` | `3` | `4` | `4` | `4` |

### Main reading from the strongest rerun

1. `Mod MTBase-C sampler 1` and `Mod MTBase-C sampler 2` still match on the
   primary pleiotropic-recovery metric.
   Across all four seeds, both plain samplers recovered exactly `3` shared
   loci inside the top-8 shared set.

2. The remaining mismatch persists inside the empty-annotated model family.
   In the strongest rerun, `MT_EmptyAnnotated_BayesC_I` recovered `4.25`
   shared loci on average inside the top-8 shared set, versus `3.75` for
   `MT_EmptyAnnotated_BayesC_II`.

3. Sampler II is materially slower on this production path.
   In the strongest rerun, sampler-II runtime was about `3.8x` the sampler-I
   runtime in both the plain and empty-annotated families.

## Marker-Level Diagnosis

The final longer run shows that the remaining empty-annotated mismatch is not a
global posterior failure. It is a ranking difference concentrated near the
top-8 shared cutoff.

Observed marker-level diagnostics from the final rerun:

- plain `MT_BayesC_I` versus `MT_BayesC_II`
  - per-seed `P11` correlation: about `0.93`
  - top-8 overlap: only `3` to `4` markers
  - interpretation: the samplers swap mostly false positives, so the shared
    count metric stays identical even though the exact top-8 sets differ
- empty-annotated `MT_EmptyAnnotated_BayesC_I` versus
  `MT_EmptyAnnotated_BayesC_II`
  - per-seed `P11` correlation: about `0.48` to `0.89`
  - top-8 overlap: `3` to `5` markers
  - interpretation: the samplers disagree more on borderline markers in this
    intercept-only annotated model

Two concrete seeds show the pattern clearly:

- seed `101`
  - sampler I kept true shared marker `m390` at `P11 = 0.2933`
  - sampler II pushed `m390` down to `P11 = 0.0133`
  - sampler II filled the boundary with false positives such as `m200`,
    `m396`, and `m662`
- seed `303`
  - sampler I kept true shared marker `m411` above the cutoff
  - sampler II assigned more mass to false positives such as `m200` and `m458`

So the persistent gap is a practical ranking effect:

- sampler I more often keeps one extra weak true shared locus above the top-8
  boundary in the empty-annotated family
- sampler II more often swaps in a near-boundary false positive

## Empty-Only Follow-Up

To check whether the empty-annotated gap was being distorted by benchmark
composition, two additional empty-family-only production reruns were run.

### 12k empty-only control

Settings:

- output directory: `/tmp/jwas_mt_empty_only_12k_20260411`
- seeds: `101`, `202`, `303`, `404`
- `chain_length = 12000`
- `burnin = 3000`
- `output_samples_frequency = 120`
- warmup: enabled

Result:

- `MT_EmptyAnnotated_BayesC_I`: shared precision/recall/F1 `0.53125`
- `MT_EmptyAnnotated_BayesC_II`: shared precision/recall/F1 `0.46875`

This matches the earlier mixed-family 12k benchmark exactly.

So the earlier 12k result was not a case-order artifact. Running only the empty
family reproduces the same answer.

### 24k empty-only longer run

Settings:

- output directory: `/tmp/jwas_mt_empty_only_longer_20260411`
- seeds: `101`, `202`, `303`, `404`
- `chain_length = 24000`
- `burnin = 6000`
- `output_samples_frequency = 120`
- warmup: enabled

Result:

- `MT_EmptyAnnotated_BayesC_I`: shared precision/recall/F1 `0.46875`
- `MT_EmptyAnnotated_BayesC_II`: shared precision/recall/F1 `0.53125`

This longer run flips the 12k ordering.

That is the strongest evidence from the overnight work:

- the empty-annotated sampler-I versus sampler-II gap is real as a finite-chain
  production phenomenon
- but it does not stabilize monotonically with the current chain lengths
- at `12k` post-benchmark settings, sampler I looks better
- at `24k`, sampler II looks better

So the correct conclusion is not that sampler I or sampler II is definitively
better in the empty-annotated family. The correct conclusion is that this
family remains chain-length sensitive on the shared top-8 metric.

## Interpretation

### 1. Plain versus empty-annotated is a model comparison, not only a sampler comparison.

`MT_BayesC_I` and `MT_BayesC_II` are ordinary multi-trait BayesC with the
standard global joint-`Pi` update.

`MT_EmptyAnnotated_BayesC_I` and `MT_EmptyAnnotated_BayesC_II` still use the
annotated BayesC prior-update path:

- intercept-only annotation design
- 3-step probit tree update
- reconstructed joint marker priors

So if plain and empty-annotated results differ, that is expected in principle.
They are not the same posterior model with a different input matrix; they are
different prior-update families.

### 2. Sampler-I versus sampler-II should be judged within a model family.

The exact tiny-case regressions show that:

- plain sampler I and II share the same target posterior
- annotated sampler I and II share the same target posterior

So the real production question is whether finite-chain behavior differs enough
to matter in practice.

### 3. The reruns narrow the issue further than the first production pass.

On this dataset:

- plain sampler I and II already look matched on the primary shared metric
- the only material remaining gap is between empty-annotated sampler I and II
- that remaining gap survives both longer chains and additional seeds
- the 12k empty-only control shows the gap is not a benchmark-order artifact
- the 24k empty-only run shows the direction of the gap can flip with longer
  chains
- the gap is concentrated in boundary-marker ranking, not a global failure of
  posterior targeting

That is why the final interpretation is:

- plain sampler I versus II: effectively matched for the headline shared metric
  on this fixture
- empty-annotated sampler I versus II: same target in the exact tiny-case test,
  but a persistent practical finite-chain ranking difference on this fixture

## Conclusion

The overnight production work now supports three concrete conclusions:

1. `Mod MTBase-C sampler 1` and `Mod MTBase-C sampler 2` match on the primary
   pleiotropic-recovery metric for this packaged dataset.
   They do not return identical top-8 marker sets, but they recover the same
   number of true shared loci in every production rerun we ran.

2. The remaining mismatch is inside the empty-annotated model family.
   At `12k`, `Empty annotated BayesC sampler 1` outperformed
   `Empty annotated BayesC sampler 2` on the shared top-8 metric.
   At `24k`, that ordering reversed.

3. The empty-annotated sampler gap looks like a practical finite-chain ranking
   effect, not evidence of different target posteriors.
   The exact tiny-case regression still shows the same target posterior, while
   the production mismatch comes from a few weak shared markers and false
   positives swapping around the top-8 boundary.

4. The empty-annotated family is not yet benchmark-stable enough for a
   definitive sampler-I versus sampler-II recommendation on this fixture.
   If we want a firmer answer, the next work should be threshold-free shared
   metrics or still longer chains / multiple independent chains.
