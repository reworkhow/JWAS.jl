# Annotated BayesR Benchmark Report

## Goal

Benchmark dense single-trait annotated BayesR on the production JWAS path and
compare it against:

- dense BayesR
- dense annotated BayesC

on a synthetic dataset where annotations truly enrich nonzero and larger-effect
markers.

## Protocol

Production benchmark script:

- `benchmarks/annotated_bayesr_comparison.jl`

Output directory used for the full run:

- `/tmp/annotated_bayesr_benchmark_20260328`

Methods:

- `BayesR`
- `Annotated_BayesR`
- `Annotated_BayesC`

Chain settings:

- `chain_length = 10000`
- `burnin = 2000`
- `output_samples_frequency = 10`
- MCMC seeds `2026, 2027, 2028, 2029, 2030`

Synthetic dataset:

- `n_obs = 200`
- `n_markers = 1000`
- two annotation columns
- target trait heritability `0.45`

## Data-Generating Truth

The dataset was built so that `annotation_1` is informative and `annotation_2`
is nuisance.

Truth summary:

- total markers: `1000`
- enriched markers (`annotation_1 = 1`): `250`
- baseline markers (`annotation_1 = 0`): `750`
- causal markers: `71`
- class counts:
  - `class1 = 929`
  - `class2 = 44`
  - `class3 = 18`
  - `class4 = 9`

Enrichment was real in the data-generating process:

- baseline causal rate: `0.0373`
- enriched causal rate: `0.1720`
- baseline large-class rate (`class3` or `class4`): `0.0107`
- enriched large-class rate: `0.0760`

So the benchmark does contain a real annotation signal, but it is still a
sparse scenario with only `27` markers in classes `3` and `4` combined.

## Method-Level Summary

Mean results across the five seeds:

| Method | cor(y, EBV) | mean PIP causal | mean PIP noncausal | enriched PIP | baseline PIP | top-causal recall |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `BayesR` | `0.9063` | `0.3891` | `0.3499` | `0.3630` | `0.3493` | `0.1437` |
| `Annotated_BayesR` | `0.7942` | `0.1127` | `0.0248` | `0.0953` | `0.0096` | `0.2282` |
| `Annotated_BayesC` | `0.8928` | `0.5114` | `0.3642` | `0.6223` | `0.2921` | `0.1775` |

## What Worked

### 1. Annotated BayesR learned the direction of the informative annotation

Posterior mean annotation coefficients for annotated BayesR:

- step 1 (`zero vs nonzero`)
  - `Annotation_1 = +1.6286`
  - `Annotation_2 = -0.2220`
- step 2 (`small vs larger`)
  - `Annotation_1 = +0.3944`
  - `Annotation_2 = +0.0405`
- step 3 (`medium vs large`)
  - `Annotation_1 = -0.1852`
  - `Annotation_2 = -0.1948`

The important point is step 1:

- the informative annotation got a clearly positive coefficient
- the nuisance annotation stayed near zero

So the model did detect the main annotation signal that was built into the
dataset.

### 2. Annotated BayesR improved PIP separation strongly

BayesR barely separated enriched from baseline markers:

- enriched mean PIP: `0.3630`
- baseline mean PIP: `0.3493`
- gap: `0.0137`

Annotated BayesR separated them much more strongly:

- enriched mean PIP: `0.0953`
- baseline mean PIP: `0.0096`
- gap: `0.0857`

It also separated causal from noncausal markers much more strongly:

- BayesR causal minus noncausal PIP gap: `0.0392`
- Annotated BayesR causal minus noncausal PIP gap: `0.0879`

### 3. Annotated BayesR improved class-order separation

Mean posterior PIP by true class:

| Method | class1 | class2 | class3 | class4 |
| --- | ---: | ---: | ---: | ---: |
| `BayesR` | `0.3499` | `0.3459` | `0.3647` | `0.6496` |
| `Annotated_BayesR` | `0.0248` | `0.0429` | `0.1057` | `0.4683` |

This is a real improvement in ranking behavior:

- ordinary BayesR gave very similar PIP to classes `1`, `2`, and `3`
- annotated BayesR produced a much more ordered progression from class `1` to class `4`

So annotated BayesR is doing something meaningful in terms of prioritization,
not just collapsing randomly.

## What Did Not Work Well

### 1. Annotated BayesR hurt EBV correlation on this benchmark

This was the clearest negative result:

- `BayesR`: `0.9063`
- `Annotated_BayesR`: `0.7942`

So even though annotated BayesR learned the direction of the annotation and
improved marker prioritization, it produced worse in-sample `cor(y, EBV)` on
this sparse synthetic scenario.

### 2. Annotated BayesR was too conservative overall

Its overall inclusion level was low:

- mean PIP across all markers: `0.0310`

Compared with:

- `BayesR`: `0.3527`

This matters because the model was not just separating better. It was also
shrinking almost everything toward zero, including many truly nonzero markers.

### 3. Upper-class conditional models were numerically weak

The step-2 and step-3 intercepts were very large:

- step 2 intercept mean: `19.28`
- step 3 intercept mean: `28.13`

That is a sign the conditional models for the rarer upper classes are close to
degenerate in this scenario. The truth only had:

- `18` markers in class `3`
- `9` markers in class `4`

So the upper conditional annotation models are being estimated from very little
effective information.

## Comparison With Annotated BayesC

Annotated BayesC behaved differently:

- better EBV correlation than annotated BayesR: `0.8928` vs `0.7942`
- very strong enriched-vs-baseline PIP gap: `0.3302`
- but much larger between-seed variability

Its annotation coefficient summary also recovered the informative annotation in
the expected direction:

- `Annotation_1` mean coefficient: `+5.48`
- `Annotation_2` mean coefficient: `-0.38`

But the between-seed instability was substantial:

- EBV correlation SD: `0.0772`
- mean PIP SD across seeds: `0.3877`

So annotated BayesC looked more aggressive than annotated BayesR and often more
useful on this benchmark, but also less stable.

## Interpretation

The production implementation itself looks sound:

- the benchmark runs on the production JWAS path
- annotated BayesR learns the informative annotation at step 1
- nuisance annotation stays near zero
- marker ranking and class-order separation improve relative to ordinary BayesR

But the current sparse synthetic scenario exposed a scientific limitation:

- annotated BayesR is trying to estimate three conditional annotation models
- the upper two conditionals are informed by very few markers
- the result is over-conservative nonzero assignment and lower EBV correlation

So the correct reading is not:

- "annotated BayesR is broken"

It is:

- "annotated BayesR is working mechanically, but this first benchmark does not
  show a predictive advantage under a sparse truth with very rare upper classes"

## Limitations

This report is based on one synthetic scenario:

- one genotype dataset
- one annotation design
- one truth regime

So it should not be over-generalized. The next benchmark that would matter most
is:

- a second synthetic scenario with a higher prevalence of nonzero and upper
  classes, so the step-2 and step-3 annotation models are less sparse

## Bottom Line

For this benchmark:

- annotated BayesR successfully learned the informative annotation direction
- annotated BayesR improved SNP prioritization and class-order separation
- annotated BayesR did **not** improve EBV correlation relative to ordinary
  BayesR
- the likely reason is that the upper-class conditional annotation models are
  too weakly identified in this sparse scenario

So the current production annotated BayesR implementation is benchmarked and
working, but it is not yet benchmark-proven to outperform ordinary BayesR on
prediction quality.
