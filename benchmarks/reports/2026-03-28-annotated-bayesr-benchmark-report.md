# Annotated BayesR Benchmark Report

## Goal

Benchmark dense single-trait annotated BayesR on the production JWAS path and
compare it against:

- dense BayesR
- dense annotated BayesC

under two synthetic regimes:

- a sparse upper-class regime
- a less-sparse upper-class regime

The main question is whether annotated BayesR improves SNP prioritization and
prediction when the annotations truly enrich nonzero and larger-effect markers.

## Protocol

Production benchmark script:

- `benchmarks/annotated_bayesr_comparison.jl`

Output directories:

- sparse scenario: `/tmp/annotated_bayesr_benchmark_20260328`
- less-sparse scenario: `/tmp/annotated_bayesr_benchmark_20260328_less_sparse`

Methods:

- `BayesR`
- `Annotated_BayesR`
- `Annotated_BayesC`

Common settings:

- single trait
- dense storage only
- `n_obs = 200`
- `n_markers = 1000`
- two annotation columns
- target trait heritability `0.45`
- `chain_length = 10000`
- `burnin = 2000`
- `output_samples_frequency = 10`
- seeds `2026, 2027, 2028, 2029, 2030`

## Data-Generating Truth

Both scenarios used the same annotation layout:

- `250` enriched SNPs (`annotation_1 = 1`)
- `750` baseline SNPs (`annotation_1 = 0`)
- `annotation_2` is nuisance

The only change was the true class probabilities assigned to enriched and
baseline SNPs.

### Scenario 1: `sparse_upper_classes`

Truth summary:

- total markers: `1000`
- causal markers: `71`
- class counts:
  - `class1 = 929`
  - `class2 = 44`
  - `class3 = 18`
  - `class4 = 9`

Rates:

- baseline causal rate: `0.0373`
- enriched causal rate: `0.1720`
- baseline large-class rate (`class3` or `class4`): `0.0107`
- enriched large-class rate: `0.0760`

This is a genuinely sparse setting. Only `27` markers are in classes `3` and
`4` combined.

### Scenario 2: `less_sparse_upper_classes`

Truth summary:

- total markers: `1000`
- causal markers: `161`
- class counts:
  - `class1 = 839`
  - `class2 = 76`
  - `class3 = 56`
  - `class4 = 29`

Rates:

- baseline causal rate: `0.1053`
- enriched causal rate: `0.3280`
- baseline large-class rate (`class3` or `class4`): `0.0560`
- enriched large-class rate: `0.1720`

This second regime gives the step-2 and step-3 annotation models much more real
signal than the sparse regime.

## Scenario 1 Results: `sparse_upper_classes`

Mean results across five seeds:

| Method | cor(y, EBV) | mean PIP causal | mean PIP noncausal | enriched PIP | baseline PIP | top-causal recall |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `BayesR` | `0.9063` | `0.3891` | `0.3499` | `0.3630` | `0.3493` | `0.1437` |
| `Annotated_BayesR` | `0.7942` | `0.1127` | `0.0248` | `0.0953` | `0.0096` | `0.2282` |
| `Annotated_BayesC` | `0.8928` | `0.5114` | `0.3642` | `0.6223` | `0.2921` | `0.1775` |

Main observations:

1. Annotated BayesR learned the informative annotation direction at step 1.
   - `Annotation_1 = +1.6286`
   - `Annotation_2 = -0.2220`
2. Annotated BayesR improved prioritization strongly.
   - causal minus noncausal mean PIP gap:
     - `BayesR = 0.0392`
     - `Annotated_BayesR = 0.0879`
   - enriched minus baseline mean PIP gap:
     - `BayesR = 0.0137`
     - `Annotated_BayesR = 0.0857`
3. Annotated BayesR improved class-order separation.
   - BayesR mean PIP by true class:
     - `class1 = 0.3499`
     - `class2 = 0.3459`
     - `class3 = 0.3647`
     - `class4 = 0.6496`
   - Annotated BayesR mean PIP by true class:
     - `class1 = 0.0248`
     - `class2 = 0.0429`
     - `class3 = 0.1057`
     - `class4 = 0.4683`
4. Annotated BayesR hurt EBV correlation in this sparse setting.
   - `BayesR = 0.9063`
   - `Annotated_BayesR = 0.7942`

Interpretation:

- The implementation was mechanically working.
- The model used the informative annotation.
- But the upper conditional models were weakly identified because the truth only
  had `18` class-3 markers and `9` class-4 markers.
- The step-2 and step-3 intercepts became very large:
  - step 2 intercept mean: `19.2796`
  - step 3 intercept mean: `28.1303`

So the sparse scenario mostly showed that annotated BayesR can improve ranking,
but not that it improves prediction.

## Scenario 2 Results: `less_sparse_upper_classes`

Mean results across five seeds:

| Method | cor(y, EBV) | mean PIP causal | mean PIP noncausal | enriched PIP | baseline PIP | top-causal recall |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `BayesR` | `0.9585` | `0.6046` | `0.5870` | `0.5904` | `0.5897` | `0.2683` |
| `Annotated_BayesR` | `0.8839` | `0.1363` | `0.1017` | `0.1431` | `0.0953` | `0.2944` |
| `Annotated_BayesC` | `0.9174` | `0.5680` | `0.5766` | `0.5473` | `0.5845` | `0.1714` |

Main observations:

1. Annotated BayesR still improved prioritization relative to BayesR.
   - causal minus noncausal mean PIP gap:
     - `BayesR = 0.0176`
     - `Annotated_BayesR = 0.0347`
   - enriched minus baseline mean PIP gap:
     - `BayesR = 0.0007`
     - `Annotated_BayesR = 0.0478`
2. Annotated BayesR again improved top-causal recall over BayesR.
   - `BayesR = 0.2683`
   - `Annotated_BayesR = 0.2944`
3. Annotated BayesR still preserved class ordering better than BayesR.
   - BayesR mean PIP by true class:
     - `class1 = 0.5870`
     - `class2 = 0.5903`
     - `class3 = 0.5941`
     - `class4 = 0.6623`
   - Annotated BayesR mean PIP by true class:
     - `class1 = 0.1017`
     - `class2 = 0.1183`
     - `class3 = 0.1236`
     - `class4 = 0.2083`
4. EBV correlation was still worse than BayesR, but the gap narrowed.
   - sparse gap: `0.1120`
   - less-sparse gap: `0.0746`

Annotation coefficient summary for annotated BayesR:

- step 1 (`zero vs nonzero`)
  - `Annotation_1 = +0.2711`
  - `Annotation_2 = -0.0038`
- step 2 (`small vs larger`)
  - `Annotation_1 = +0.1240`
  - `Annotation_2 = -0.0726`
- step 3 (`medium vs large`)
  - `Annotation_1 = -0.0500`
  - `Annotation_2 = -0.0133`

Interpretation:

- The step-1 informative annotation effect remained in the expected direction,
  but it was much smaller than in the sparse scenario.
- The upper-step coefficients stayed weak and noisy.
- The step-2 and step-3 intercepts were still very large:
  - step 2 intercept mean: `13.1159`
  - step 3 intercept mean: `18.9376`

So the less-sparse regime helped, but it did not fully stabilize the upper
conditional annotation models.

## What This Says About Annotated BayesR

Across both scenarios, the pattern was consistent.

### 1. The production implementation is working

The benchmark used only the production `runMCMC` path, and annotated BayesR
behaved coherently in both regimes:

- it learned the informative annotation at step 1
- it kept the nuisance annotation near zero
- it changed SNP prioritization in the expected direction

This is not a mechanics failure.

### 2. Annotated BayesR is better at prioritization than at prediction

In both scenarios:

- enriched vs baseline PIP separation improved
- causal vs noncausal PIP separation improved
- top-causal recall improved
- class-order separation improved

But in both scenarios:

- `cor(y, EBV)` stayed below ordinary BayesR

So the current implementation is stronger as a fine-mapping or prioritization
tool than as a predictor on these synthetic regimes.

### 3. The main difficulty is still the upper conditional models

Even after increasing class-3 and class-4 prevalence, the upper-step intercepts
were still large and the informative step-2/step-3 coefficients were weak.

That means the central scientific issue is not the first step. It is the
stability and usefulness of the deeper sequential annotation models.

## What This Says About Annotated BayesC

Annotated BayesC behaved very differently across the two regimes.

In the sparse regime:

- strong enriched-vs-baseline PIP separation: `0.3302`
- reasonable EBV correlation: `0.8928`

In the less-sparse regime:

- enriched-vs-baseline gap turned negative: `-0.0372`
- causal-vs-noncausal gap also turned negative: `-0.0086`
- informative annotation coefficient mean became negative:
  - `Annotation_1 = -1.0116`

So in this benchmark family, annotated BayesC is not the more reliable
annotation baseline. It looked aggressive in the sparse regime but unstable
across scenarios.

## Bottom Line

The two-scenario benchmark supports these conclusions.

1. Annotated BayesR is implemented correctly enough to learn the annotation and
   change posterior prioritization on the production JWAS path.
2. Annotated BayesR consistently improves marker prioritization relative to
   ordinary BayesR.
3. Annotated BayesR is not yet benchmark-proven to improve prediction
   (`cor(y, EBV)`) relative to ordinary BayesR.
4. The remaining scientific weakness is the upper conditional annotation models,
   not the basic annotated BayesR implementation.
5. Annotated BayesC is not a uniformly stronger benchmark comparator here; it
   became unstable when the truth regime changed.

So the current implementation should be viewed as:

- a sound v1 annotated BayesR sampler
- promising for SNP prioritization
- not yet validated as a superior predictor

The next benchmark that would matter most is a regime where the upper classes
carry more signal and where the annotation effects for steps 2 and 3 are
stronger by design. That would test whether the current limitations are due to
weak truth or to the model structure itself.
