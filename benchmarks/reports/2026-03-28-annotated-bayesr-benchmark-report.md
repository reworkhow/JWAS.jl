# Annotated BayesR Benchmark Report

## Goal

Benchmark dense single-trait annotated BayesR on the production JWAS path and
compare it against:

- dense BayesR
- dense annotated BayesC

under three synthetic regimes:

- a sparse upper-class regime
- a less-sparse upper-class regime
- a stepwise conditional-signal regime

The main question is whether annotated BayesR improves SNP prioritization and
prediction when the annotations truly enrich nonzero and larger-effect markers.

## Protocol

Production benchmark script:

- `benchmarks/annotated_bayesr_comparison.jl`

Output directories:

- sparse scenario: `/tmp/annotated_bayesr_benchmark_20260328`
- less-sparse scenario: `/tmp/annotated_bayesr_benchmark_20260328_less_sparse`
- stepwise-signal scenario: `/tmp/annotated_bayesr_benchmark_20260328_stepwise_signal`

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

All three scenarios used the same annotation layout:

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

### Scenario 3: `stepwise_annotation_signal`

This regime was designed directly in the sequential conditional parameterization
used by annotated BayesR.

Target conditional probabilities:

- baseline SNPs:
  - `P(delta > 1) = 0.10`
  - `P(delta > 2 | delta > 1) = 0.20`
  - `P(delta > 3 | delta > 2) = 0.20`
- enriched SNPs:
  - `P(delta > 1) = 0.30`
  - `P(delta > 2 | delta > 1) = 0.60`
  - `P(delta > 3 | delta > 2) = 0.60`

Realized truth summary:

- total markers: `1000`
- causal markers: `144`
- class counts:
  - `class1 = 856`
  - `class2 = 92`
  - `class3 = 31`
  - `class4 = 21`

Realized rates:

- baseline causal rate: `0.1053`
- enriched causal rate: `0.2600`
- baseline large-class rate (`class3` or `class4`): `0.0213`
- enriched large-class rate: `0.1440`
- baseline realized `P(delta > 2 | delta > 1)`: `0.2025`
- enriched realized `P(delta > 2 | delta > 1)`: `0.5538`
- baseline realized `P(delta > 3 | delta > 2)`: `0.1875`
- enriched realized `P(delta > 3 | delta > 2)`: `0.5000`

This is the first regime in this benchmark family where steps 2 and 3 were
given strong signal by design rather than only becoming less rare indirectly.

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

## Scenario 3 Results: `stepwise_annotation_signal`

Mean results across five seeds:

| Method | cor(y, EBV) | mean PIP causal | mean PIP noncausal | enriched PIP | baseline PIP | top-causal recall |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `BayesR` | `0.9164` | `0.4983` | `0.4840` | `0.4910` | `0.4844` | `0.2222` |
| `Annotated_BayesR` | `0.7409` | `0.1027` | `0.0606` | `0.1388` | `0.0426` | `0.2986` |
| `Annotated_BayesC` | `0.8595` | `0.3156` | `0.2340` | `0.4695` | `0.1712` | `0.2472` |

Main observations:

1. Annotated BayesR finally learned positive signal in all three steps.
   - step 1:
     - `Annotation_1 = +1.3349`
     - `Annotation_2 = -0.9746`
   - step 2:
     - `Annotation_1 = +0.2365`
     - `Annotation_2 = -0.2444`
   - step 3:
     - `Annotation_1 = +0.2518`
     - `Annotation_2 = +0.0275`
2. The deeper-step intercepts dropped materially relative to the earlier
   scenarios.
   - sparse regime:
     - step 2 intercept `19.2796`
     - step 3 intercept `28.1303`
   - less-sparse regime:
     - step 2 intercept `13.1159`
     - step 3 intercept `18.9376`
   - stepwise-signal regime:
     - step 2 intercept `11.0214`
     - step 3 intercept `8.3824`
3. Annotated BayesR again improved prioritization over BayesR.
   - causal minus noncausal mean PIP gap:
     - `BayesR = 0.0143`
     - `Annotated_BayesR = 0.0422`
   - enriched minus baseline mean PIP gap:
     - `BayesR = 0.0067`
     - `Annotated_BayesR = 0.0963`
   - top-causal recall:
     - `BayesR = 0.2222`
     - `Annotated_BayesR = 0.2986`
4. True-class ordering improved again.
   - BayesR mean PIP by true class:
     - `class1 = 0.4840`
     - `class2 = 0.4829`
     - `class3 = 0.4938`
     - `class4 = 0.5721`
   - Annotated BayesR mean PIP by true class:
     - `class1 = 0.0606`
     - `class2 = 0.0679`
     - `class3 = 0.0962`
     - `class4 = 0.2651`
5. Prediction still stayed below ordinary BayesR.
   - `BayesR = 0.9164`
   - `Annotated_BayesR = 0.7409`

Interpretation:

- This scenario answered the main open scientific question from the first two
  benchmarks: the annotated BayesR implementation can move in the correct
  direction on steps 2 and 3 when the truth explicitly drives those steps.
- The deeper sequential models are therefore not purely dead code or purely a
  numerical artifact.
- But the prediction penalty did not disappear, even after giving the deeper
  steps signal by design.
- The nuisance annotation still moved away from zero, especially in step 1 and
  step 2, even though it was generated independently of the truth. That is a
  real caution sign about overfitting in the annotation submodel.

## What This Says About Annotated BayesR

Across all three scenarios, the pattern was consistent.

### 1. The production implementation is working

The benchmark used only the production `runMCMC` path, and annotated BayesR
behaved coherently in all regimes:

- it learned the informative annotation at step 1
- it moved step-2 and step-3 informative coefficients in the expected direction
  once the truth actually supplied that signal
- it changed SNP prioritization in the expected direction

This is not a mechanics failure.

### 2. Annotated BayesR is better at prioritization than at prediction

In all three scenarios:

- enriched vs baseline PIP separation improved
- causal vs noncausal PIP separation improved
- top-causal recall improved
- class-order separation improved

But in all three scenarios:

- `cor(y, EBV)` stayed below ordinary BayesR

So the current implementation is stronger as a fine-mapping or prioritization
tool than as a predictor on these synthetic regimes.

### 3. The main difficulty is still the upper conditional models

Even after increasing class-3 and class-4 prevalence, the upper-step intercepts
were still large and the informative step-2/step-3 coefficients were weak.

The stepwise-signal scenario improved that picture:

- informative step-2 and step-3 coefficients became positive
- upper-step intercepts dropped substantially

But the remaining issue did not disappear:

- EBV correlation still lagged ordinary BayesR
- nuisance coefficients were still not cleanly suppressed

So the central scientific issue is no longer “can the deeper steps move at
all?” It is whether those deeper steps can improve posterior quality without
overfitting the annotation layer.

## What This Says About Annotated BayesC

Annotated BayesC remained a useful comparator, but it is not testing the same
scientific question as annotated BayesR in the stepwise-signal regime because
it only models step 1.

Across the three regimes:

- sparse regime:
  - strong enriched-vs-baseline separation: `0.3302`
  - reasonable EBV correlation: `0.8928`
- less-sparse regime:
  - enriched-vs-baseline gap turned negative: `-0.0372`
  - informative annotation coefficient mean became negative:
    - `Annotation_1 = -1.0116`
- stepwise-signal regime:
  - enriched-vs-baseline gap became positive again: `0.2983`
  - top-causal recall improved over BayesR: `0.2472` vs `0.2222`
  - `cor(y, EBV)` stayed below BayesR: `0.8595` vs `0.9164`

So annotated BayesC is still a useful comparator, but it is not a stable
scientific baseline across this benchmark family.

In the new stepwise-signal regime specifically:

- it achieved stronger enriched-vs-baseline separation than BayesR
- it improved top-causal recall over BayesR
- it stayed well below BayesR on `cor(y, EBV)`

So it is still a useful baseline for “does annotation help prioritization,” but
it is not a substitute for testing the deeper sequential logic of annotated
BayesR.

## Bottom Line

The three-scenario benchmark supports these conclusions.

1. Annotated BayesR is implemented correctly enough to learn the annotation and
   change posterior prioritization on the production JWAS path.
2. Annotated BayesR consistently improves marker prioritization relative to
   ordinary BayesR.
3. Annotated BayesR is not yet benchmark-proven to improve prediction
   (`cor(y, EBV)`) relative to ordinary BayesR.
4. The deeper sequential annotation models can respond correctly when the truth
   explicitly drives them, so the current limitation is not a basic
   implementation failure.
5. The remaining scientific weakness is the tradeoff between stronger
   prioritization and weaker prediction, together with signs of annotation-layer
   overfitting on nuisance features.
6. Annotated BayesC remains useful as a comparator, but it does not answer the
   same step-2/step-3 question as annotated BayesR.

So the current implementation should be viewed as:

- a sound v1 annotated BayesR sampler
- promising for SNP prioritization
- not yet validated as a superior predictor
- scientifically stronger after the stepwise-signal benchmark, but still not
  settled enough to claim the deeper annotation structure improves prediction

The next useful work is no longer another benchmark of the same kind. It is
code and model review:

- review the annotated BayesR sampler together
- focus on the annotation submodel regularization and nuisance-feature behavior
- decide whether v1 should be merged as a prioritization-oriented method or
  whether the annotation layer needs further restraint first
