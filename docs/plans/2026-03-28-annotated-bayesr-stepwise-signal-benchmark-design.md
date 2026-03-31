# Annotated BayesR Stepwise-Signal Benchmark Design

## Goal

Add a third production benchmark scenario for annotated BayesR where the
annotation truth is specified directly in the sequential BayesR conditional
probabilities:

- `p1 = P(delta > 1)`
- `p2 = P(delta > 2 | delta > 1)`
- `p3 = P(delta > 3 | delta > 2)`

The purpose is to test whether the annotated BayesR implementation can recover
signal in steps 2 and 3, not just step 1.

## Motivation

The current two scenarios established:

- annotated BayesR works mechanically
- annotated BayesR improves prioritization relative to ordinary BayesR
- step-2 and step-3 annotation effects remain weak

But both existing scenarios were still driven mostly by step-1 signal. Making
upper classes less rare helped only indirectly. It did not explicitly create
large annotation differences in the deeper sequential conditionals.

## Proposed Scenario

New scenario name:

- `stepwise_annotation_signal`

Truth is defined through conditional probabilities.

Baseline SNPs:

- `p1_base = 0.10`
- `p2_base = 0.20`
- `p3_base = 0.20`

Enriched SNPs:

- `p1_enriched = 0.30`
- `p2_enriched = 0.60`
- `p3_enriched = 0.60`

Convert to 4-class joint probabilities:

- `pi1 = 1 - p1`
- `pi2 = p1 * (1 - p2)`
- `pi3 = p1 * p2 * (1 - p3)`
- `pi4 = p1 * p2 * p3`

That gives:

Baseline:

- `pi = [0.90, 0.08, 0.016, 0.004]`

Enriched:

- `pi = [0.70, 0.12, 0.072, 0.108]`

This benchmark deliberately creates annotation signal at every step:

- step 1: `0.10` vs `0.30`
- step 2: `0.20` vs `0.60`
- step 3: `0.20` vs `0.60`

## Scope

Keep the existing benchmark structure:

- production `runMCMC` path only
- methods:
  - `BayesR`
  - `Annotated_BayesR`
  - `Annotated_BayesC`
- single trait
- dense only
- `n_obs = 200`
- `n_markers = 1000`
- `chain_length = 10000`
- `burnin = 2000`
- seeds `2026:2030`

## Expected Readout

If annotated BayesR is using the deeper annotation hierarchy well, this scenario
should show:

- positive `Annotation_1` coefficients in all three steps
- step-2 and step-3 coefficients that are materially larger and more stable than
  in the first two scenarios
- stronger true-class ordering than ordinary BayesR
- potentially a smaller EBV penalty, or at least a better explanation for where
  the remaining penalty comes from

## Outputs

No new output files are needed. Extend the existing benchmark outputs with the
new scenario label:

- `comparison_runs.csv`
- `comparison_summary.csv`
- `pip_group_summary.csv`
- `annotation_coefficients.csv`
- `truth_metadata.csv`

Update the existing benchmark report to compare all three scenarios.

## Risks

- Annotated BayesC may remain unstable or uninformative in this regime
- Even with explicit stepwise truth, prediction may still favor ordinary BayesR
- If step-2 and step-3 remain weak here, the issue is more likely structural in
  the model or the available sample size rather than just poor benchmark design
