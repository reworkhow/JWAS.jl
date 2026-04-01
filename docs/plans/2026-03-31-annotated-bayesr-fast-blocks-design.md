# Annotated BayesR Fast Blocks Design

## Goal

Enable `fast_blocks` for annotated BayesR on the production single-trait dense
path and verify that the existing block sampler correctly consumes the
annotation-driven per-marker class priors.

## Motivation

Annotated BayesR was released as dense-only v1, but the current BayesR block
sampler already accepts per-marker prior rows from `annotations.snp_pi`. That
means the core block kernel is already annotation-aware; the remaining blockers
are validation, tests, and evidence that the production path runs correctly.

## Design

Keep the current BayesR block sampler unchanged. The feature work is:

- remove the validation guard that rejects `fast_blocks` when BayesR
  annotations are present
- replace the current rejection test with a real annotated BayesR block run
- add one production benchmark comparing dense vs block annotated BayesR on the
  existing synthetic benchmark harness
- update docs to state that annotated BayesR supports dense `fast_blocks`

## Scope

In scope:

- dense single-trait annotated BayesR with `fast_blocks`
- unit coverage for annotated BayesR block runs
- one production-path benchmark/report for dense vs block annotated BayesR
- manual and implementation notes

Out of scope:

- streaming annotated BayesR
- multi-trait annotated BayesR
- block-specific algorithm redesign
- changing the existing BayesR block repetition schedule
