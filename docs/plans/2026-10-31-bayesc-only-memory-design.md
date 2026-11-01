# BayesC-only memory and residual optimization

Base: latest fetched reworkhow master 038c6c8 (JWAS 2.4.0).
Scope: BayesC ordinary and fast-block, single- and multi-trait; include both
existing MT samplers and sequential/independent blocks. No BayesR sampler edits.
BayesC defaults to no-copy matrix input and contiguous shared genotype views.
Other methods retain upstream defaults and block storage.

Use one in-place block residual path, with BLAS multipliers matching the
genotype block element type. Remove the residual-mode keyword and MT product
workspaces. Keep block-sized alpha copies. Rounding can change seeded chains;
compare predictions against upstream and report any differences explicitly.

Preserve priors, annotation models, sampling schedules and prediction formulas.
Support one dense BayesC term with matching precision, contiguous unique rows
and unit weights; keep opt-outs for unsupported layouts/modes. File/DataFrame
ownership stays unchanged. Validate default and copied storage, partial trait
missingness, held-out predictions, sampler I/II, both precisions and block modes.
Compare keyed predictions directly with upstream across multiple seeds; test
the unified residual path. Run full regressions and documentation. Keep master and
upstream unchanged until review.
