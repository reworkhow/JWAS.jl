# Annotated BayesR Jian-Style Startup Design

## Goal

Align dense annotated BayesR startup with Jian's `sbayesrc.R` startup structure:

- start from the supplied common class prior
- start annotation coefficients at zero
- do not fit the annotation model before the first phenotype-informed marker sweep

The motivation is the same as for annotated BayesC: remove a startup-only source
of seed dependence and make JWAS match the intended sampler order more closely.

## Design

### 1. Startup prior

Annotated BayesR should keep the supplied starting mixture prior as a constant
marker-level `snp_pi` matrix:

`snp_pi[j, :] = Pi` for all markers `j`.

This is already the object consumed by the BayesR marker sampler.

### 2. Startup annotation state

The three annotation-step coefficient vectors should all start at zero.

This means the annotation linear predictors `mu` also start at zero. The startup
annotation state is therefore just a placeholder; it does not need to reproduce
the supplied `snp_pi` yet.

### 3. MCMC startup flow

JWAS should not run the annotated BayesR annotation update before the first marker
sweep. The first BayesR marker update should use the supplied `snp_pi` directly,
and the first annotation fit should happen only inside the normal `estimatePi`
section after phenotype-informed class assignments have been sampled once.

### 4. Scope

This pass is narrow:

- dense annotated BayesR startup only
- no fast-block algorithm changes
- no annotation regularization redesign
- no BayesC changes beyond preserving the already-landed Jian-style startup there

## Expected Outcome

If startup pre-fit was contributing materially to BayesR instability, the large
two-seed reproducer should show noticeably improved cross-seed agreement after
this change.
