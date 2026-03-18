# BayesR Design

**Date:** 2026-03-17

## Goal

Add a first BayesR implementation to JWAS for single-trait dense genotype analysis, with code structure aligned to the existing BayesABC path and with clear extension points for later `stream`, `fast_blocks`, and annotation support.

## Scope

### In scope
- Add `method="BayesR"` for single-trait dense genotype analyses.
- Reuse the current Bayesian Alphabet dense data flow and MCMC integration pattern.
- Implement a BayesR-specific Gibbs sampler with fixed mixture multipliers `gamma = [0.0, 0.01, 0.1, 1.0]`.
- Add concise setup printouts that help users verify priors and initialization.
- Add BayesR-specific validation, initialization, posterior bookkeeping, and output handling.
- Add focused unit and integration tests for the new BayesR path.

### Out of scope
- Multi-trait BayesR.
- `storage=:stream` support.
- `fast_blocks` support.
- Annotation-aware BayesR.
- RRM support.
- BayesABC wrapper refactoring or sampler interface cleanup unrelated to BayesR.
- User-configurable BayesR mixture classes in v1.

## Architecture

BayesR should be integrated as a new method in the existing Bayesian Alphabet pipeline rather than as a separate MCMC driver. The implementation should match BayesABC structurally wherever that does not obscure the statistical differences between the models.

Planned structure:
- Add BayesR method validation in `src/1.JWAS/src/input_data_validation.jl`.
- Add a BayesR dense sampler file at `src/1.JWAS/src/markers/BayesianAlphabet/BayesR.jl`.
- Include the new file from `src/1.JWAS/src/JWAS.jl`.
- Add a BayesR branch in `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`.
- Keep BayesR state inside the existing `Genotypes` object:
  - `Mi.G.val` stores the shared BayesR base variance `sigmaSq`.
  - `Mi.π` stores the global 4-element mixture probability vector.
  - `Mi.δ[1]` stores per-marker class labels `1:4`.
  - `Mi.α[1]` stores the sampled marker effects.

## Sampler Design

BayesR will follow the dense single-trait BayesABC calling pattern:
- High-level wrapper: `BayesR!(genotypes, ycorr, vare)`
- Dense kernel: `BayesR!(xArray, xRinvArray, xpRinvx, yCorr, α, δ, vare, sigmaSq, π, gamma)`

The Gibbs sweep for each marker will:
1. Reconstruct the right-hand side using the current effect, following the same in-place residual update pattern used by BayesABC.
2. Compute class log posterior values for:
   - class 1: exact zero effect
   - classes 2-4: normal components with variances `gamma[k] * sigmaSq`
3. Normalize class probabilities using a stable log-sum-exp calculation.
4. Sample the class label `δ[j]`.
5. Sample `α[j]` from the class-specific normal posterior for nonzero classes, or set it to zero for the zero class.
6. Update `yCorr` in place with the `oldAlpha - newAlpha` BLAS pattern.

The implementation should remain visually and structurally close to `BayesABC.jl`, but BayesR-specific helpers should keep the categorical mixture logic explicit instead of forcing it into binary-inclusion abstractions.

## Priors And Initialization

BayesR will use the fixed mixture multipliers:
- `gamma = [0.0, 0.01, 0.1, 1.0]`

Initialization rules:
- `Pi` for BayesR must be a 4-element vector summing to 1, or `0.0` to request the BayesR default initialization path.
- `Mi.G.val` represents the shared BayesR base variance `sigmaSq`, not a per-marker variance vector.
- When JWAS derives marker-effect variance from genetic variance, BayesR should use:
  - `sigmaSq = genetic_variance / (sum2pq * sum(gamma .* pi))`
- `Mi.δ[1]` should be initialized by sampling class labels from `Mi.π`, with a small safeguard against a degenerate all-zero start.
- `Mi.α[1]` can remain zero-initialized before the first Gibbs sweep.

Useful setup information should be printed once, concisely, during model setup. The printout should include only short, high-value items such as:
- method name
- mixture multipliers
- starting `pi`
- implied shared marker variance prior mean
- initial class counts

## Validation And Restrictions

BayesR should fail early and clearly for unsupported configurations:
- multi-trait models
- `storage=:stream`
- `fast_blocks != false`
- annotation inputs
- RRM

BayesR-specific input validation should also reject:
- `Pi` values that are not length 4
- `Pi` values with negative entries
- `Pi` values that do not sum to 1 within a small tolerance
- derived priors with nonpositive `sum(gamma .* pi)` or nonpositive `sigmaSq`

Numerical behavior should be more robust than the tutorial R example:
- handle the zero-effect class explicitly instead of relying on `log(0)`
- use stable log-probability normalization
- throw clear errors for impossible states rather than silently masking them

## Posterior Updates And Output

After each BayesR sweep:
- if `estimatePi=true`, update `Mi.π` by Dirichlet sampling from the 4 class counts
- if `estimate_variance=true`, update `Mi.G.val` as the shared base variance `sigmaSq`

The shared variance update should reuse existing JWAS scalar variance machinery by transforming nonzero effects according to their current mixture class.

Output requirements:
- BayesR `pi` output must preserve all 4 mixture probabilities
- single-trait BayesR must not collapse the `pi` vector to a scalar
- the existing marker-effects table may stay in place
- `Model_Frequency` in the marker-effects output should mean posterior nonzero frequency `Pr(δ > 1)`
- `Mi.meanVara` should report posterior mean of `sigmaSq`

## Testing Strategy

Testing should stay narrow and implementation-focused for v1.

Planned tests:
- validation failures for unsupported BayesR configurations
- validation failures for malformed `Pi`
- a small dense single-trait integration test using `method="BayesR"`
- output checks confirming BayesR returns a 4-row `pi` result
- sampler checks on a tiny synthetic problem confirming:
  - class labels stay in `1:4`
  - mixture probabilities stay valid
  - BayesR marker effects output is produced

Likely first test file:
- `test/unit/test_bayesr.jl`

## Future Compatibility

The first BayesR slice should stay small, but the code should leave a clean path for later extensions:
- `storage=:stream`
- `fast_blocks`
- annotation-aware BayesR

That means v1 should keep method-specific logic localized, avoid hardwiring binary-inclusion assumptions into shared helpers, and preserve clear separation between:
- dense BayesR marker updates
- BayesR-specific validation
- BayesR-specific output behavior

## Risks And Mitigations

- Risk: BayesR code drifts too far from BayesABC structure.
  - Mitigation: mirror BayesABC wrapper/kernel organization and in-place residual update style.
- Risk: BayesR output is silently forced into BayesC single-scalar `pi` assumptions.
  - Mitigation: add method-aware output handling and direct tests for 4-class `pi`.
- Risk: later stream/block/annotation support becomes awkward.
  - Mitigation: keep v1 restrictions explicit and isolate BayesR-specific mixture logic instead of spreading it across shared code paths.
