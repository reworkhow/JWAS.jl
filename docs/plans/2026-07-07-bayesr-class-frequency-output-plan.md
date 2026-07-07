# BayesR Class Frequency Output Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add exact BayesR class-specific posterior marker frequencies so users can plot small-, medium-, large-, and medium-or-large-effect PIPs without changing existing BayesR `Model_Frequency` behavior. Also make the large marker-level MCMC sample files optional.

**Architecture:** Keep the current BayesR `Model_Frequency = mean(delta > 1)` output unchanged for backward compatibility. Add BayesR-only posterior accumulators for class frequencies `Pr(delta == 1)`, `Pr(delta == 2)`, `Pr(delta == 3)`, `Pr(delta == 4)`, and expose them as extra columns in the existing marker-effects output table. Do not save raw per-iteration delta class samples by default because that can become large for production marker counts and long chains. Add a user-facing `runMCMC` keyword to suppress marker-level MCMC sample files while keeping posterior summaries and smaller downstream files available.

**Tech Stack:** Julia, JWAS internal MCMC output code, `Test`, `DataFrames.jl`, `CSV.jl`, shell commands (`rg`, `julia --project=.`), git.

---

## Recommendation

Implement final class-frequency columns first, not raw delta sample files. Add an explicit marker-level sample switch at the same time so production BayesR runs can avoid the largest sample files.

This is the smallest production-friendly change that gives the exact plot we need:

- `Class1_Frequency = mean(delta == 1)`
- `Class2_Frequency = mean(delta == 2)`
- `Class3_Frequency = mean(delta == 3)`
- `Class4_Frequency = mean(delta == 4)`
- `MediumLarge_Frequency = mean(delta >= 3)`

It avoids large raw class sample files while preserving the existing output meaning:

- `Model_Frequency` remains `mean(delta > 1)`, the overall nonzero BayesR PIP.
- Users can plot medium/large PIP directly from `MediumLarge_Frequency`.
- Users can diagnose inflated overall PIP by comparing `Class2_Frequency` with `Class3_Frequency + Class4_Frequency`.

Raw delta sample output can be added later behind an explicit opt-in keyword if trace-level diagnostics become necessary.

For MCMC sample control, use a narrow keyword instead of changing `output_samples_frequency` semantics:

- Add `output_marker_effect_samples::Bool=true` to `runMCMC`.
- When `false`, do not create or write `MCMC_samples_marker_effects_*` files.
- Continue updating posterior means, variances, `Model_Frequency`, and the new BayesR class-frequency accumulators.
- Continue writing smaller sample files that existing downstream summaries may depend on, such as variance, pi, EBV, and heritability outputs.
- Reject `output_marker_effect_samples=false` with `causal_structure`, because SEM post-processing requires marker-effect MCMC sample files.

Avoid using `output_samples_frequency=0` as the public off switch. The current validation rejects non-positive frequencies, and the current MCMC loop uses this frequency to decide when posterior summaries are accumulated. Reusing it as both "sampling interval" and "do not write files" would be more error-prone than a boolean output switch.

## Task 1: Add Optional Marker-Level MCMC Sample Output

**Files:**
- Modify: `src/1.JWAS/src/JWAS.jl`
- Modify: `src/1.JWAS/src/types.jl`
- Modify: `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`
- Modify: `src/1.JWAS/src/output.jl`
- Modify: `test/unit/test_bayesr.jl`
- Modify: relevant docs page where `output_samples_frequency` is documented

**Step 1: Add the public keyword**

Add a keyword to both `runMCMC` methods:

```julia
output_marker_effect_samples::Bool = true
```

Store it in `MCMCinfo` so lower-level MCMC and output code can branch on the setting.

**Step 2: Keep posterior accumulation independent from marker sample file output**

In the MCMC loop, continue calling `output_posterior_mean_variance` at the existing `output_samples_frequency` cadence. Only the marker-effect sample file setup and writes should be skipped when `output_marker_effect_samples=false`.

The intended behavior is:

```julia
if output_samples_frequency != 0
    output_posterior_mean_variance(...)
    output_MCMC_samples(...; output_marker_effect_samples=...)
end
```

Inside `output_MCMC_samples_setup` and `output_MCMC_samples`, gate only these file families:

- `MCMC_samples_marker_effects_*`
- any future raw marker class sample files

Do not gate smaller files that other JWAS output summaries may rely on, including pi, marker-effect variance, residual variance, EBV, and heritability sample files.

Add a clear validation error if `output_marker_effect_samples=false` is combined with `causal_structure`, because SEM post-processing uses the direct marker-effect sample files to generate indirect and overall marker-effect samples.

**Step 3: Add a focused regression test**

Extend the BayesR output test to call:

```julia
runMCMC(...;
        output_marker_effect_samples=false,
        output_samples_frequency=5,
        outputEBV=false,
        output_heritability=false)
```

Assert:

- final marker-effects output is still returned
- final marker-effects output still has `Model_Frequency`
- no `MCMC_samples_marker_effects_*` file exists in the output folder
- `MCMC_samples_pi_*` still exists for BayesR

**Step 4: Update documentation**

Document the distinction:

- `output_samples_frequency` controls how often retained posterior samples are accumulated/written.
- `output_marker_effect_samples=false` suppresses the large marker-effect sample files.
- Final marker-effect summaries remain available.

## Task 2: Add A BayesR Class-Frequency Regression Test

**Files:**
- Modify: `test/unit/test_bayesr.jl`

**Step 1: Add a focused output test**

Extend the existing BayesR output testset or add a new one:

```julia
@testset "BayesR class frequency output" begin
    global geno = get_genotypes(genofile, 1.0, separator=',',
                                method="BayesR",
                                Pi=Float64[0.94, 0.049, 0.01, 0.001],
                                estimatePi=true,
                                estimate_variance=true)
    model = build_model("y1 = intercept + geno", 1.0)
    outdir = tempname()
    output = runMCMC(model, phenotypes,
                     chain_length=30,
                     burnin=5,
                     output_samples_frequency=5,
                     output_folder=outdir,
                     seed=321,
                     printout_model_info=false,
                     outputEBV=false,
                     output_heritability=false,
                     fast_blocks=false)

    marker_effects = output["marker effects geno"]
    for col in [:Class1_Frequency, :Class2_Frequency, :Class3_Frequency,
                :Class4_Frequency, :MediumLarge_Frequency]
        @test col in propertynames(marker_effects)
        @test all(0.0 .<= marker_effects[!, col] .<= 1.0)
    end

    @test all(isapprox.(
        marker_effects[!, :Model_Frequency],
        marker_effects[!, :Class2_Frequency] .+
        marker_effects[!, :Class3_Frequency] .+
        marker_effects[!, :Class4_Frequency];
        atol=1e-12,
    ))
    @test all(isapprox.(
        marker_effects[!, :MediumLarge_Frequency],
        marker_effects[!, :Class3_Frequency] .+
        marker_effects[!, :Class4_Frequency];
        atol=1e-12,
    ))

    isdir(outdir) && rm(outdir, recursive=true)
end
```

**Step 2: Run the test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected: FAIL because the marker-effects table does not yet include the class-frequency columns.

## Task 3: Add BayesR Class-Frequency State

**Files:**
- Modify: `src/1.JWAS/src/types.jl`
- Modify: `src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl`

**Step 1: Add a field to `Genotypes`**

Add a field near `meanDelta` in `Genotypes`:

```julia
meanDeltaClass
```

Update the inner constructor default so non-BayesR analyses initialize it to `false`.

**Step 2: Initialize BayesR class-frequency arrays**

In `MCMC_BayesianAlphabet`, after the existing BayesR delta initialization:

```julia
if Mi.method == "BayesR"
    Mi.delta = [ones(Int, Mi.nMarkers) for traiti = 1:Mi.ntraits]
    Mi.meanDelta = [zeros(Float64, Mi.nMarkers) for traiti = 1:Mi.ntraits]
    Mi.meanDeltaClass = [
        zeros(Float64, Mi.nMarkers, length(BAYESR_GAMMA))
        for traiti = 1:Mi.ntraits
    ]
end
```

Use the repository's actual field name `Mi.δ` in code, not `Mi.delta`; the snippet above is written with ASCII in prose for clarity.

**Step 3: Run the focused test**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected: still FAIL because the accumulators are not updated or output yet.

## Task 4: Accumulate BayesR Class Frequencies During Posterior Sampling

**Files:**
- Modify: `src/1.JWAS/src/output.jl`

**Step 1: Update posterior accumulation**

In `output_posterior_mean_variance`, keep the current BayesR `meanDelta` update:

```julia
Mi.meanDelta[trait] += (Float64.(Mi.δ[trait] .> 1) - Mi.meanDelta[trait]) / nsamples
```

Add class-frequency updates immediately after it:

```julia
for klass in 1:length(BAYESR_GAMMA)
    Mi.meanDeltaClass[trait][:, klass] +=
        (Float64.(Mi.δ[trait] .== klass) - Mi.meanDeltaClass[trait][:, klass]) / nsamples
end
```

**Step 2: Run the focused test**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected: still FAIL because the output table does not expose the new accumulators.

## Task 5: Add Class-Frequency Columns To Marker Effects Output

**Files:**
- Modify: `src/1.JWAS/src/output.jl`

**Step 1: Append BayesR-only columns in `output_result`**

After constructing:

```julia
output["marker effects " * Mi.name] = DataFrame(...)
```

add BayesR-specific columns. Build one vector per output row, matching the trait-by-trait row order already used for `whichtrait`, `whichmarker`, and `whichdelta`:

```julia
if Mi.method == "BayesR"
    class1 = Mi.meanDeltaClass[1][:, 1]
    class2 = Mi.meanDeltaClass[1][:, 2]
    class3 = Mi.meanDeltaClass[1][:, 3]
    class4 = Mi.meanDeltaClass[1][:, 4]
    for traiti in 2:ntraits_geno
        class1 = vcat(class1, Mi.meanDeltaClass[traiti][:, 1])
        class2 = vcat(class2, Mi.meanDeltaClass[traiti][:, 2])
        class3 = vcat(class3, Mi.meanDeltaClass[traiti][:, 3])
        class4 = vcat(class4, Mi.meanDeltaClass[traiti][:, 4])
    end

    marker_output = output["marker effects " * Mi.name]
    marker_output[!, :Class1_Frequency] = class1
    marker_output[!, :Class2_Frequency] = class2
    marker_output[!, :Class3_Frequency] = class3
    marker_output[!, :Class4_Frequency] = class4
    marker_output[!, :MediumLarge_Frequency] = class3 .+ class4
end
```

**Step 2: Run the focused BayesR test**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected: PASS for the new class-frequency checks.

## Task 6: Update Documentation

**Files:**
- Modify: `docs/src/manual/bayesc_bayesr_comparison.md`

**Step 1: Document BayesR output meanings**

Add a short section near the BayesR prior/output discussion:

```markdown
## BayesR Marker Frequency Output

For BayesR, `Model_Frequency` reports the posterior nonzero frequency
`Pr(delta > 1)`, so it includes small-, medium-, and large-effect classes.

BayesR marker-effect output also reports class-specific frequencies:

- `Class1_Frequency = Pr(delta == 1)`, the zero-effect class
- `Class2_Frequency = Pr(delta == 2)`, the small-effect class
- `Class3_Frequency = Pr(delta == 3)`, the medium-effect class
- `Class4_Frequency = Pr(delta == 4)`, the large-effect class
- `MediumLarge_Frequency = Pr(delta >= 3)`

Use `MediumLarge_Frequency` when the plot should ignore small class-2 effects.
```

**Step 2: Document optional marker-effect sample output**

Near the existing `runMCMC` output options, document:

```markdown
Set `output_marker_effect_samples=false` to skip writing the large
`MCMC_samples_marker_effects_*` files. Final marker-effect summaries,
including `Estimate`, `SD`, `Model_Frequency`, and BayesR class-frequency
columns, are still computed.
```

Keep the documentation clear that `output_samples_frequency` controls the retained posterior sampling cadence; it is not the on/off switch for marker-effect sample files.

**Step 3: Run docs if feasible**

Run:

```bash
julia --project=docs --startup-file=no docs/make.jl
```

Expected: PASS. If this is too slow for the current environment, record that docs were not run and why.

## Task 7: Run Verification

**Files:**
- No edits.

**Step 1: Run focused tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_bayesr.jl")'
```

Expected: PASS.

**Step 2: Run relevant output tests**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_annotated_bayesr.jl")'
```

Expected: PASS. This checks that annotated BayesR output still works.

**Step 3: Run full test suite if time permits**

Run:

```bash
julia --project=. --startup-file=no test/runtests.jl
```

Expected: PASS.

## Task 8: Write Implementation Record

**Files:**
- Create: `docs/plans/2026-07-07-bayesr-class-frequency-output-implementation.md`

**Step 1: Record what changed**

Include:

- why `Model_Frequency` remains unchanged
- the new BayesR class-frequency columns
- the new `output_marker_effect_samples=false` option for suppressing large marker-level sample files
- the test commands run
- any docs build result
- any known limitation, especially that raw per-iteration delta class samples are still not saved

**Step 2: Commit**

```bash
git add \
  src/1.JWAS/src/JWAS.jl \
  src/1.JWAS/src/types.jl \
  src/1.JWAS/src/MCMC/MCMC_BayesianAlphabet.jl \
  src/1.JWAS/src/output.jl \
  test/unit/test_bayesr.jl \
  docs/src/manual/bayesc_bayesr_comparison.md \
  docs/plans/2026-07-07-bayesr-class-frequency-output-plan.md \
  docs/plans/2026-07-07-bayesr-class-frequency-output-implementation.md
git commit -m "feat: report BayesR class frequencies"
```
