# Annotated BayesR Jian-Style Startup Report

## Benchmark

Dense production-path benchmark:

- `n_obs = 200`
- `n_markers = 1000`
- `chain_length = 10000`
- `burnin = 2000`
- data seed `20260402`
- MCMC seeds `2026, 2027`
- scenario `less_sparse_upper_classes`

Output:

- `/tmp/annotated_bayesr_seed_stability_jian_startup`

## Change Under Test

Annotated BayesR startup was changed to match Jian's startup structure more
closely:

- start from constant marker-level `snp_pi`
- start annotation coefficients at zero
- skip pre-MCMC annotation fitting

## Cross-Seed Results

### Before This Change

Using the earlier BayesC-only Jian-style startup output:

- marker correlation: `0.4532361462863509`
- model-frequency correlation: `-0.37499214409597487`
- EBV correlation: `0.7438258830716546`
- `pi` correlation: `0.301881351255632`
- annotation-coefficient correlation: `0.13284637371772`
- annotation-coefficient max absolute difference: `25.371907060952722`

### After This Change

- marker correlation: `0.9050738036570479`
- model-frequency correlation: `0.289040252413319`
- EBV correlation: `0.984134868900681`
- `pi` correlation: `0.941265550119911`
- annotation-coefficient correlation: `0.21841901001773492`
- annotation-coefficient max absolute difference: `43.03721540860477`

### Controls

Annotated BayesC on the same benchmark remains:

- marker correlation: `0.979755604995099`
- model-frequency correlation: `0.9963342868373964`
- EBV correlation: `0.9974783833696612`

Ordinary BayesR on the same benchmark remains:

- marker correlation: `0.9907735822811`
- model-frequency correlation: `0.74585582219425`
- EBV correlation: `0.9983051625434268`

## Conclusion

The Jian-style startup fix helps annotated BayesR substantially:

- marker agreement improved a lot
- EBV agreement improved a lot
- `pi` agreement improved a lot

But annotated BayesR is still not as stable as annotated BayesC or ordinary
BayesR. Startup mismatch was a real part of the problem, but not the whole
problem. The remaining instability is now in the annotation sampler behavior
after startup.
