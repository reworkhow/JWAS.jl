# Annotated BayesC Jian-Style Startup Report

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

- `/tmp/annotated_bayesc_seed_stability_jian_startup`

## Change Under Test

Annotated BayesC startup was changed to match Jian's startup structure more
closely:

- start from common marker-level `π_j`
- start annotation coefficients at zero
- skip pre-MCMC annotation fitting

## Cross-Seed Results

### Annotated BayesC

- marker correlation: `0.979755604995099`
- model-frequency correlation: `0.9963342868373964`
- EBV correlation: `0.9974783833696612`
- annotation-coefficient correlation: `0.9921264474320386`
- annotation-coefficient max absolute difference: `2.1879835564566097`

### Annotated BayesR

Unchanged observation on the same benchmark:

- marker correlation: `0.4532361462863509`
- model-frequency correlation: `-0.37499214409597487`
- EBV correlation: `0.7438258830716546`

### Ordinary BayesR

Unchanged control on the same benchmark:

- marker correlation: `0.9907735822811`
- model-frequency correlation: `0.74585582219425`
- EBV correlation: `0.9983051625434268`

## Comparison To The Earlier Startup Patch

The previous deterministic-intercept startup patch gave:

- annotated BayesC marker correlation: `0.6969476462530532`
- model-frequency correlation: `-0.9451556319728573`
- EBV correlation: `0.9371368410542033`
- annotation-coefficient correlation: `-0.9765307153717362`

After switching to the Jian-style zero-alpha startup:

- annotated BayesC marker correlation improved to `0.979755604995099`
- model-frequency correlation improved to `0.9963342868373964`
- EBV correlation improved to `0.9974783833696612`
- annotation-coefficient correlation improved to `0.9921264474320386`

## Conclusion

The earlier BayesC startup pre-fit was a real source of seed instability.

Matching Jian's startup structure largely fixes dense annotated BayesC seed
stability on this large reproducer. Annotated BayesR remains unstable on the same
benchmark, but that is now a separate problem.
