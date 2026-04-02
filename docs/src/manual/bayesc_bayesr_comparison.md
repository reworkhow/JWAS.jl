# BayesC and BayesR Comparison (Single-Trait)

This page compares the current **single-trait** JWAS implementations of:

- `BayesC`
- `BayesR`
- `Annotated_BayesC`
- `Annotated_BayesR`

The comparison is organized around:

1. starting values
2. prior structure
3. sampler order

## Starting Values

| Method | Mixture state at iteration 1 | Marker prior at iteration 1 | Annotation state at iteration 1 |
| :-- | :-- | :-- | :-- |
| `BayesC` | Marker effects start at `0`. Indicators are binary, `δ_j ∈ {0,1}`. | One shared exclusion probability `π`. In typical JWAS BayesC runs this starts from the supplied scalar `Pi`. | None. |
| `BayesR` | Marker effects start at `0`. Indicators are four-class, `δ_j ∈ {1,2,3,4}`. The shared BayesR variance factor `sigmaSq` starts either from the supplied marker-effect variance input, or from a genetic-variance-to-`sigmaSq` conversion if the user supplied genetic variance instead. | One shared four-class probability vector `π = (π_1, π_2, π_3, π_4)`. If `Pi=0.0`, JWAS uses the default BayesR start `(0.95, 0.03, 0.015, 0.005)`. | None. |
| `Annotated_BayesC` | Same binary marker state as BayesC. Marker effects start at `0`. | JWAS expands the supplied BayesC starting `Pi` to a marker-level exclusion vector `π_j`, so all markers start from the same BayesC prior. | Annotation coefficients start at `0`. Annotation linear predictor `μ` starts at `0`. JWAS does not fit the annotation model before the first phenotype-informed marker sweep. |
| `Annotated_BayesR` | Same four-class marker state as BayesR. Marker effects start at `0`. | JWAS expands the supplied BayesR starting vector `π` to a constant marker-level class-probability matrix `snp_pi`, so all markers start from the same four-class BayesR prior. | Step-specific annotation coefficients start at `0`. Annotation linear predictors `μ` start at `0`. JWAS does not fit the annotation model before the first phenotype-informed marker sweep. |

## Starting Value Sources

This section states whether a starting value is actually used, and if so whether
it is:

- supplied by the user
- filled from a JWAS default
- computed from other supplied quantities
- not used at all

| State | `BayesC` | `BayesR` | `Annotated_BayesC` | `Annotated_BayesR` |
| :-- | :-- | :-- | :-- | :-- |
| Marker effects `α` | If `starting_value` is provided to `get_genotypes`, JWAS uses it. Otherwise marker effects start at `0`. | If `starting_value` is provided to `get_genotypes`, JWAS uses it. Otherwise marker effects start at `0`. | Same rule as `BayesC`: user-supplied `starting_value` if provided, otherwise `0`. | Same rule as `BayesR`: user-supplied `starting_value` if provided, otherwise `0`. |
| Marker indicators `δ` | No separate user starting value is used. The first sampled `δ` values come from the first BayesC marker sweep. | No separate user starting value is used. The first sampled class labels come from the first BayesR marker sweep. | No separate user starting value is used. The first sampled `δ` values come from the first phenotype-informed annotated BayesC marker sweep. | No separate user starting value is used. The first sampled class labels come from the first phenotype-informed annotated BayesR marker sweep. |
| Mixture prior `π` | JWAS uses the supplied scalar `Pi`. If no `Pi` is supplied, the default argument is `Pi=0.0`. | JWAS uses the supplied 4-class `Pi`. If `Pi=0.0`, JWAS replaces it with the BayesR default `(0.95, 0.03, 0.015, 0.005)`. | JWAS uses the supplied BayesC `Pi`, then expands it to a marker-level vector `π_j`. If no `Pi` is supplied, the default argument is `Pi=0.0`. | JWAS uses the supplied BayesR `Pi`, then expands it to a constant marker-level matrix `snp_pi`. If `Pi=0.0`, JWAS first replaces it with the BayesR default `(0.95, 0.03, 0.015, 0.005)`. |
| Annotation coefficients | Not used. | Not used. | Used. If no user start is provided, they start at `0`. | Used. If no user start is provided, they start at `0`. |
| Annotation linear predictor `μ` | Not used. | Not used. | Used. If no user start is provided, it starts at `0`. | Used. If no user start is provided, it starts at `0`. |
| Marker-effect variance start | A positive starting value is used. It comes from the supplied genomic variance input `G`, or from a value computed from phenotypes if `G` is omitted. | A positive starting value is used. Two cases exist. If `G_is_marker_variance=true`, JWAS uses the supplied value directly as the BayesR shared marker-effect scale `sigmaSq`. If `G_is_marker_variance=false`, JWAS treats `G` as genetic variance and converts it to `sigmaSq`. | Same as `BayesC`. | Same as `BayesR`. |
| Marker-effect variance prior scale | Used. Once the starting marker-effect variance mean is available, JWAS computes `G.scale` from it. | Used. Once the starting `sigmaSq` mean is available, JWAS computes `G.scale` from it. | Same as `BayesC`. | Same as `BayesR`. |
| Residual variance start | A positive starting value is used. It comes from the `R` argument in `build_model`, or from a value computed from phenotypes if `R` is omitted. | Same as `BayesC`. | Same as `BayesC`. | Same as `BayesC`. |
| Residual variance prior scale | Used. Once the starting residual variance mean is available, JWAS computes `R.scale` from it. | Same as `BayesC`. | Same as `BayesC`. | Same as `BayesC`. |

For the BayesR genetic-variance conversion, JWAS uses the production rule in
`genetic2marker`:

```math
\sigma^2 = \frac{V_g}{\mathrm{sum2pq} \cdot \sum_k \gamma_k \pi_k}
```

where:

- `V_g` is the **genetic variance**, not the marker-effect variance
- `sum2pq = \sum_j 2 p_j (1-p_j)` for the current marker set
- `gamma_k` are the BayesR class variance multipliers
- `pi_k` are the current BayesR class probabilities

If the user instead supplies marker-effect variance directly by setting
`G_is_marker_variance=true`, JWAS skips this conversion and uses that supplied
value as `sigmaSq`.

## Variance-Scale Starting Values

JWAS initializes the variance-component scale parameters from the corresponding
starting variance means.

For single-trait residual variance, JWAS uses:

```math
\mathrm{scale}_R = R \cdot \frac{df - 2}{df}
```

where `R` is the starting residual variance mean.

For single-trait marker-effect variance, JWAS uses:

```math
\mathrm{scale}_G = G \cdot \frac{df - 2}{df}
```

where `G` is the starting marker-effect variance mean.

For BayesR, this same marker-scale formula is applied after the starting shared
marker-effect variance `sigmaSq` has been determined. So if JWAS is given
genetic variance first, the order is:

1. convert genetic variance `V_g` to BayesR `sigmaSq`
2. compute `G.scale` from that `sigmaSq`

In that step, `sigmaSq` is the shared **marker-effect variance** scale, while
`V_g` is the **genetic variance**. They are not the same quantity.

## Prior Structure

| Method | Marker-effect prior | Prior on mixture state | Annotation role |
| :-- | :-- | :-- | :-- |
| `BayesC` | Two-component mixture: one zero-effect class and one nonzero normal class. The nonzero effect variance is the current marker-effect variance `G`. | One shared exclusion probability `π`, updated from the total number of included markers. | None. |
| `BayesR` | Four-component mixture: class 1 is zero; classes 2 to 4 are normal with variances `gamma[k] * sigmaSq`. | One shared four-class vector `π`, updated from class counts with a Dirichlet draw. | None. |
| `Annotated_BayesC` | Same BayesC effect prior as above. | Marker-specific exclusion probabilities `π_j`. JWAS uses the probit rule `Pr(δ_j = 1 | a_j, γ) = Φ(a_j'γ)`, so the stored BayesC exclusion prior is `π_j = 1 - Φ(a_j'γ)`. | One binary annotation model with latent error variance fixed to `1`. |
| `Annotated_BayesR` | Same BayesR effect prior as above. | Marker-specific four-class prior matrix `π_j`. JWAS builds it from three conditional probabilities: `p1_j = Pr(δ_j > 1)`, `p2_j = Pr(δ_j > 2 | δ_j > 1)`, and `p3_j = Pr(δ_j > 3 | δ_j > 2)`. | Three nested annotation models, each with latent error variance fixed to `1`. |

## Sampler Order

| Method | Marker update | Mixture-prior update | Variance updates |
| :-- | :-- | :-- | :-- |
| `BayesC` | Sample binary inclusion `δ_j` and marker effect for each marker. | Sample one shared exclusion probability `π` from the current number of included markers. | Update marker-effect variance `G`, then update residual variance. |
| `BayesR` | Sample one four-way class assignment `δ_j` and then the corresponding marker effect for each marker. | Sample one shared four-class vector `π` from the current class counts. | Update shared BayesR `sigmaSq`, then update residual variance. |
| `Annotated_BayesC` | Run the ordinary BayesC marker sweep using the current marker-level `π_j`. | After the marker sweep, sample annotation liabilities and annotation coefficients, then rewrite the full marker-level exclusion vector `π_j`. | Update BayesC marker-effect variance `G`, then update residual variance. |
| `Annotated_BayesR` | Run the ordinary BayesR marker sweep using the current marker-level four-class prior matrix `snp_pi`. | After the marker sweep, construct `z1 = 1(δ > 1)`, `z2 = 1(δ > 2)`, `z3 = 1(δ > 3)`, fit the three annotation submodels, then rebuild the full marker-level four-class prior matrix `snp_pi`. | Update shared BayesR `sigmaSq`, then update residual variance. |

## Key Distinctions

| Comparison | Main difference |
| :-- | :-- |
| `BayesC` vs `BayesR` | `BayesC` is a binary inclusion model with one shared exclusion probability. `BayesR` is a four-class mixture model with one shared class-probability vector. |
| `Annotated_BayesC` vs `Annotated_BayesR` | `Annotated_BayesC` uses one annotation-driven binary probit model to produce marker-specific `π_j`. `Annotated_BayesR` uses three nested annotation probit models to produce marker-specific four-class probabilities `snp_pi`. |
| Plain vs annotated methods | The plain methods learn one shared mixture prior. The annotated methods learn marker-specific priors from the supplied annotations. |

## Related Pages

- [Annotated BayesC](annotated_bayesc.md)
- [Annotated BayesR](annotated_bayesr.md)
- [Block BayesC](block_bayesc.md)
