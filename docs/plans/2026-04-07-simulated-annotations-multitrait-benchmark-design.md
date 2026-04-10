# Simulated Annotations Multi-Trait Benchmark Design

## Goal

Extend the packaged `simulated_annotations` fixture with a native 2-trait
scenario, then benchmark the new dense annotated 2-trait BayesC implementation
against relevant production baselines.

## Scope

This work adds:

- native 2-trait simulated phenotype and truth files under
  `src/4.Datasets/data/simulated_annotations/`
- documentation of the 2-trait simulation in the generator and README
- a benchmark script in `benchmarks/` that runs the requested methods on the
  packaged fixture
- a benchmark report under `benchmarks/reports/`

It does not change the production multi-trait BayesC algorithm itself.

## 2-Trait Simulation Design

The existing packaged dataset remains the single-trait fixture. The 2-trait
extension is added as new files so current single-trait tests and reports stay
unchanged.

New generated files:

- `phenotypes_mt.csv`
- `annotations_mt.csv`
- `truth_mt.csv`

### Marker states

Use a fixed 4-state truth over the filtered marker set:

- `00`: inactive on both traits
- `10`: trait 1 only
- `01`: trait 2 only
- `11`: active on both traits

Choose nonzero markers in all three active states so the benchmark actually
exercises:

- active vs null
- pleiotropic vs singleton
- trait-1-only vs trait-2-only

### Effect-size generation

For active markers, draw effects with the same MAF-scaled style used by the
single-trait generator.

- singleton markers: one nonzero trait effect, the other fixed to zero
- pleiotropic markers: two nonzero effects with shared positive correlation

### Annotation design

Use a marker-level annotation matrix that contains:

- an activity-enrichment signal
- a pleiotropy-enrichment signal
- a singleton-direction signal
- a noise annotation

This is deliberate: the new 2-trait annotated BayesC prior has three binary
steps, so the simulation should supply signal for all three steps rather than
only the first one.

### Phenotype generation

Construct two genetic values from the same genotype matrix and add correlated
residual noise so the saved phenotype file contains:

- `ID`
- `y1`
- `y2`

## Benchmark Method Matrix

The benchmark should run the dense production path only.

Required baselines:

- `MT_BayesC`
- `MT_Annotated_BayesC`
- `BayesC_y1`
- `Annotated_BayesC_y1`
- `BayesC_y2`
- `Annotated_BayesC_y2`
- `BayesR_y1`
- `BayesR_y2`

Recommended additional reference rows:

- `Annotated_BayesR_y1`
- `Annotated_BayesR_y2`

The additional annotated BayesR runs are useful because they are the closest
existing annotation baseline in JWAS, even though they were not explicitly
required.

## Benchmark Metrics

For each run, compute:

- runtime
- trait-wise `cor(y, EBV)`
- trait-wise correlation between posterior mean marker effects and true marker
  effects
- trait-wise top-k recall using PIP ranking, with `k` equal to the true number
  of active markers for that trait
- any-active top-k recall using `max(PIP_y1, PIP_y2)` for multi-trait methods
  and the corresponding single-trait runs where applicable
- annotation coefficient summaries for annotated methods

All marker-level comparisons must join on normalized marker IDs rather than row
order.

## Reproducibility

Because this is a stochastic benchmark, use multiple seeds rather than a single
short chain. The benchmark script should expose seeds and MCMC length via
environment variables, then write:

- per-run summaries
- aggregate method summaries

## Deliverables

- generator and README updates
- packaged 2-trait fixture files
- dataset access test coverage
- benchmark script
- benchmark report
- implementation note
