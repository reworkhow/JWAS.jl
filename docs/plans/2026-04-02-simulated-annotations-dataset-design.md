# Simulated Annotations Dataset Design

## Goal

Add one permanent JWAS dataset fixture based on Jian's Practical-style annotation simulation so the annotation methods can be tested and benchmarked against a fixed saved dataset instead of regenerating data ad hoc.

## Scope

This change adds:

- a new packaged dataset directory `src/4.Datasets/data/simulated_annotations/`
- the packaged raw genotype input used by the simulation:
  - `raw_genotypes.txt`
- the saved dataset files:
  - `genotypes.csv`
  - `phenotypes.csv`
  - `annotations.csv`
  - `truth.csv`
- the exact generator script used to create those files
- a short dataset README
- a small test that the new dataset is discoverable through `Datasets.dataset(...)`

This change does not modify any production annotation sampler logic.

## Dataset Definition

The dataset should be the fixed output of Jian's Practical-style simulation setup:

- genotype source copied into the dataset folder as `raw_genotypes.txt`
- simulation seed: `123`
- `h2 = 0.5`
- `ncv = 10`
- `S = -1`
- training sample size: `400`
- annotation columns:
  - `functional`
  - `random_anno`

The saved files should match the already-used Practical-style export shape:

- `raw_genotypes.txt`: raw genotype matrix used as simulation input
- `genotypes.csv`: `ID` plus marker columns
- `phenotypes.csv`: `ID`, `y1`
- `annotations.csv`: `marker_id`, `functional`, `random_anno`
- `truth.csv`: `marker_id`, `is_causal`, `true_effect`

## Location Choice

The correct permanent home is the packaged datasets tree under `src/4.Datasets/data/`, not `benchmarks/debug/out/`.

Reason:

- this dataset is now intended as a reusable fixture
- `Datasets.dataset(...)` already resolves files from that tree
- tests and benchmarks can reuse one canonical saved dataset
- it avoids accidental drift from repeated regeneration in temporary folders

## Testing

The minimum test is:

- `Datasets.dataset("genotypes.csv", dataset_name="simulated_annotations")` resolves
- same for `phenotypes.csv`, `annotations.csv`, `truth.csv`
- and the self-contained inputs `raw_genotypes.txt` and `generate_dataset.R`
- the files exist and have the expected key columns

That is enough for this pass because the change is a data fixture addition, not a sampler change.
