# simulated_annotations

Packaged annotation benchmark fixture derived from Jian's Practical-style
simulation setup for annotated Bayes methods.

This folder is self-contained. It includes:

- the raw genotype source used by the simulation
- the exact generator script
- the generated benchmark CSV files

## Source

- Raw genotype source: `raw_genotypes.txt`
- Generator: [generate_dataset.R](generate_dataset.R)

## Fixed simulation settings

- Simulation seed: `123`
- Heritability: `0.5`
- Number of causal SNPs: `10`
- Effect-size exponent `S`: `-1`
- Training sample size: `400`
- Minor-allele-frequency filter: `0.01 < MAF < 0.99`

## Multi-trait extension

This packaged fixture now also includes a native 2-trait scenario for dense
multi-trait annotation benchmarks.

The 2-trait simulation uses the same filtered genotype matrix and training
sample size as the single-trait fixture, but writes separate files so the
single-trait benchmark remains unchanged.

### 2-trait truth

Markers are assigned to one of four inclusion states:

- `00`: inactive on both traits
- `10`: active on trait 1 only
- `01`: active on trait 2 only
- `11`: active on both traits

The generator currently uses:

- `8` shared (`11`) markers
- `6` trait-1-only (`10`) markers
- `6` trait-2-only (`01`) markers

### 2-trait effects

- singleton markers get one nonzero MAF-scaled effect and one zero effect
- shared markers get two MAF-scaled effects with positive correlation
- residual noise is sampled with positive cross-trait covariance

### 2-trait annotations

The multi-trait annotation file contains four marker-level columns:

- `active_signal`: enriched for any active marker (`10`, `01`, `11`)
- `pleiotropy_signal`: enriched for shared markers (`11`)
- `direction_signal`: positive for `10`, negative for `01`
- `random_signal`: uninformative noise

This matches the 3-step tree prior used by dense 2-trait annotated BayesC:

1. `00` vs active
2. `11` vs singleton
3. `10` vs `01`

## Saved files

- `raw_genotypes.txt`
  - raw genotype matrix used as the simulation input
  - copied into this folder so regeneration does not depend on `temp/`
- `genotypes.csv`
  - `400` individuals by `964` markers
  - first column `ID`, remaining columns `m1` to `m964`
- `phenotypes.csv`
  - training phenotypes with columns `ID`, `y1`
- `annotations.csv`
  - marker-level annotations with columns:
    - `marker_id`
    - `functional`
    - `random_anno`
- `truth.csv`
  - marker-level simulation truth with columns:
    - `marker_id`
    - `is_causal`
    - `true_effect`
- `phenotypes_mt.csv`
  - two-trait phenotype file with columns `ID`, `y1`, `y2`
- `annotations_mt.csv`
  - two-trait marker annotations with columns:
    - `marker_id`
    - `active_signal`
    - `pleiotropy_signal`
    - `direction_signal`
    - `random_signal`
- `truth_mt.csv`
  - two-trait marker-level truth with columns:
    - `marker_id`
    - `state`
    - `is_active_y1`
    - `is_active_y2`
    - `is_shared`
    - `true_effect_y1`
    - `true_effect_y2`

## Regeneration

From the repository root:

```bash
Rscript src/4.Datasets/data/simulated_annotations/generate_dataset.R
```

The script reads `raw_genotypes.txt` in this directory and overwrites both the
single-trait and 2-trait generated CSV files in the same directory.
