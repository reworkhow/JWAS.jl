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

## Regeneration

From the repository root:

```bash
Rscript src/4.Datasets/data/simulated_annotations/generate_dataset.R
```

The script reads `raw_genotypes.txt` in this directory and overwrites the four
generated CSV files in the same directory.
