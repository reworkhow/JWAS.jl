# Simulated Annotations Dataset Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a permanent packaged `simulated_annotations` dataset fixture based on Jian's Practical-style annotation simulation.

**Architecture:** Save one canonical generated dataset under `src/4.Datasets/data/simulated_annotations/`, keep the generator script and raw genotype source alongside the saved files for reproducibility, and add a small dataset-access test through `JWAS.Datasets`. This change is data-only plus a retrieval test; no sampler code changes are needed.

**Tech Stack:** Julia, R, CSV, DataFrames, JWAS.Datasets, packaged repo data files.

---

### Task 1: Add the failing dataset-access test

**Files:**
- Modify: `test/unit/test_misc_coverage.jl`

**Step 1: Write the failing test**

Add a new testset under `Datasets module` that checks:

- `Datasets.dataset("genotypes.csv", dataset_name="simulated_annotations")`
- `Datasets.dataset("phenotypes.csv", dataset_name="simulated_annotations")`
- `Datasets.dataset("annotations.csv", dataset_name="simulated_annotations")`
- `Datasets.dataset("truth.csv", dataset_name="simulated_annotations")`
- `Datasets.dataset("raw_genotypes.txt", dataset_name="simulated_annotations")`
- `Datasets.dataset("generate_dataset.R", dataset_name="simulated_annotations")`

and asserts the files exist.

**Step 2: Run test to verify it fails**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected:

- failure because the new dataset files do not exist yet

### Task 2: Add the packaged dataset fixture

**Files:**
- Create: `src/4.Datasets/data/simulated_annotations/raw_genotypes.txt`
- Create: `src/4.Datasets/data/simulated_annotations/genotypes.csv`
- Create: `src/4.Datasets/data/simulated_annotations/phenotypes.csv`
- Create: `src/4.Datasets/data/simulated_annotations/annotations.csv`
- Create: `src/4.Datasets/data/simulated_annotations/truth.csv`
- Create: `src/4.Datasets/data/simulated_annotations/README.md`
- Create: `src/4.Datasets/data/simulated_annotations/generate_dataset.R`

**Step 1: Generate the dataset**

Use the existing Practical-style export logic with fixed parameters:

- copy the source genotypes into `raw_genotypes.txt`
- seed `123`
- `h2 = 0.5`
- `ncv = 10`
- `S = -1`
- `n_gwas = 400`

**Step 2: Save the exact generated files**

Place the generated CSVs in `src/4.Datasets/data/simulated_annotations/`.

**Step 3: Save the generator script**

Copy the generation logic into `generate_dataset.R` in the same directory with the fixed parameters documented directly in the file. The script must read the local `raw_genotypes.txt` so the dataset folder is self-contained.

**Step 4: Add a short README**

Document:

- where the source genotypes came from
- that `raw_genotypes.txt` is the packaged local input used for regeneration
- the fixed simulation parameters
- the meaning of each saved file

### Task 3: Verify the dataset-access test passes

**Files:**
- Test: `test/unit/test_misc_coverage.jl`

**Step 1: Run the same test again**

Run:

```bash
julia --project=. --startup-file=no -e 'include("test/unit/test_misc_coverage.jl")'
```

Expected:

- the new dataset access test passes

### Task 4: Verify docs/build remains clean enough for a data-only change

**Files:**
- No doc content changes required

**Step 1: Skip docs build unless docs are touched**

No docs build is required for this dataset-only change.

### Task 5: Commit intentionally

**Files:**
- All files above

**Step 1: Commit**

```bash
git add docs/plans/2026-04-02-simulated-annotations-dataset-design.md \
        docs/plans/2026-04-02-simulated-annotations-dataset-implementation.md \
        src/4.Datasets/data/simulated_annotations \
        test/unit/test_misc_coverage.jl
git commit -m "data: add simulated annotations dataset fixture"
```
