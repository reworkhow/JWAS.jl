# JWAS.jl Testing Guide

Complete guide for testing JWAS.jl development version.

---

## Quick Start

```bash
# Run all unit tests (20 seconds)
julia --project=. test/runtests.jl

# That's it! ✅
```

---

## Table of Contents
1. [Test Structure](#test-structure)
2. [Running Tests](#running-tests)
3. [Writing New Tests](#writing-new-tests)
4. [CI/CD Setup](#cicd-setup)

---

## Test Structure

```
test/
├── runtests.jl          ⭐ Main test file (63 unit tests)
├── integration/         📁 Integration tests (minutes)
├── long/                📁 Long-running tests (hours)
└── docs/                📁 This documentation
```

### Unit Tests (runtests.jl)
- **Time:** ~20 seconds
- **Tests:** 63 automated tests with @test assertions
- **Coverage:** All core JWAS functionality
- **Use:** Daily development

### Integration Tests (integration/)
- **Time:** Minutes to hours
- **Tests:** Complete workflow validation
- **Coverage:** All Bayesian methods, single/multi-trait, complete/incomplete data
- **Use:** Before major commits

### Long Tests (long/)
- **Time:** Hours
- **Tests:** Accuracy validation, RRM, single-step analysis
- **Coverage:** Real-world scenarios with prediction accuracy
- **Use:** Before releases

---

## Running Tests

### Mode 1: Unit Tests Only (Default - Fast)
```bash
cd /Users/haocheng/Github/JWAS.jl
julia --project=. test/runtests.jl
```

**Output:**
```
Test Summary:           | Pass  Total   Time
JWAS.jl Full Test Suite |   63     63  22.2s
✓ All tests completed successfully!
```

**Use:** Daily development (20 seconds)

---

### Mode 2: With Integration Tests
```bash
RUN_INTEGRATION_TESTS=true julia --project=. test/runtests.jl
```

**What it does:**
- Runs 63 unit tests
- Plus integration tests from `test/integration/`
- Tests all Bayesian methods comprehensively

**Use:** Before commits (minutes to hours)

**Note:** Currently placeholders - uncomment the `include()` lines in `runtests.jl` to enable.

---

### Mode 3: With Long-Running Tests
```bash
RUN_LONG_TESTS=true julia --project=. test/runtests.jl
```

**What it does:**
- Runs 63 unit tests
- Plus long tests from `test/long/`
- Includes accuracy validation, RRM, single-step analysis

**Use:** Before releases (may take hours)

**Note:** Currently placeholders - uncomment the `include()` lines in `runtests.jl` to enable.

---

### Mode 4: Everything
```bash
RUN_INTEGRATION_TESTS=true RUN_LONG_TESTS=true julia --project=. test/runtests.jl
```

**What it does:**
- All unit tests
- All integration tests
- All long-running tests

**Use:** Final validation before major releases

---

### Running Individual Tests (Alternative)

You can also run test files directly:

```bash
# Individual integration test
julia --project=. test/integration/test_BayesianAlphabet.jl

# Individual long test
julia --project=. test/long/Unitest.jl
```

---

## What Gets Tested

### Unit Tests (63 tests in runtests.jl)

| Component | Tests | What's Tested |
|-----------|-------|---------------|
| **Model Building** | 11 | Single/multi-trait models, covariates, random effects |
| **Genotype Loading** | 19 | All formats, all methods, QC |
| **Pedigree** | 16 | Loading, inbreeding, relationship matrices |
| **MCMC** | 8 | BayesC, RR-BLUP, reproducibility, output |
| **GWAS** | 4 | Model frequency, window-based GWAS |
| **Edge Cases** | 3 | Error handling, invalid inputs |
| **Data Types** | 2 | Single/double precision |

**Methods Tested:** BayesA, BayesB, BayesC, RR-BLUP, BayesL, GBLUP

---

## Unique Test Folder System

Each test run creates an isolated environment:

```
Run 1: test_run_73058/
Run 2: test_run_90754/
Run 3: test_run_52447/
       └─ Random 5-digit number
```

**Benefits:**
- ✅ No file conflicts between runs
- ✅ Safe for parallel execution
- ✅ Clean workspace (auto-deleted after tests)
- ✅ Can keep for debugging (set `CLEANUP_ON_SUCCESS = false`)

---

## Writing New Tests

### Add to Unit Tests (runtests.jl)

Add new testsets following this pattern:

```julia
@testset "My New Feature" begin
    @testset "Basic functionality" begin
        result = my_new_function(input)
        @test result > 0
        @test typeof(result) == Float64
    end
    
    @testset "Edge cases" begin
        @test_throws ErrorException my_new_function(-1)
    end
end
```

See `docs/TESTING_PATTERNS.md` for more examples.

---

## Recommended Workflow

### Daily Development:
```bash
# Make code changes
julia --project=. test/runtests.jl
# 20 seconds - fast feedback!
```

### Before Committing:
```bash
# Run unit tests
julia --project=. test/runtests.jl

# If pass → commit
# If fail → fix and repeat
```

### Before Major Release:
```bash
# 1. Unit tests
julia --project=. test/runtests.jl

# 2. Integration tests
RUN_INTEGRATION_TESTS=true julia --project=. test/runtests.jl

# 3. Long-running tests
RUN_LONG_TESTS=true julia --project=. test/runtests.jl

# All pass → Release!
```

---

## CI/CD Setup

For automated testing on GitHub, see `docs/GITHUB_ACTIONS.md`.

Quick setup:
```bash
mkdir -p .github/workflows
cp test/docs/github_actions_example.yml .github/workflows/CI.yml
git add .github/workflows/CI.yml
git commit -m "Add CI"
git push
```

Then tests run automatically on every push! 🤖

---

## Debugging Failed Tests

### Keep test artifacts:
Edit `test/runtests.jl` line ~392:
```julia
CLEANUP_ON_SUCCESS = false
```

Test outputs will be kept in `test_run_XXXXX/` for inspection.

### Run specific test:
```julia
# In Julia REPL
include("test/runtests.jl")

# Or filter specific testset
# (requires TestSetExtensions package)
```

---

## Test Data

All test data is in `src/4.Datasets/data/example/`:
- `phenotypes.txt` - Phenotype data
- `genotypes.txt` - Genotype matrix
- `pedigree.txt` - Pedigree information
- `map.txt` - Marker positions

Access via:
```julia
using JWAS.Datasets
phenofile = Datasets.dataset("phenotypes.txt", dataset_name="example")
```

---

## Common Test Patterns

### Testing Models
```julia
model = build_model("y = intercept + x1", 1.0)
@test model.nModels == 1
@test :y in model.lhsVec
```

### Testing Genotypes
```julia
geno = get_genotypes("file.txt", 1.0, method="BayesC")
@test geno.nMarkers > 0
@test geno.method == "BayesC"
```

### Testing MCMC
```julia
output = runMCMC(model, data, chain_length=50, seed=123)
@test haskey(output, "location parameters")
@test haskey(output, "residual variance")
```

### Testing Errors
```julia
@test_throws ErrorException build_model("y = intercept", -1.0)
```

See `docs/TESTING_PATTERNS.md` for complete reference.

---

## Summary

| Test Type | File | Time | When |
|-----------|------|------|------|
| **Unit** | `runtests.jl` | 20s | Daily ⚡ |
| **Integration** | `integration/*.jl` | Minutes | Before commits |
| **Long** | `long/*.jl` | Hours | Before releases |

**For most development: Just run `test/runtests.jl`** 🚀

---

## See Also

- `docs/TESTING_PATTERNS.md` - Code examples for writing tests
- `docs/GITHUB_ACTIONS.md` - CI/CD setup guide
- `docs/github_actions_example.yml` - CI workflow template

