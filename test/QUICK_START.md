# JWAS.jl Testing - Quick Start

## Run Tests (3 Modes)

### Fast (20 seconds) - Daily Development ⚡
```bash
julia --project=. test/runtests.jl
```

### With Integration Tests (minutes-hours)
```bash
RUN_INTEGRATION_TESTS=true julia --project=. test/runtests.jl
```

### With Long Tests (hours)
```bash
RUN_LONG_TESTS=true julia --project=. test/runtests.jl
```

### Everything
```bash
RUN_INTEGRATION_TESTS=true RUN_LONG_TESTS=true julia --project=. test/runtests.jl
```

---

## Test Folder Structure

```
test/
├── runtests.jl           ⭐ Main test (63 unit tests)
├── integration/          Integration tests (3 files)
├── long/                 Long-running tests (3 files)
└── docs/                 Documentation (5 files)
```

---

## Documentation

- **`docs/TESTING_GUIDE.md`** - Complete testing guide
- **`docs/TESTING_PATTERNS.md`** - Code examples
- **`docs/GITHUB_ACTIONS.md`** - CI/CD setup

---

## That's It!

For most development: `julia --project=. test/runtests.jl` 🚀

