# GitHub Actions for JWAS.jl - Visual Guide

## What Happens When You Push Code?

```
┌─────────────────────────────────────────────────────────────────┐
│                    YOU (Local Machine)                          │
│                                                                 │
│  1. Make code changes to JWAS.jl                               │
│  2. git add .                                                   │
│  3. git commit -m "Add new feature"                            │
│  4. git push origin master                                      │
│                                                                 │
└─────────────────────┬───────────────────────────────────────────┘
                      │
                      │ Push triggers GitHub Actions
                      ↓
┌─────────────────────────────────────────────────────────────────┐
│                    GITHUB.COM                                   │
│                                                                 │
│  GitHub receives your push and checks:                         │
│  "Is there a .github/workflows/CI.yml file?"                   │
│                                                                 │
│  YES → Start automated testing                                 │
│                                                                 │
└─────────────────────┬───────────────────────────────────────────┘
                      │
                      │ Spawns multiple test runners
                      ↓
┌─────────────────────────────────────────────────────────────────┐
│              GITHUB ACTIONS RUNNERS (Cloud VMs)                 │
│                                                                 │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐           │
│  │  Ubuntu     │  │   macOS     │  │  Windows    │           │
│  │  Julia 1.8  │  │  Julia 1.8  │  │  Julia 1.8  │           │
│  └─────────────┘  └─────────────┘  └─────────────┘           │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐           │
│  │  Ubuntu     │  │   macOS     │  │  Windows    │           │
│  │  Julia 1.11 │  │  Julia 1.11 │  │  Julia 1.11 │           │
│  └─────────────┘  └─────────────┘  └─────────────┘           │
│                                                                 │
│  Each runner executes:                                         │
│  1. Checkout your code                                         │
│  2. Install Julia                                              │
│  3. Install JWAS dependencies                                  │
│  4. Run unit tests                                             │
│  5. (Optional) Run integration tests                           │
│  6. Report results                                             │
│                                                                 │
└─────────────────────┬───────────────────────────────────────────┘
                      │
                      │ Results sent back
                      ↓
┌─────────────────────────────────────────────────────────────────┐
│                    GITHUB.COM - Results                         │
│                                                                 │
│  ✅ Ubuntu + Julia 1.8:  All tests passed                      │
│  ✅ Ubuntu + Julia 1.11: All tests passed                      │
│  ✅ macOS + Julia 1.8:   All tests passed                      │
│  ✅ macOS + Julia 1.11:  All tests passed                      │
│  ✅ Windows + Julia 1.8: All tests passed                      │
│  ❌ Windows + Julia 1.11: 2 tests failed                       │
│                                                                 │
│  Overall Status: ❌ FAILED                                      │
│                                                                 │
└─────────────────────┬───────────────────────────────────────────┘
                      │
                      │ Email notification sent
                      ↓
┌─────────────────────────────────────────────────────────────────┐
│                    YOUR EMAIL                                   │
│                                                                 │
│  Subject: [JWAS.jl] Tests failed on master                     │
│                                                                 │
│  Your push to master failed CI checks:                         │
│  - Windows + Julia 1.11: FAILED                                │
│                                                                 │
│  Click here to see details: [link]                             │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## Step-by-Step: What Each Part Does

### **1. Trigger Events (Lines 6-14)**

```yaml
on:
  push:
    branches:
      - master
      - main
  pull_request:
    branches:
      - master
      - main
```

**What it means:**
- Run tests when pushing to `master` or `main` branch
- Run tests when someone creates a Pull Request to those branches
- Doesn't run on feature branches (optional - you can add them)

---

### **2. Test Matrix (Lines 20-29)**

```yaml
matrix:
  version:
    - '1.8'   # LTS version
    - '1'     # Latest stable
  os:
    - ubuntu-latest
    - macOS-latest
    - windows-latest
```

**What it means:**
- Tests run on **6 combinations**: 
  - 2 Julia versions × 3 operating systems = 6 parallel jobs
- If one fails, you know exactly where (e.g., "Julia 1.11 on Windows")

**Visual:**
```
         Ubuntu    macOS    Windows
Julia 1.8   ✅       ✅        ✅
Julia 1.11  ✅       ✅        ❌  ← This one failed!
```

---

### **3. Caching (Lines 40-46)**

```yaml
- name: Cache Julia artifacts
  uses: actions/cache@v3
  with:
    path: ~/.julia
    key: ${{ runner.os }}-julia-${{ matrix.version }}-${{ hashFiles('**/Project.toml') }}
```

**What it means:**
- First run: Downloads all packages (~5 minutes)
- Later runs: Reuses cached packages (~30 seconds)
- **2-10x faster!**

**Time Comparison:**
```
Without cache:
  Install packages: 5 min
  Run tests:       2 min
  Total:           7 min

With cache:
  Load cache:      30 sec
  Run tests:       2 min
  Total:           2.5 min
```

---

### **4. Conditional Tests (Lines 58-62)**

```yaml
- name: Run integration tests (only on master)
  if: github.ref == 'refs/heads/master' || github.ref == 'refs/heads/main'
  env:
    RUN_INTEGRATION_TESTS: "true"
  run: julia --project=. test/runtests_example.jl
```

**What it means:**
- **Feature branches**: Only run fast unit tests (< 1 minute)
- **Master/Main branch**: Run full integration tests (might take hours)

**Decision Tree:**
```
Is this a push to master/main?
├─ YES → Run ALL tests (unit + integration)
└─ NO  → Run only unit tests (fast feedback)
```

---

## Real-World Example: Pull Request Workflow

Let's say you want to add a new feature to JWAS:

### **Scenario: Adding BayesD Method**

```
Day 1: Create feature branch
├─ You: git checkout -b add-bayesD
├─ You: Write code for BayesD
├─ You: git push origin add-bayesD
└─ GitHub Actions:
    ├─ Runs unit tests only (fast!)
    ├─ Takes 2 minutes
    └─ Result: ✅ All pass

Day 2: Fix a bug
├─ You: Fix typo in BayesD
├─ You: git push origin add-bayesD
└─ GitHub Actions:
    ├─ Runs unit tests again
    ├─ Takes 2 minutes
    └─ Result: ✅ All pass

Day 3: Create Pull Request
├─ You: Create PR: "Add BayesD method"
└─ GitHub Actions:
    ├─ Runs unit tests (6 jobs in parallel)
    ├─ Result: ✅ Ubuntu: Pass, ✅ macOS: Pass, ✅ Windows: Pass
    └─ PR shows green checkmark ✅
    
Reviewer: "Looks good! ✅"

Day 4: Merge to master
├─ You: Click "Merge Pull Request"
└─ GitHub Actions:
    ├─ Runs FULL test suite (unit + integration)
    ├─ Takes 30 minutes
    └─ Result: ✅ All tests pass
    
Master branch updated with new feature! 🎉
```

---

## Benefits for JWAS.jl

### **Before GitHub Actions:**
```
You: Make changes
You: Manually run tests locally
You: Hope you ran all tests
You: Push to GitHub
Collaborator: Pulls code
Collaborator: Code breaks on their machine
Collaborator: "This doesn't work on Windows!"
You: "Works on my Mac..." 😞
```

### **After GitHub Actions:**
```
You: Make changes
You: Push to GitHub
GitHub: Automatically tests on Ubuntu, macOS, Windows
GitHub: Automatically tests on Julia 1.8, 1.11
GitHub: "✅ All tests pass on all platforms!"
You: Confident merge
Collaborator: Pulls code
Collaborator: Works perfectly!
Everyone: Happy! 😊
```

---

## Setting It Up (3 Steps!)

### **Step 1: Create the file**
```bash
cd /Users/haocheng/Github/JWAS.jl
mkdir -p .github/workflows
cp test/github_actions_example.yml .github/workflows/CI.yml
```

### **Step 2: Commit and push**
```bash
git add .github/workflows/CI.yml
git commit -m "Add GitHub Actions CI"
git push origin master
```

### **Step 3: Watch it run!**
1. Go to https://github.com/reworkhow/JWAS.jl/actions
2. See your tests running in real-time
3. Get automatic emails if tests fail

---

## Viewing Results

### **On GitHub:**
```
https://github.com/reworkhow/JWAS.jl/actions

┌────────────────────────────────────────┐
│  Workflows                             │
├────────────────────────────────────────┤
│  ✅ CI - Add BayesD method             │
│     #123 - 2 minutes ago               │
│     ├─ ✅ Julia 1.8 - Ubuntu           │
│     ├─ ✅ Julia 1.8 - macOS            │
│     ├─ ✅ Julia 1.8 - Windows          │
│     ├─ ✅ Julia 1.11 - Ubuntu          │
│     ├─ ✅ Julia 1.11 - macOS           │
│     └─ ✅ Julia 1.11 - Windows         │
│                                        │
│  ❌ CI - Fix bug                       │
│     #122 - 1 hour ago                  │
│     ├─ ✅ Julia 1.8 - Ubuntu           │
│     ├─ ✅ Julia 1.8 - macOS            │
│     ├─ ❌ Julia 1.8 - Windows          │
│     └─ ... (3 more)                    │
└────────────────────────────────────────┘
```

### **On Pull Requests:**
```
Pull Request #45: Add BayesD method

┌────────────────────────────────────────┐
│  Checks                                │
├────────────────────────────────────────┤
│  ✅ CI / Julia 1.8 - Ubuntu            │
│  ✅ CI / Julia 1.8 - macOS             │
│  ✅ CI / Julia 1.8 - Windows           │
│  ✅ CI / Julia 1.11 - Ubuntu           │
│  ✅ CI / Julia 1.11 - macOS            │
│  ✅ CI / Julia 1.11 - Windows          │
│                                        │
│  All checks passed!                    │
│                                        │
│  [Merge Pull Request]  ← Safe to merge│
└────────────────────────────────────────┘
```

---

## Advanced Features

### **1. Add a Badge to README**
```markdown
[![CI](https://github.com/reworkhow/JWAS.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/reworkhow/JWAS.jl/actions/workflows/CI.yml)
```

Shows: ![CI](https://img.shields.io/badge/build-passing-brightgreen)

### **2. Block Merges if Tests Fail**
In GitHub repository settings:
```
Settings → Branches → Branch protection rules
✅ Require status checks to pass before merging
✅ Require branches to be up to date before merging
Select: CI / Julia 1.8 - Ubuntu
        CI / Julia 1.8 - macOS
        CI / Julia 1.8 - Windows
        ... (all 6 checks)
```

Now you **cannot** merge if tests fail!

### **3. Test on Pull Requests from Forks**
```yaml
on:
  pull_request_target:  # Allows testing PRs from forks
    branches:
      - master
```

---

## Cost

**Good news: GitHub Actions is FREE for public repositories!**

Free tier includes:
- ✅ Unlimited minutes for public repos
- ✅ 2,000 minutes/month for private repos
- ✅ 20 concurrent jobs

For JWAS.jl (public repo): **$0/month**

---

## Summary

```
┌─────────────────────────────────────────────────────────────┐
│  Without GitHub Actions         With GitHub Actions         │
├─────────────────────────────────────────────────────────────┤
│  Manual testing                 Automatic testing           │
│  Test on your machine only      Test on 3 OS, 2 Julia vers │
│  Miss edge cases                Catch issues early          │
│  Hope nothing breaks            Know immediately if broken  │
│  Slow feedback                  Fast feedback               │
│  Manual effort                  Zero manual effort          │
│  Costs: Your time               Costs: $0                   │
└─────────────────────────────────────────────────────────────┘
```

**Bottom Line:**
1. You push code
2. GitHub automatically tests it on 6 different configurations
3. You get results in 2-5 minutes
4. Completely free for public repos
5. Zero configuration after initial setup

It's like having 6 robot assistants testing your code 24/7! 🤖✨

