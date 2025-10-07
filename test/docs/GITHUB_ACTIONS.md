# GitHub Actions CI/CD for JWAS.jl

Automated testing on every push - free for public repos!

---

## Quick Setup (3 Commands)

```bash
# 1. Copy workflow file
mkdir -p .github/workflows
cp test/docs/github_actions_example.yml .github/workflows/CI.yml

# 2. Commit
git add .github/workflows/CI.yml
git commit -m "Add GitHub Actions CI"

# 3. Push and watch it run!
git push origin master
# Go to: https://github.com/reworkhow/JWAS.jl/actions
```

**Done!** Tests now run automatically on every push. 🎉

---

## What You Get

✅ **Automatic testing** on every push/PR  
✅ **6 parallel test jobs:**
- Ubuntu + Julia 1.8
- Ubuntu + Julia 1.11
- macOS + Julia 1.8
- macOS + Julia 1.11
- Windows + Julia 1.8
- Windows + Julia 1.11

✅ **Results in 2-5 minutes**  
✅ **Email notifications** if tests fail  
✅ **Free forever** (public repos)  
✅ **Zero maintenance**

---

## How It Works

```
You push code to GitHub
         ↓
GitHub Actions triggers
         ↓
Spawns 6 test environments in parallel
         ↓
Each environment:
  1. Checks out your code
  2. Installs Julia
  3. Installs JWAS dependencies
  4. Runs test/runtests.jl
  5. Reports pass/fail
         ↓
Results shown on GitHub
         ↓
✅ All pass → Safe to merge
❌ Some fail → Shows which OS/version failed
```

**Time:** 2-5 minutes for all 6 jobs (run in parallel)

---

## What It Tests

The workflow runs your unit tests (`test/runtests.jl`):
- All 63 unit tests
- On 3 operating systems
- On 2 Julia versions
- Total: 63 tests × 6 configurations = 378 test runs per push!

---

## Viewing Results

### On GitHub Actions Page
```
https://github.com/reworkhow/JWAS.jl/actions

✅ CI - Add new feature (#123)
   ├─ ✅ Julia 1.8 - Ubuntu
   ├─ ✅ Julia 1.8 - macOS
   ├─ ✅ Julia 1.8 - Windows
   ├─ ✅ Julia 1.11 - Ubuntu
   ├─ ✅ Julia 1.11 - macOS
   └─ ✅ Julia 1.11 - Windows
```

### On Pull Requests
```
Pull Request #45: Add BayesD method

Checks:
✅ CI / Julia 1.8 - Ubuntu
✅ CI / Julia 1.8 - macOS
✅ CI / Julia 1.8 - Windows
✅ All checks passed!

[Merge Pull Request] ← Safe to merge
```

---

## Customization

Edit `.github/workflows/CI.yml`:

### Test Only on Ubuntu
```yaml
matrix:
  os:
    - ubuntu-latest
```

### Test Only Julia 1.11
```yaml
matrix:
  version:
    - '1'
```

### Add More Julia Versions
```yaml
matrix:
  version:
    - '1.8'
    - '1.10'
    - '1'
```

### Run Integration Tests on Master
```yaml
- name: Integration tests
  if: github.ref == 'refs/heads/master'
  run: julia --project=. test/integration/test_BayesianAlphabet.jl
```

---

## Cost

**FREE for public repositories!**

GitHub Actions free tier:
- ✅ Unlimited minutes for public repos
- ✅ 2,000 minutes/month for private repos
- ✅ 20 concurrent jobs

For JWAS.jl (public): **$0/month** 💰

---

## Add Status Badge

Add to your README.md:

```markdown
[![CI](https://github.com/reworkhow/JWAS.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/reworkhow/JWAS.jl/actions/workflows/CI.yml)
```

Shows: ![CI](https://img.shields.io/badge/build-passing-brightgreen)

---

## Benefits

| Without CI | With CI |
|------------|---------|
| Manual testing | Automatic |
| Test on your OS only | Test on 3 OS × 2 Julia versions |
| Hope nothing breaks | Know immediately |
| Find issues after release | Catch before merge |
| Your time | Robot time |

---

## Workflow File Structure

The workflow file (`github_actions_example.yml`) has these sections:

### 1. Triggers
```yaml
on:
  push:
    branches: [master]
  pull_request:
    branches: [master]
```
Runs on pushes to master and on all pull requests.

### 2. Test Matrix
```yaml
matrix:
  version: ['1.8', '1']
  os: [ubuntu-latest, macOS-latest, windows-latest]
```
Creates 6 jobs (2 versions × 3 OS).

### 3. Steps
```yaml
- Checkout code
- Setup Julia
- Cache dependencies (2-10x faster!)
- Install packages
- Run tests
- Upload coverage (optional)
```

### 4. Caching
Speeds up subsequent runs from ~5 min to ~2 min by caching Julia packages.

---

## Example: Pull Request Workflow

```
Day 1: Create PR
├─ Push code
└─ GitHub Actions: ✅ All tests pass (2 min)

Day 2: Make changes
├─ Push update
└─ GitHub Actions: ❌ Windows + Julia 1.11 failed
    ├─ You: Fix Windows issue
    └─ Push fix

Day 3: Final review
├─ Push last change
└─ GitHub Actions: ✅ All tests pass
    └─ Reviewer: "LGTM! ✅" → Merge
```

**Benefit:** Caught Windows issue before merge!

---

## Protecting Master Branch

In GitHub settings:
```
Settings → Branches → Add rule
✅ Require status checks before merging
Select: All CI checks
```

**Result:** Cannot merge if tests fail! Safety guaranteed.

---

## FAQ

**Q: Does it cost money?**  
A: No! Free for public repos.

**Q: How long does it take?**  
A: 2-5 minutes (all 6 jobs run in parallel).

**Q: Can I run tests manually?**  
A: Yes! Actions tab → Select workflow → "Run workflow"

**Q: What if I only want Ubuntu?**  
A: Edit the workflow file, remove other OS from matrix.

**Q: Will it slow down my pushes?**  
A: No! Tests run in background. You can keep working.

**Q: What if tests fail?**  
A: GitHub emails you + shows which OS/version failed.

---

## Summary

GitHub Actions = **Free robot team** that:
- Tests your code on 6 configurations
- Every time you push
- Reports in 2-5 minutes
- Catches issues before they reach users
- Zero effort after 3-command setup

**Highly recommended for JWAS.jl!** 🚀

---

## Next Steps

1. Copy workflow file to `.github/workflows/CI.yml`
2. Push to GitHub
3. Watch tests run automatically
4. Add status badge to README
5. Enable branch protection (optional)

See `github_actions_example.yml` in this folder for the complete workflow template.

