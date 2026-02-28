# GitHub Tutorial - The Modern Standard for Collaboration
## 1. Main branch is always
- **main** (or sometimes **master**) contains production-ready, stable code.
- Never commit broken code directly to **main**

## 2. Create a new branch for everychange
(feature, bug fix, documentation, experiment)
- Name it meaningful: **feature/2D_fix_support**, **doc/update-readme**
- Do this locally
```bash
git checkout main
git pull origin main    # alway start from latest main
git checkout -b feature/my-change
```
## 3. Work locally on your branch
- Edit files, add commits as you go (small & frequent is best)
``` bash
git add .
git commit -m "Add 2D plane strain formulation"
git commit -m "Update boundary conditions for 2D cantilever"
```

## 4. Push your branch to Github
```bash
git push -u origin feature/my-change
```
- **-u** sets tracking so future **git push** works without arguments.

## 5. Create a Pull Request (PR) on GitHub

- Go to your repo → Pull requests → New pull request
- Base: **main** ← Compare: **feature/my-change**
- Write a good title & description:

"Add 2D support to geometrically nonlinear SIMP code"

Explain what changed, why, how to test, screenshot if visual (evolution of topology results.)

- Create PR → now teammates see your work.

## 6. Review & discuss
- Collaborators comment on code, suggest improvements.
- You push more commit to the same the same branch → they appear automatically in the PR.
- Run tests/CI (GitHub Actions) if set up - many repos auto-check linting. MATLAB syntax, etc.
## 7. Merge when ready
- Once approved (and tests pass):
Click **Merge pull request** → choose method:
`hightlight` the title of method you want to merge
    - *Merge commit** (default - keeps full history, most common)
    - **Squash and merge** (clean single commit on main - nice for small changes)
    - **Rebase and merge** (linear history)
- Delete the branch after merge (GitHub button makes it ease)
## 8. Pull lastest main locally
```bash
git checkout main
git pull origin main
git branch -d feature/my-change   # clean up local branch
```

## Quick Starter Commands Cheat Sheet (for your next change)

```bash
# Start work
git checkout main
git pull
git checkout -b docs/improve-readme-for-2d

# Work + commit
# ... edit files ...
git add README.md
git commit -m "Clarify 3D → 2D adaptation objective"

# Share
git push -u origin docs/improve-readme-for-2d

# Then: go to GitHub → create PR → merge when happy
```










