# Pushing GenomeSight to GitHub

## First-Time Setup

```bash
git init
git remote add origin https://github.com/YOUR_USERNAME/genomesight.git
```

## Add, Commit, and Push

```bash
git add -A
git commit -m "Initial commit: GenomeSight bioinformatics analyzer"
git branch -M main
git push -u origin main
```

## Subsequent Pushes

```bash
git add -A
git commit -m "Your commit message here"
git push
```

## Verify .gitignore Is Working

Before pushing, confirm no cache files are staged:

```bash
git status
```

You should not see any `__pycache__/`, `*.pyc`, `.pytest_cache/`, `.ruff_cache/`, or `.mypy_cache/` entries.

If cached files were previously tracked, remove them:

```bash
git rm -r --cached __pycache__/ app/core/__pycache__/ app/ui/__pycache__/ app/__pycache__/ tests/__pycache__/ 2>/dev/null
git commit -m "Remove tracked cache files"
git push
```
