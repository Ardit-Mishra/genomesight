# Deployment — GenomeSight

| | |
|---|---|
| **Production domain** | https://genomesight.arditmishra.com |
| **Repository** | `github.com/Ardit-Mishra/genomesight` |
| **Default branch** | `main` |
| **Hosting** | Replit Autoscale (live responses carry `Server: Google Frontend`) |
| **Run command** | `streamlit run main.py --server.port 5000` |
| **Deploy trigger** | Replit Deployment. ⚠️ **Not yet verified** whether a GitHub push redeploys automatically or a manual redeploy is required — confirm before relying on a push to ship. |
| **Environment variables** | None required. |
| **Last verified** | 2026-08-28 |

## Rollback

`main` is the deploying branch, so rollback is `git revert <sha>` (preferred — keeps history) or
resetting `main` to the last known-good commit and redeploying. The pre-reconciliation state of the
local copy is preserved on the branch `backup/pre-reconcile-2026-08-28` (`ae73663`).

## Run locally

The system Python is 3.14, which has **no prebuilt `pyarrow` wheel** — `uv run` on 3.14 tries to
compile it from source and fails with `error: command 'cmake' failed`. Pin 3.12:

```bash
uv run --python 3.12 streamlit run main.py --server.port 8502
```

Tests (pytest is not in the lockfile, so bring it in for the run):

```bash
uv run --python 3.12 --with pytest python -m pytest tests/ -q
```

## Reconciliation note (2026-08-28)

The local copy had diverged from `origin/main`: **13 commits behind, 1 ahead**. The local-only commit
was a Replit auto-commit ("Published your App") that added a **complete nested duplicate of the project**
under `genomesight/` — byte-identical to the root files, and the "duplicate unused package" the portfolio
audit flagged. Local was reset to `origin/main`, which removed the duplicate and picked up the upstream
cleanup (README updates, banner asset, deletion of `.replit` and `replit.md`).

One genuinely valuable local change was **re-applied** on top: `st.set_page_config()` moved to sit
directly beneath `import streamlit as st`. Streamlit requires it to be the first Streamlit command
executed; upstream had it *after* the `app.*` imports, which works only for as long as none of those
modules makes a module-level `st.*` call. The page title was also corrected to **GenomeSight**.

## Known state

- `app/core/alignment.py` is fully implemented (`align_sequences`, `calculate_identity`,
  `format_alignment_display`) but is **not imported by `main.py`** — the feature is unreachable from the
  UI. `main.py`'s docstring deliberately does not list it.
- In-page branding still reads "Genome Sequencing Analyzer" (heading, footer, `app/core/export.py`
  report metadata) while the tab title and domain say **GenomeSight**.
- Hosted on an autoscaling tier → a first request after idle carries a cold-start delay.
