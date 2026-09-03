# Deployment — GenomeSight

| | |
|---|---|
| **Frontend** | https://genomesight-frontend.vercel.app — Vercel, project `genomesight-frontend` |
| **Backend** | https://genomesight-api.onrender.com — Render web service `genomesight-api`, Docker, free tier |
| **Repository** | `Ardit-Mishra/genomesight` |
| **Deploying branch** | `main` |
| **Frontend root** | `frontend/` · build `npm run build` · output `dist/` |
| **Backend root** | `backend/` · `backend/Dockerfile` |
| **Frontend → backend** | `frontend/.env.production` sets `VITE_API_URL` |
| **Last verified** | 2026-09-03 |

## Is production current?

One request answers it:

```
curl -s https://genomesight-api.onrender.com/api/health
```

`commit` is the SHA the running image was built from. Compare it to `git rev-parse HEAD`.
Equal means current; different means the deploy is behind, regardless of what any
dashboard says.

This endpoint exists because of a real incident (below). Before it, telling a stale
deploy from a code defect meant comparing behaviour endpoint by endpoint.

## Deploying

Both platforms deploy on push to `main`.

**Frontend**, if a manual deploy is ever needed:

```
cd frontend && vercel --prod
```

**Backend**: Render rebuilds on commit. A manual redeploy is available in the Render
dashboard under **Manual Deploy → Deploy latest commit**.

## Rollback

- **Frontend** — Vercel keeps every deployment; promote a previous one from the
  dashboard, or `vercel rollback`.
- **Backend** — Render's Events tab has a **Rollback** button on each past successful
  deploy.
- **Source** — `git revert <sha>` on `main`, which is preferred over resetting because
  it keeps the history that explains what happened.

## Incident, 2026-09-03 — read this before debugging a deployment

Production served alignment identity of **27.27%** where the correct answer was
**90.0%**, and returned `count / total` in place of RSCU. Both had been fixed in the
repository hours earlier. Every symptom pointed at the code.

The code was fine. **Render was watching the branch `fastapi-react-rewrite`, pinned at
`e304fc0`, while all work was being pushed to `main`** — by then five commits ahead.
Manual redeploys "worked" and changed nothing, because they faithfully rebuilt the same
stale commit. Auto-Deploy was on; it was simply watching a different branch.

Compounding it, the commit Render was serving carried the message *"complete v2 FastAPI
+ React full-stack refactor with correct RSCU"* — while the code it contained still had
`return {k: v / total for k, v in counts.items()}`. The message asserted a fix the diff
did not make.

Three changes came out of it:

1. `/api/health` reports the deployed commit, branch, and k-mer engine. The question
   "is production current?" is now a lookup rather than an investigation.
2. Both branches were fast-forwarded to the same commit, and Render should be pointed
   at `main` so there is one branch of truth.
3. This file records the real hosts, because the version it replaced described a host
   the project no longer used and admitted its deploy trigger was unverified.

## Known operational limits

- **Free-tier cold start.** The Render instance spins down when idle; the first request
  after a quiet period can take ~50 seconds. The frontend retries with exponential
  backoff, so it degrades gracefully rather than failing — but a first visitor waits.
  `backend/app/core/pinger.py` exists and is not yet wired to anything.
- **Native k-mer accelerator.** Built inside the Docker image. If the build fails the
  image still ships and the service runs the pure-Python path — `/api/health` reports
  `kmer_engine`, so which path is live is always visible rather than assumed.
