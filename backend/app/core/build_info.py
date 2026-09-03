"""What is actually running, so a stale deploy is a lookup rather than an investigation.

A deployment once served a commit five behind the repository for hours. Every
symptom pointed at the code being wrong -- alignment returning 27% where the
truth was 90% -- when the code was right and the deploy was old. Nothing served
by the API said which commit it was built from, so the only way to tell was to
compare behaviour against local runs, endpoint by endpoint.

This module removes that entire class of confusion. /api/health now reports the
commit, so "is production current?" is answered by reading one field.

Resolution order, most to least trustworthy:
  1. GIT_COMMIT baked in at image build time (survives anywhere the image runs)
  2. RENDER_GIT_COMMIT / RENDER_GIT_BRANCH, which Render injects at runtime
  3. git rev-parse, for local development
  4. "unknown" -- reported honestly rather than guessed
"""
from __future__ import annotations

import os
import subprocess
from functools import lru_cache


def _from_env() -> tuple[str | None, str | None]:
    commit = os.environ.get("GIT_COMMIT") or os.environ.get("RENDER_GIT_COMMIT")
    branch = os.environ.get("GIT_BRANCH") or os.environ.get("RENDER_GIT_BRANCH")
    return (commit or None), (branch or None)


def _from_git() -> tuple[str | None, str | None]:
    """Local development only. Never raises -- a missing git is not an error here."""
    def run(args: list[str]) -> str | None:
        try:
            out = subprocess.run(args, capture_output=True, text=True, timeout=5)
            return out.stdout.strip() or None if out.returncode == 0 else None
        except Exception:
            return None

    return run(["git", "rev-parse", "HEAD"]), run(["git", "rev-parse", "--abbrev-ref", "HEAD"])


@lru_cache(maxsize=1)
def build_info() -> dict:
    commit, branch = _from_env()
    if not commit:
        commit, branch = _from_git()

    return {
        # Full SHA, not shortened: this is compared against `git rev-parse HEAD`,
        # and an abbreviation invites a mismatch that looks like a difference.
        "commit": commit or "unknown",
        "commit_short": commit[:7] if commit else "unknown",
        "branch": branch or "unknown",
        "source": "environment" if _from_env()[0] else ("git" if commit else "unavailable"),
    }
