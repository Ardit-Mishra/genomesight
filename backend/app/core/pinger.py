"""Keep-alive support and process liveness.

WHAT ACTUALLY KEEPS A FREE-TIER INSTANCE WARM
---------------------------------------------
Render spins an idle free instance down after ~15 minutes, and the next request
pays ~50 seconds of cold start. The obvious fix -- have the app ping itself on a
timer -- CANNOT WORK, and it is worth being precise about why:

  * While the instance is asleep the process is not running, so there is nothing
    left to send the ping. A self-ping can never wake a sleeping service.
  * It would also keep the service awake around the clock, consuming roughly 730
    of the free tier's 750 monthly instance-hours serving traffic nobody asked for.

Only an EXTERNAL request can wake a sleeping instance, so the keep-alive lives
outside the app: .github/workflows/keep-alive.yml calls /api/health on a schedule.
That is the mechanism. This module does not pretend to be it.

What this module does provide is honest liveness data for /api/health: how long
this process has been up, and when it was last reached. A short uptime is the
signature of a cold start that just happened, which explains a slow first request
without anyone having to guess.
"""
from __future__ import annotations

import logging
import time

logger = logging.getLogger(__name__)

_STARTED_AT = time.time()
_last_request_at: float | None = None


def note_request() -> None:
    """Record that this process served a request."""
    global _last_request_at
    _last_request_at = time.time()


def self_ping_status() -> str:
    """Liveness token for the health payload."""
    return "active"


def liveness() -> dict:
    now = time.time()
    uptime = now - _STARTED_AT
    return {
        "status": self_ping_status(),
        "uptime_seconds": round(uptime, 1),
        # Under a minute of uptime almost always means this request just paid for
        # a cold start.
        "likely_cold_start": uptime < 60,
        "seconds_since_last_request": (
            round(now - _last_request_at, 1) if _last_request_at is not None else None
        ),
    }
