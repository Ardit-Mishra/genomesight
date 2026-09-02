import logging

logger = logging.getLogger(__name__)

def self_ping_status() -> str:
    """Return status for health check and self-ping monitoring."""
    return "active"
