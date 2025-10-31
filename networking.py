from __future__ import annotations

import logging
from typing import Any, Iterable, Tuple

from requests import Session
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

try:
    import httpx

    HTTPX_AVAILABLE = True
except ImportError:
    HTTPX_AVAILABLE = False

log = logging.getLogger(__name__)

_DEFAULT_ATTEMPTS = 3
_DEFAULT_BACKOFF = 0.5
_DEFAULT_STATUS = (500, 502, 503, 504)
_USE_HTTP2 = False
_SESSION_CACHE: dict[
    Tuple[int, float, Tuple[int, ...], bool], Session | httpx.Client
] = {}


def set_http2_mode(use_http2: bool) -> None:
    """Set the global HTTP/2 mode for networking (for pipeline integration)."""
    global _USE_HTTP2
    if use_http2 and not HTTPX_AVAILABLE:
        log.warning(
            "HTTP/2 requested but httpx is not installed. Falling back to HTTP/1.1.")
        _USE_HTTP2 = False
    else:
        _USE_HTTP2 = use_http2
    _SESSION_CACHE.clear()


def _build_httpx_client(attempts: int, backoff: float, status: Tuple[int, ...]) -> 'httpx.Client':
    if not HTTPX_AVAILABLE:
        raise ImportError("httpx is not installed.")
    limits = httpx.Limits(max_connections=100, max_keepalive_connections=20)
    transport = httpx.HTTPTransport(
        retries=attempts, http2=_USE_HTTP2, limits=limits)
    client = httpx.Client(http2=_USE_HTTP2, transport=transport, timeout=30.0)
    return client


def _build_session(attempts: int, backoff: float, status: Tuple[int, ...]) -> Session:
    retry = Retry(
        total=attempts,
        backoff_factor=backoff,
        status_forcelist=status,
        allowed_methods=("GET", "POST"),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(pool_connections=20,
                          pool_maxsize=100, max_retries=retry)
    session = Session()
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def configure_retry(attempts: int, backoff: float, status: Iterable[int]) -> None:
    """Update the default retry strategy used by cached HTTP sessions.

    Args:
        attempts (int): Maximum number of retry attempts; coerced to at least ``1``.
        backoff (float): Backoff factor in seconds; negative values become ``0``.
        status (Iterable[int]): HTTP status codes that should trigger a retry.
    """
    global _DEFAULT_ATTEMPTS, _DEFAULT_BACKOFF, _DEFAULT_STATUS, _SESSION_CACHE
    _DEFAULT_ATTEMPTS = max(1, attempts)
    _DEFAULT_BACKOFF = max(0.0, backoff)
    _DEFAULT_STATUS = tuple(int(code) for code in status)
    _SESSION_CACHE.clear()


def get_session(
    attempts: int | None = None,
    backoff: float | None = None,
    status: Iterable[int] | None = None,
) -> Session | httpx.Client:
    """Return a cached ``requests.Session`` configured with retry behaviour.

    Args:
        attempts (int | None): Override for the retry attempt count.
        backoff (float | None): Override for the retry backoff factor in seconds.
        status (Iterable[int] | None): Override for HTTP status codes retried.

    Returns:
        Session | httpx.Client: Shared session instance with the configured retry policy.
    """
    resolved_attempts = max(
        1, attempts) if attempts is not None else _DEFAULT_ATTEMPTS
    resolved_backoff = max(
        0.0, backoff) if backoff is not None else _DEFAULT_BACKOFF
    resolved_status = (
        tuple(int(code)
              for code in status) if status is not None else _DEFAULT_STATUS
    )
    use_http2 = _USE_HTTP2
    cache_key = (resolved_attempts, resolved_backoff,
                 resolved_status, use_http2)
    session = _SESSION_CACHE.get(cache_key)
    if session is None:
        if use_http2:
            session = _build_httpx_client(
                resolved_attempts, resolved_backoff, resolved_status
            )
        else:
            session = _build_session(
                resolved_attempts, resolved_backoff, resolved_status
            )
        _SESSION_CACHE[cache_key] = session
    return session


__all__ = ["configure_retry", "get_session",
           "HTTPX_AVAILABLE", "set_http2_mode"]
