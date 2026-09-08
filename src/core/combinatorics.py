"""Canonical combinatorial primitives shared across the repo's packages.

This is the single source of truth for binomial coefficients and the
classical Johnson/Sperner bound family. cwc and cgt previously carried
independent copies of these (some exact, some float-based via scipy, some
camelCase relics of older code) that could silently diverge on large
inputs; everything should import from here instead of redefining them.

This package has no dependency on cwc or cgt -- keep it that way so it can
stay a shared leaf both packages import from, not something they
cross-reference each other through.
"""

import logging
import math

logger = logging.getLogger(__name__)


def binom(n, k):
    """Exact binomial coefficient. Returns 0 for out-of-range k."""
    if k < 0 or k > n or n < 0:
        return 0
    return math.comb(n, k)


def binom_float(n, k):
    """Float binomial coefficient via lgamma, for n/k too large for math.comb."""
    if k < 0 or k > n or n < 0:
        return 0.0
    try:
        return math.exp(math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1))
    except Exception:
        return 0.0


def johnson_bound(m, D, w):
    """Classical Johnson bound: upper bound on A(m, D, w) via the floor recursion.

    Returns the maximum number of codewords for a constant-weight code with
    length m, minimum distance D, and constant weight w.
    """
    ret = 1
    for delta in range(w - D // 2, -1, -1):
        ret = math.floor(ret * (m - delta) / (w - delta))
    return ret


def johnson_bound_without_floor(m, D, w, no_floor=False):
    """Johnson bound computed without intermediate flooring.

    With no_floor=True, returns the exact rational value (as a float) before
    the final floor is applied -- useful when the caller needs to compare the
    bound against a non-integer threshold rather than just the integer bound.
    """
    ret = 1.0
    for delta in range(w - D // 2, -1, -1):
        ret = ret * (m - delta) / (w - delta)
    if no_floor:
        return ret
    return math.floor(ret)


def sperner_bound(n, d):
    """Smallest m such that C(n, d) <= C(m, floor(m/2)) (Sperner's theorem).

    Used by johnson_for_n as a lower-bound starting point for its search.
    """
    nd = binom(n, d)
    for m in range(d, n):
        if nd <= binom(m, m // 2):
            logger.debug("sperner_bound(n=%s, d=%s) -> %s", n, d, m)
            return m
    logger.debug("sperner_bound: no m found for n=%s, d=%s (nd=%s)", n, d, nd)


def johnson(m, d):
    """Smallest n such that a d-cover-free-style Johnson recursion peaks, with its t.

    Returns (n, t) where n is the best bound found by sweeping t and t is the
    t at which the sweep's optimum was attained.
    """
    nmax = None
    t_ = None
    prev = None
    for t in range(m // (2 * d)):
        n = binom(m, t + 1) // binom(t * d + 1, t + 1)
        if nmax is None or n > nmax:
            nmax = n
            t_ = t
        if prev is not None and prev >= n:
            break
        prev = n
    logger.debug("johnson(m=%s, d=%s) -> (%s, %s)", m, d, nmax, t_)
    return nmax, t_


def inverse_q_pow_q(n, d):
    """Smallest q such that n <= q^(q/d)."""
    for q in range(2, n):
        if n <= pow(q, q / d):
            logger.debug("inverse_q_pow_q(n=%s, d=%s) -> %s", n, d, q)
            return q
    logger.debug("inverse_q_pow_q: no q found for n=%s, d=%s", n, d)


def is_isolated_pair(x, y, d):
    r"""Pairwise d-isolation check: |x \cap y| <= floor((min(|x|, |y|) - 1) / d).

    x and y are set-like (any object supporting len() and intersection()).
    """
    intersection_size = len(x.intersection(y))
    threshold = (min(len(x), len(y)) - 1) // d
    return intersection_size <= threshold


def johnson_for_n(n, d):
    """Binary-search the smallest m for which johnson(m, d) yields >= n codewords.

    Returns (m, t) -- the bound m and the t at which johnson(m, d) attained it.
    """
    mlb = sperner_bound(n, d)
    mub = inverse_q_pow_q(n, d) ** 2  # basic RS upper bound
    t = None
    for _ in range(n):
        m = (mub + mlb) // 2
        n_, t_ = johnson(m, d)
        if n_ >= n and mub > m:
            mub = m
            t = t_
        if n_ < n and mlb < m:
            mlb = m
        if mlb + 1 == mub:
            logger.debug("johnson_for_n(n=%s, d=%s) -> (%s, %s)", n, d, mub, t)
            return mub, t
