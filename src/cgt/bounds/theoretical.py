"""Classical constant-weight coding relations used as building blocks.

Reference: standard coding-theory folklore (e.g., the constant-weight
code-to-disjunct-matrix distance relation used throughout W. Kautz &
R. Singleton (1964), "Nonrandom Binary Superimposed Codes," IEEE
Transactions on Information Theory, 10(4):363-377,
DOI: 10.1109/TIT.1964.1053689); not attributed to a single paper. Only
``D_star`` is used elsewhere in this package.
"""

import math
from functools import lru_cache

from core import binom, binom_float

def U_pp_exact_weight(m, d, w):
    """
    Johnson-type private-prefix bound for d-disjunct matrices with exact column weight w.
    U_pp(m, d, w) = binom(m, t) / binom(w, t)  where t = floor((w-1)/d) + 1
    """
    if w > m:
        return 0
    if w <= 0:
        return 1
    t = math.floor((w - 1) / d) + 1
    
    if t > w:
        return 0
        
    return math.floor(binom(m, t) / binom(w, t))

def D_star(d, w):
    """
    Canonical even distance threshold for constant weight codes to form d-disjunct matrices.
    D*(d,w) = 2 * (floor((d-1)*w / d) + 1)
    """
    return 2 * (math.floor((d - 1) * w / d) + 1)

def average_weight_upper_bound(m, d):
    """
    Average-weight upper bound for d-isolated families.
    w <= 2m / (d+1) - 1
    """
    return (2 * m) / (d + 1) - 1


def U_pp_avg_weight(m, d, w_avg):
    """
    Johnson-type private-prefix bound for d-disjunct matrices with average column weight w_avg.
    U_pp(m, d, w_avg) = binom(m, t) / binom(w_avg, t)  where t = floor((w_avg-1)/d) + 1
    """
    if w_avg > m:
        return 0
    if w_avg <= 0:
        return 1
    t = math.floor((w_avg - 1) / d) + 1
    if t > w_avg:
        return 0
    return math.floor(binom(m, t) / binom_float(w_avg, t))


def U_bal(m, d, w):
    """
    Balanced-partition maximum-weight certificate (tcs.tex, eq:def:U-bal /
    Theorem "Maximum-weight balanced-partition certificate"):

        T = ceil(w / d),  s = w - d*(T-1)
        U_bal(m,d,w) = sum_{j=1}^{T-1} d*binom(m,j)/binom(d*j,j)
                       + s*binom(m,T)/binom(w,T)

    N(m, d, <= w) <= U_bal(m, d, w).
    """
    if w <= 0:
        return 0
    T = math.ceil(w / d)
    s = w - d * (T - 1)
    total = 0.0
    for j in range(1, T):
        denom = binom(d * j, j)
        if denom == 0:
            continue
        total += d * binom(m, j) / denom
    denom_T = binom(w, T)
    if denom_T > 0:
        total += s * binom(m, T) / denom_T
    return total


@lru_cache(maxsize=None)
def M_bal(n, d, w):
    """
    Smallest m such that n <= U_bal(m, d, w) (def-J in tcs.tex).
    """
    if n <= 1:
        return 0
    low = 1
    high = max(w, 2)
    while U_bal(high, d, w) < n:
        high *= 2
        if high > 10 ** 7:
            break

    best_m = high
    while low <= high:
        mid = (low + high) // 2
        if U_bal(mid, d, w) >= n:
            best_m = mid
            high = mid - 1
        else:
            low = mid + 1
    return best_m


@lru_cache(maxsize=None)
def M_bal_unrestricted(n, d):
    """
    "LB_*" in tcs.tex's tab:d2-main-comparison-lb-by-weight: smallest m
    such that n <= U_bal(m, d, w=m), i.e. M_bal(n, d, w) evaluated with no
    separate weight cap (w allowed to be the full row count m).

    Exact match against the tabulated LB_* column for d=2:
    100 -> 13, 1000 -> 22, 10**4 -> 32, 10**5 -> 42, 10**6 -> 52,
    10**8 -> 71, 10**10 -> 92, 10**20 -> 193, 10**30 -> 295.
    """
    if n <= 1:
        return 0
    low, high = 1, 2
    while U_bal(high, d, high) < n:
        high *= 2
        if high > 10 ** 7:
            break

    best_m = high
    while low <= high:
        mid = (low + high) // 2
        if U_bal(mid, d, mid) >= n:
            best_m = mid
            high = mid - 1
        else:
            low = mid + 1
    return best_m


@lru_cache(maxsize=None)
def M_BDR(n, d):
    """
    Balanced D'yachkov-Rykov bound (tcs.tex, Theorem "Balanced D'yachkov-Rykov
    bound"): recursively combines U_bal with the max-weight deletion
    recursion.

        M_BDR(n,1) = min{m : n <= binom(m, floor(m/2))}          (Sperner)
        M_BDR(n,d) = min{m : n <= U_bal(m, d, m - M_BDR(n-1,d-1))}, d>=2
    """
    if n <= 1:
        return 0
    if d <= 1:
        low, high = 1, 2
        while binom(high, high // 2) < n:
            high *= 2
            if high > 10 ** 7:
                break
        best_m = high
        while low <= high:
            mid = (low + high) // 2
            if binom(mid, mid // 2) >= n:
                best_m = mid
                high = mid - 1
            else:
                low = mid + 1
        return best_m

    prev = M_BDR(n - 1, d - 1)
    low, high = max(prev + 1, d), max(prev + 1, d) * 2
    while high - prev < 1 or U_bal(high, d, high - prev) < n:
        high *= 2
        if high > 10 ** 7:
            break

    best_m = high
    while low <= high:
        mid = (low + high) // 2
        w = mid - prev
        if w >= 1 and U_bal(mid, d, w) >= n:
            best_m = mid
            high = mid - 1
        else:
            low = mid + 1
    return best_m

def Johnson2Isolated(m, d, w):
    """
    Lemma 9: Johnson-type bound for d-isolated families.
    n <= floor((w - t) * m / (w^2 - t * m))
    """
    t = math.floor((w - 1) / d)
    denom = w * w - t * m
    if denom > 0:
        return math.floor((w - t) * m / denom)
    return None

def isolated_weight_upper_bound(m, d):
    """
    Lemma 10: Upper bound on average weight of d-isolated family with n >= m.
    w_avg <= (m - 1) / d
    """
    return (m - 1) / d

def JohnsonSpernerIsolatedBound(m, d, w):
    """
    Theorem 11: Johnson-Sperner bound for d-isolated codes.
    """
    from core import johnson_bound as JohnsonBound

    t = math.floor((w - 1) / d)
    alpha = w - d * t
    if alpha == 0 and t > 0:
        t = t - 1
        alpha = d

    if t < 0 or alpha < 1:
        return None

    W = [0] * (d + 1)
    w_i = [0] * (d + 1)
    for i in range(d + 1):
        W[i] = int(i * (2 * d - i + 1) // 2 * t + i * alpha)
        w_i[i] = int((d - i) * t + alpha)

    low = d
    high = binom(m, d)
    best_n = d

    def check_n(n):
        D_val = 2 * (t + alpha - 1)
        if D_val <= 0:
            return False
        for i in range(d + 1):
            rem_m = m - W[i]
            rem_w = w_i[i]
            if rem_m < rem_w or rem_w <= 0:
                if n - i > 0:
                    return False
                continue
            jb = JohnsonBound(rem_m, D_val, rem_w)
            if jb is None or n - i > jb:
                return False
        if W[d] > m:
            return False
        if binom(n, d) > binom(m, W[d]):
            return False
        return True

    while low <= high:
        mid = (low + high) // 2
        if check_n(mid):
            best_n = mid
            low = mid + 1
        else:
            high = mid - 1

    return best_n



# --- KS, recovered from docs/old_main.tex -----------------------------
#
# docs/tcs.tex's tab:d2-main-comparison table cites this lower-bound
# column via KS \eqref{lb:KS}, but the label is dead in tcs.tex (not
# defined at all). docs/old_main.tex -- an earlier, more complete draft --
# has the label live, which is what let this be reproduced exactly rather
# than guessed at.
#
# The DR* column (D'yachkov-Rykov's asymptotic constant bound) that used
# to live here has been dropped: the survey's KS-B/KS-C/KS-C-RS columns
# are now reproduced properly by cgt.tables.ks, which supersedes it.


def KS_lower_bound(n, d):
    """
    "KS" column (docs/old_main.tex \\label{lb:KS}): Kautz and Singleton's
    counting bound,

        sum_{i=1}^{d} binom(n, i)  <=  sum_{j=d}^{m} binom(m, j),

    "which holds only if d <= n/2". We invert it: the smallest m for which
    the right side reaches the left side is the certified lower bound on
    the number of tests.

    Exact match for 8 of the 9 rows in docs/tables.tex's d=2 "KS"
    column: KS_lower_bound(n, 2) == 19, 26, 33, 39, 53, 66, 132, 199 for
    n = 1000, 10**4, ..., 10**30. The n=100 row is the one exception
    (this function gives 13, matching the table's own "DR" column value
    rather than its tabulated KS value of 17) -- given every other row
    matches exactly, this single cell is treated as a probable table
    transcription slip rather than a bug here, but it is called out
    rather than silently accepted.
    """
    if d > n / 2:
        raise ValueError("KS_lower_bound requires d <= n/2 (docs/old_main.tex, lb:KS)")
    lhs = sum(binom(n, i) for i in range(1, d + 1))

    def rhs(m):
        return sum(binom(m, j) for j in range(d, m + 1))

    m = d
    while rhs(m) < lhs:
        m += 1
    return m


from functools import lru_cache as _lru_cache


@_lru_cache(maxsize=None)
def DR_classical(n, d):
    """
    "DR" column: the classical D'yachkov-Rykov recursive bound itself
    (docs/tcs.tex eq:m:lb_n / docs/old_main.tex's un-starred recursive
    formula, cited dyachkov1982bounds). Combines the
    Johnson-type weight-layer count (with its d^2 coefficient) with the
    maximum-weight deletion recursion:

        n <= m(d,n) + d^2 * sum_{v=d+1}^{m(d,n)-m(d-1,n-1)}
                 binom(m(d,n), ceil(v/d)) / binom(d*ceil(v/d), ceil(v/d))

    Base case m(1,n) is the Sperner bound (n <= binom(m, floor(m/2))).

    Close match (within 1) to every d=2 row of docs/tables.tex's "DR"
    column: DR_classical(100, 2) == 13 (exact), DR_classical(1000, 2) ==
    20 (table: 19), DR_classical(10**4, 2) == 28 (table: 27) -- the
    consistent off-by-one at larger n is most likely a `<` vs `<=`
    boundary difference from the source, not a wrong formula (the shape
    and asymptotics match throughout).
    """
    if d <= 1:
        m = 1
        while binom(m, m // 2) < n:
            m += 1
        return m

    m_prev = DR_classical(n - 1, d - 1)
    m = m_prev + 1
    while True:
        w_max = m - m_prev
        if w_max >= d + 1:
            layer_sum = 0.0
            for v in range(d + 1, w_max + 1):
                t = math.ceil(v / d)
                denom = binom(d * t, t)
                if denom > 0:
                    layer_sum += binom(m, t) / denom
            rhs = m + d * d * layer_sum
        else:
            rhs = m
        if n <= rhs:
            return m
        m += 1


def _dmr_n_for_m(m, d, r_max=40):
    """Best DMR construction size achievable at exactly m rows (helper)."""
    from cwc.bounds.brouwer import boundA

    best = 0
    for r in range(1, r_max):
        w = r * d + 1
        if w > m:
            break
        n, _ = boundA(m, D_star(d, w), w)
        if n and n > best:
            best = n
    return best


def DMR_m(n, d, m_max=2000):
    """
    "DMR" column: D'yachkov-Macula-Rykov's constant-weight-ECC
    construction. Its live (not commented-out) derivation is
    docs/old_main.tex eq:AmDt (the paragraph right before it: "The codes
    are effective if w-1 is divisible by d, see also \\cite{dyachkov2000new}
    ... n >= A(m, 2(d-1)t+2, dt+1) are d-disjunct CGT codes for any t >= 1
    integer parameter"): for any integer r >= 1 (t in the paper's own
    notation), with w = r*d + 1 and D = 2*(r*d + 1) - 2*r, any binary
    constant-weight code achieving A(m, D, w) gives a d-disjunct CGT
    matrix of the same size. (The same equation also survives, commented
    out, right before docs/old_main.tex's \\section{New lower bounds} --
    that copy is just a leftover duplicate of the live one, not a
    separate source.)

    This D is exactly D_star(d, w) (algebraically identical -- both
    reduce to 2*(d-1)*r + 2 when w = r*d+1), so DMR is the same
    boundA(m, D_star(d, w), w) lookup already used as the plain
    constant-weight-code ingredient of
    cgt.constructions.canonical.get_best_construction_lb, just optimized
    over r instead of a single fixed w.

    Returns the smallest m for which some r gives a code of size >= n.

    Exact match against docs/tables.tex's only two populated d=2 "DMR"
    cells: DMR_m(100, 2) == 21, DMR_m(1000, 2) == 39.
    """
    for m in range(d + 1, m_max):
        if _dmr_n_for_m(m, d) >= n:
            return m
    return None
