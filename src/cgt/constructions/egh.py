"""Eppstein--Goodrich--Hirschberg Chinese Remainder Sieve.

Reference: D. Eppstein, M. T. Goodrich & D. S. Hirschberg (2007).
"Improved Combinatorial Group Testing Algorithms for Real-World Problem
Sizes." SIAM Journal on Computing, 36(5):1360-1375.
DOI: 10.1137/050631847.

The construction uses one residue block for every selected prime power.  Item
``x`` has a one in row ``x mod modulus`` of each block.  The prime-power search
follows Figures 1--2 of the EGH paper, using an equivalent Pareto-pruned
dynamic program and Python 3 integer arithmetic.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from math import prod

from cgt.constructions.ks64 import SuperimposedCode


def primes() -> Iterator[int]:
    """Yield the primes in increasing order using an incremental sieve."""
    composites: dict[int, list[int]] = {}
    candidate = 2
    while True:
        witnesses = composites.pop(candidate, None)
        if witnesses is None:
            yield candidate
            composites[candidate * candidate] = [candidate]
        else:
            for prime in witnesses:
                composites.setdefault(candidate + prime, []).append(prime)
        candidate += 1


def initial_primes(target: int) -> tuple[int, ...]:
    """Smallest initial sequence of primes whose product reaches ``target``."""
    if target < 1:
        raise ValueError("target must be positive")
    selected: list[int] = []
    product = 1
    for prime in primes():
        selected.append(prime)
        product *= prime
        if product >= target:
            return tuple(selected)
    raise AssertionError("the prime generator is infinite")


def backtrack_exponents(
    prime_list: tuple[int, ...], max_power: int, target: int
) -> tuple[int, tuple[int, ...]] | None:
    """Minimize the sum of selected prime powers, as in EGH Figure 1.

    Every selected power is at most ``max_power``.  The returned first value is
    the sum of the non-unit powers and the second is the exponent vector.
    """
    if max_power < 2 or target < 1:
        raise ValueError("max_power must be at least 2 and target positive")

    # A direct transcription of the recursive pseudocode becomes slow for the
    # 10^20 and 10^30 survey rows.  This is the equivalent multiplicative
    # knapsack DP.  After each prime, discard a state whenever another state
    # has at least its product at no greater sum.
    states: dict[int, tuple[int, tuple[int, ...]]] = {1: (0, ())}
    for prime in prime_list:
        choices = [(0, 1)]
        exponent, power = 1, prime
        while power <= max_power:
            choices.append((exponent, power))
            exponent += 1
            power *= prime

        candidates: dict[int, tuple[int, tuple[int, ...]]] = {}
        for old_product, (old_sum, old_exponents) in states.items():
            for exponent, power in choices:
                new_product = min(target, old_product * power)
                candidate = (
                    old_sum + (power if exponent else 0),
                    old_exponents + (exponent,),
                )
                previous = candidates.get(new_product)
                if previous is None or candidate < previous:
                    candidates[new_product] = candidate

        states = {}
        best_sum: int | None = None
        for product_value in sorted(candidates, reverse=True):
            candidate = candidates[product_value]
            if best_sum is None or candidate[0] < best_sum:
                states[product_value] = candidate
                best_sum = candidate[0]

    return states.get(target)


@dataclass(frozen=True)
class EGHParameters:
    n: int
    d: int
    variant: str
    primes: tuple[int, ...]
    exponents: tuple[int, ...]
    moduli: tuple[int, ...]
    target: int
    modulus_product: int
    m: int
    w: int
    maximum_pair_overlap: int
    minimum_distance: int

    def as_dict(self) -> dict[str, object]:
        return {
            "n": self.n,
            "d": self.d,
            "variant": self.variant,
            "primes": list(self.primes),
            "exponents": list(self.exponents),
            "moduli": list(self.moduli),
            "target": self.target,
            "modulus_product": self.modulus_product,
            "m": self.m,
            "w": self.w,
            "maximum_pair_overlap": self.maximum_pair_overlap,
            "minimum_distance": self.minimum_distance,
        }


def _maximum_pair_overlap(n: int, moduli: tuple[int, ...]) -> int:
    """Exact maximum overlap between two residue columns among ``0..n-1``."""
    product = 1
    overlap = 0
    for modulus in sorted(moduli):
        if product * modulus >= n:
            break
        product *= modulus
        overlap += 1
    return overlap


def egh_parameters(n: int, d: int, variant: str = "backtrack") -> EGHParameters:
    """Return the EGH construction parameters without materializing columns."""
    if n < 2 or d < 1:
        raise ValueError("n must be at least 2 and d must be positive")
    if variant not in {"general", "backtrack"}:
        raise ValueError("variant must be 'general' or 'backtrack'")

    target = n**d
    prime_list = initial_primes(target)
    if variant == "general":
        exponents = (1,) * len(prime_list)
    else:
        result = backtrack_exponents(prime_list, prime_list[-1], target)
        if result is None:
            raise RuntimeError("EGH backtracking unexpectedly found no solution")
        _, exponents = result

    moduli = tuple(
        prime**exponent
        for prime, exponent in zip(prime_list, exponents)
        if exponent
    )
    overlap = _maximum_pair_overlap(n, moduli)
    return EGHParameters(
        n=n,
        d=d,
        variant=variant,
        primes=prime_list,
        exponents=exponents,
        moduli=moduli,
        target=target,
        modulus_product=prod(moduli),
        m=sum(moduli),
        w=len(moduli),
        maximum_pair_overlap=overlap,
        minimum_distance=2 * (len(moduli) - overlap),
    )


def egh_crt_zfd(n: int, d: int, variant: str = "backtrack") -> SuperimposedCode:
    """Materialize the EGH binary residue matrix as column supports."""
    parameters = egh_parameters(n, d, variant)
    offsets: list[int] = []
    offset = 0
    for modulus in parameters.moduli:
        offsets.append(offset)
        offset += modulus
    columns = [
        frozenset(
            base + item % modulus
            for base, modulus in zip(offsets, parameters.moduli)
        )
        for item in range(n)
    ]
    return SuperimposedCode(
        columns,
        parameters.m,
        d,
        "disjunct",
        f"egh-crt-{variant}",
        parameters.as_dict(),
    )
