"""Explicit constructions from Kautz and Singleton.

Reference: W. Kautz & R. Singleton (1964). "Nonrandom Binary
Superimposed Codes." IEEE Transactions on Information Theory,
10(4):363-377. DOI: 10.1109/TIT.1964.1053689.

The public convention in this module is the group-testing convention:
``m`` is the number of tests, ``N`` the number of items/columns, and ``d``
the maximum number of defectives.  A codeword is stored as the set of row
indices containing a one.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations, product
from typing import Any, Iterable, Iterator, Sequence


@dataclass
class SuperimposedCode:
    columns: list[frozenset[int]]
    m: int
    d: int
    guarantee: str
    construction: str
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.columns = [frozenset(column) for column in self.columns]
        if any(not column.issubset(range(self.m)) for column in self.columns):
            raise ValueError("A column contains a row outside 0, ..., m-1")
        if self.guarantee not in {"disjunct", "separable"}:
            raise ValueError("guarantee must be 'disjunct' or 'separable'")

    @property
    def N(self) -> int:
        return len(self.columns)

    @property
    def weights(self) -> tuple[int, ...]:
        return tuple(len(column) for column in self.columns)

    @property
    def w(self) -> int | None:
        weights = set(self.weights)
        return next(iter(weights)) if len(weights) == 1 else None

    @property
    def historical_name(self) -> str:
        return f"{'ZFD' if self.guarantee == 'disjunct' else 'UD'}_{self.d}"

    def binary_columns(self) -> list[list[int]]:
        return [[int(row in column) for row in range(self.m)] for column in self.columns]

    def binary_matrix(self) -> list[list[int]]:
        return [[int(row in column) for column in self.columns] for row in range(self.m)]


def _unions(columns: Sequence[frozenset[int]], d: int) -> Iterator[tuple[tuple[int, ...], frozenset[int]]]:
    yield (), frozenset()
    for size in range(1, min(d, len(columns)) + 1):
        for indices in combinations(range(len(columns)), size):
            yield indices, frozenset().union(*(columns[index] for index in indices))


def is_d_disjunct(code: SuperimposedCode, d: int | None = None) -> bool:
    order = code.d if d is None else d
    for excluded, union in _unions(code.columns, order):
        excluded_set = set(excluded)
        for index, column in enumerate(code.columns):
            if index not in excluded_set and column <= union:
                return False
    return True


def is_d_separable(code: SuperimposedCode, d: int | None = None) -> bool:
    order = code.d if d is None else d
    seen: dict[frozenset[int], tuple[int, ...]] = {}
    for indices, union in _unions(code.columns, order):
        if union in seen and seen[union] != indices:
            return False
        seen[union] = indices
    return True


def identity_zfd(m: int, d: int | None = None) -> SuperimposedCode:
    if m < 1:
        raise ValueError("m must be positive")
    order = m - 1 if d is None else d
    if not 0 <= order < m:
        raise ValueError("identity_zfd requires 0 <= d < m")
    return SuperimposedCode(
        [{row} for row in range(m)], m, order, "disjunct", "ks64-identity"
    )


def hamming_weight3_zfd(extension_degree: int) -> SuperimposedCode:
    """The weight-three layer of the binary Hamming code (KS64 V-A)."""
    if extension_degree < 2:
        raise ValueError("extension_degree must be at least 2")
    m = 2**extension_degree - 1
    triples: set[frozenset[int]] = set()
    for a in range(1, m + 1):
        for b in range(a + 1, m + 1):
            c = a ^ b
            if c and c != a and c != b:
                triples.add(frozenset((a - 1, b - 1, c - 1)))
    return SuperimposedCode(
        sorted(triples, key=lambda x: tuple(sorted(x))),
        m,
        2,
        "disjunct",
        "ks64-hamming-weight3",
        {"extension_degree": extension_degree},
    )


def constant_weight_layer_zfd(
    m: int, words: Iterable[Iterable[int]], weight: int, d: int
) -> SuperimposedCode:
    """Select and certify a constant-weight layer of a binary code (KS64 V-A)."""
    columns = [frozenset(word) for word in words if len(frozenset(word)) == weight]
    code = SuperimposedCode(
        columns, m, d, "disjunct", "ks64-constant-weight-layer", {"weight": weight}
    )
    if not is_d_disjunct(code):
        raise ValueError("the selected constant-weight layer is not d-disjunct")
    return code


def _qary_mds_words(q: int, dimension: int, length: int) -> list[tuple[int, ...]]:
    import galois

    if dimension < 1 or length < dimension:
        raise ValueError("MDS length must be at least its dimension")
    if length > q:
        raise ValueError("This implementation supports MDS length <= q")
    field = galois.GF(q)
    points = list(field.elements[:length])
    words: list[tuple[int, ...]] = []
    for raw_coefficients in product(range(q), repeat=dimension):
        coefficients = field(raw_coefficients)
        values = []
        for point in points:
            value = field(0)
            power = field(1)
            for coefficient in coefficients:
                value += coefficient * power
                power *= point
            values.append(int(value))
        words.append(tuple(values))
    return words


def _concatenate_symbols(
    words: Iterable[Sequence[int]], inner_columns: Sequence[frozenset[int]], inner_m: int
) -> list[frozenset[int]]:
    result = []
    for word in words:
        result.append(
            frozenset(
                block * inner_m + row
                for block, symbol in enumerate(word)
                for row in inner_columns[symbol]
            )
        )
    return result


def qary_mds_zfd(q: int, dimension: int, d: int) -> SuperimposedCode:
    """One-hot q-ary MDS/Reed-Solomon Kautz-Singleton construction."""
    if d < 1 or dimension < 1:
        raise ValueError("d and dimension must be positive")
    outer_length = 1 + d * (dimension - 1)
    words = _qary_mds_words(q, dimension, outer_length)
    inner = [frozenset((symbol,)) for symbol in range(q)]
    columns = _concatenate_symbols(words, inner, q)
    return SuperimposedCode(
        columns,
        q * outer_length,
        d,
        "disjunct",
        "ks64-qary-mds",
        {"q": q, "dimension": dimension, "outer_length": outer_length},
    )


def latin_square_zfd(q: int) -> SuperimposedCode:
    """The k=2, d=2 construction from one cyclic Latin square."""
    if q < 3:
        raise ValueError("q must be at least 3")
    columns = [
        frozenset((row, q + column, 2 * q + ((row + column) % q)))
        for row in range(q)
        for column in range(q)
    ]
    return SuperimposedCode(
        columns, 3 * q, 2, "disjunct", "ks64-latin-square", {"q": q}
    )


def iterated_qary_composition(q: int, dimension: int, d: int, levels: int) -> SuperimposedCode:
    """Repeated KS64 V-C composition, starting from the q-word identity code."""
    if levels < 1:
        raise ValueError("levels must be positive")
    inner = identity_zfd(q, min(d, q - 1))
    for level in range(1, levels + 1):
        outer_length = 1 + d * (dimension - 1)
        words = _qary_mds_words(inner.N, dimension, outer_length)
        columns = _concatenate_symbols(words, inner.columns, inner.m)
        inner = SuperimposedCode(
            columns,
            inner.m * outer_length,
            d,
            "disjunct",
            "ks64-iterated-qary-composition",
            {"q": q, "dimension": dimension, "levels": level},
        )
    return inner


def bibd_zfd(v: int, blocks: Iterable[Iterable[int]]) -> SuperimposedCode:
    block_list = [frozenset(block) for block in blocks]
    if not block_list:
        raise ValueError("at least one block is required")
    if any(not block or not block.issubset(range(v)) for block in block_list):
        raise ValueError("blocks must be nonempty subsets of 0, ..., v-1")
    weights = {len(block) for block in block_list}
    if len(weights) != 1:
        raise ValueError("all BIBD blocks must have equal size")
    max_intersection = max(
        (len(left & right) for left, right in combinations(block_list, 2)), default=0
    )
    if max_intersection == 0:
        order = next(iter(weights))
    else:
        order = (next(iter(weights)) - 1) // max_intersection
    return SuperimposedCode(
        block_list,
        v,
        order,
        "disjunct",
        "ks64-block-design",
        {"max_block_intersection": max_intersection},
    )


def projective_plane_zfd(q: int) -> SuperimposedCode:
    """Projective-plane BIBD: v=q^2+q+1, w=q+1, d=q."""
    import galois

    field = galois.GF(q)

    def normalized_vectors() -> list[tuple[int, int, int]]:
        representatives: set[tuple[int, int, int]] = set()
        for raw in product(range(q), repeat=3):
            vector = field(raw)
            if all(int(value) == 0 for value in vector):
                continue
            first = next(value for value in vector if int(value) != 0)
            normalized = tuple(int(value / first) for value in vector)
            representatives.add(normalized)
        return sorted(representatives)

    points = normalized_vectors()
    blocks = []
    for line in normalized_vectors():
        line_vector = field(line)
        blocks.append(
            frozenset(
                index
                for index, point in enumerate(points)
                if int(sum(field(point) * line_vector, field(0))) == 0
            )
        )
    code = bibd_zfd(len(points), blocks)
    code.construction = "ks64-projective-plane"
    code.metadata["q"] = q
    return code


def bch_parity_ud2(extension_degree: int) -> SuperimposedCode:
    rows = _primitive_bch_rows(extension_degree, (1, 3))
    columns = []
    for row in rows:
        columns.append(
            frozenset(2 * index + (0 if bit == 0 else 1) for index, bit in enumerate(row))
        )
    return SuperimposedCode(
        columns, 4 * extension_degree, 2, "separable", "ks64-bch-parity-ud2",
        {"extension_degree": extension_degree},
    )


def bch_parity_ud3(extension_degree: int) -> SuperimposedCode:
    rows = _primitive_bch_rows(extension_degree, (1, 3, 5))
    pairs = list(combinations(range(3 * extension_degree), 2))
    columns = []
    for row in rows:
        columns.append(
            frozenset(4 * index + 2 * row[left] + row[right] for index, (left, right) in enumerate(pairs))
        )
    return SuperimposedCode(
        columns, 4 * len(pairs), 3, "separable", "ks64-bch-parity-ud3",
        {"extension_degree": extension_degree},
    )


def _primitive_bch_rows(extension_degree: int, odd_powers: Sequence[int]) -> list[tuple[int, ...]]:
    import galois

    if extension_degree < 2:
        raise ValueError("extension_degree must be at least 2")
    field = galois.GF(2**extension_degree)
    alpha = field.primitive_element
    rows = []
    for exponent in range(2**extension_degree - 1):
        rows.append(
            tuple(
                int(bit)
                for power in odd_powers
                for bit in (alpha ** (power * exponent)).vector()
            )
        )
    return rows


def graph_ud2(graph: Any) -> SuperimposedCode:
    """Weight-two UD2 code from a simple graph of girth at least five."""
    import networkx as nx

    simple = nx.Graph(graph)
    if simple.number_of_edges() != graph.number_of_edges() or nx.number_of_selfloops(simple):
        raise ValueError("the graph must be simple")
    try:
        girth = nx.girth(simple)
    except nx.NetworkXError:
        girth = float("inf")
    if girth < 5:
        raise ValueError("the graph must have girth at least five")
    nodes = list(simple.nodes())
    index = {node: position for position, node in enumerate(nodes)}
    columns = [frozenset((index[left], index[right])) for left, right in simple.edges()]
    return SuperimposedCode(
        columns, len(nodes), 2, "separable", "ks64-graph-ud2", {"girth": girth}
    )


def moore_graph_ud2(degree: int) -> SuperimposedCode:
    import networkx as nx

    if degree == 2:
        graph = nx.cycle_graph(5)
    elif degree == 3:
        graph = nx.petersen_graph()
    elif degree == 7:
        graph = nx.hoffman_singleton_graph()
    else:
        raise ValueError("implemented Moore degrees are 2, 3, and 7")
    code = graph_ud2(graph)
    code.construction = "ks64-moore-graph-ud2"
    code.metadata["degree"] = degree
    return code


def split_field_bibd_ud2(v: int, blocks: Iterable[Iterable[int]]) -> SuperimposedCode:
    """Weight-two split-field UD2 construction from a lambda=1 BIBD."""
    block_list = [frozenset(block) for block in blocks]
    if any(len(left & right) > 1 for left, right in combinations(block_list, 2)):
        raise ValueError("block intersections must have size at most one")
    columns = [
        frozenset((point, v + block_index))
        for block_index, block in enumerate(block_list)
        for point in block
    ]
    return SuperimposedCode(
        columns, v + len(block_list), 2, "separable", "ks64-split-field-bibd-ud2"
    )


def split_projective_plane_ud2(q: int) -> SuperimposedCode:
    design = projective_plane_zfd(q)
    code = split_field_bibd_ud2(design.m, design.columns)
    code.construction = "ks64-split-projective-plane-ud2"
    code.metadata["q"] = q
    return code


def pairwise_ud2_composition(base: SuperimposedCode) -> SuperimposedCode:
    """Correct general pairwise composition using a conservative connector field.

    KS64 uses a shorter recursively constructed connector for its special seed.
    Assigning one connector coordinate to each ordered pair gives the same
    composition principle for any UD2 input and keeps the implementation fully
    explicit and verifiable.
    """
    if base.guarantee != "separable" or base.d < 2:
        raise ValueError("base must be at least 2-separable")
    connector_offset = 2 * base.m
    columns = []
    for left in range(base.N):
        for right in range(base.N):
            columns.append(
                frozenset(base.columns[left])
                | frozenset(base.m + row for row in base.columns[right])
                | frozenset((connector_offset + left * base.N + right,))
            )
    return SuperimposedCode(
        columns,
        2 * base.m + base.N**2,
        2,
        "separable",
        "ks64-pairwise-composition",
        {"base_construction": base.construction},
    )


def ks64_ud2_seed() -> SuperimposedCode:
    """The paper's 7-test, 9-item first pairwise-composition example."""
    columns = []
    for left in range(3):
        for right in range(3):
            columns.append(frozenset((left, 3 + right, 6)) if left == right else frozenset((left, 3 + right)))
    return SuperimposedCode(columns, 7, 2, "separable", "ks64-pairwise-seed-ud2")
