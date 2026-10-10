"""Independent scalar oracle for global distance, including wildcard semantics."""

import random
from itertools import product

import pytest

from carmack.barcode.barcode_utils import edit_distance


def scalar_distance(first, second, marker="N", wildcard=True):
    """Row-wise Wagner–Fischer recurrence, independent of the bit-vector code."""
    previous = list(range(len(second) + 1))
    for i, left in enumerate(first, 1):
        current = [i]
        for j, right in enumerate(second, 1):
            equal = left == right or (wildcard and (left == marker or right == marker))
            current.append(min(previous[j] + 1, current[j - 1] + 1, previous[j - 1] + (not equal)))
        previous = current
    return previous[-1]


@pytest.mark.parametrize("wildcard", [False, True])
def test_exhaustive_short_sequences_match_independent_oracle(wildcard):
    sequences = ["".join(chars) for size in range(5) for chars in product("ACN", repeat=size)]
    for first in sequences:
        for second in sequences:
            assert edit_distance(first, second, n_matches_any=wildcard) == scalar_distance(
                first, second, wildcard=wildcard
            ), (first, second, wildcard)


@pytest.mark.parametrize("marker", ["N", "?", "é", "", "two"])
@pytest.mark.parametrize("wildcard", [False, True])
def test_arbitrary_width_and_symbol_semantics(marker, wildcard):
    rng = random.Random(20261010)
    lengths = [0, 1, 9, 10, 16, 31, 32, 63, 64, 65, 127, 128, 129]
    for _ in range(50):
        first = "".join(rng.choices("ACGTN?é", k=rng.choice(lengths)))
        second = "".join(rng.choices("ACGTN?é", k=rng.choice(lengths)))
        expected = scalar_distance(first, second, marker, wildcard)
        assert edit_distance(first, second, marker, wildcard) == expected
        assert edit_distance(second, first, marker, wildcard) == expected


@pytest.mark.parametrize(
    "first,second,expected",
    [
        ("", "ACT", 3),
        ("AC", "TTACGG", 4),
        ("ACGT", "N", 3),
        ("NNNN", "ACGT", 0),
        ("ACGT", "ACNT", 0),
        ("kitten", "sitting", 3),
    ],
)
def test_global_boundaries_and_wildcards(first, second, expected):
    assert edit_distance(first, second) == expected
