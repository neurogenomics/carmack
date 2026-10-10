from collections.abc import Mapping
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from matplotlib.figure import Figure


def hamming_distance(a: np.ndarray, b: np.ndarray) -> int:
    """
    Calculate the Hamming distance between two arrays.

    Args:
        a (np.ndarray): First array (1D, same length as b).
        b (np.ndarray): Second array (1D, same length as a).

    Returns:
        int: Number of positions where a and b differ.
    """
    if len(a) != len(b):
        raise ValueError("Sequences must have equal length.")
    if a.ndim != 1 or b.ndim != 1:
        raise ValueError("Arrays must be 1D.")
    return int(np.count_nonzero(a != b))


def edit_distance(seq1, seq2, n_char="N", n_matches_any=True):
    """
    Calculate the Levenshtein (edit) distance between two sequences.

    Parameters:
    -----------
    seq1 : str
        First sequence
    seq2 : str
        Second sequence
    n_char : str, optional
        Character to treat as wildcard (default: 'N')
    n_matches_any : bool, optional
        If True, n_char matches any character with cost 0 (default: True)
        If False, n_char is treated as a regular character

    Returns:
    --------
    int or float
        The Levenshtein distance between seq1 and seq2

    Examples:
    ---------
    >>> edit_distance("ACGT", "ACNT")
    0  # N matches T with no penalty

    >>> edit_distance("ACGT", "ACNT", n_matches_any=False)
    1  # N is treated as a different character from T

    >>> edit_distance("kitten", "sitting")
    3
    """
    # Myers' bit-vector recurrence computes the same global unit-cost distance.
    # Equality masks include the declared wildcard in either sequence. Python's
    # arbitrary-width integers avoid a machine-word length limit.
    # https://doi.org/10.1145/316542.316550
    if not seq1:
        return len(seq2)
    if not seq2:
        return len(seq1)
    # Symmetric substitution rule permits the shorter pattern/mask.
    if len(seq1) > len(seq2):
        seq1, seq2 = seq2, seq1
    length = len(seq1)
    all_bits = (1 << length) - 1
    high_bit = 1 << (length - 1)
    masks = {}
    for i, symbol in enumerate(seq1):
        masks[symbol] = masks.get(symbol, 0) | (1 << i)
    wildcard = masks.get(n_char, 0) if n_matches_any else 0
    positive, negative, distance = all_bits, 0, length
    for symbol in seq2:
        equality = (
            all_bits if n_matches_any and symbol == n_char else masks.get(symbol, 0) | wildcard
        )
        vertical = equality | negative
        diagonal = (((equality & positive) + positive) ^ positive) | equality
        horizontal_positive = negative | ~(diagonal | positive)
        horizontal_negative = positive & diagonal
        if horizontal_positive & high_bit:
            distance += 1
        elif horizontal_negative & high_bit:
            distance -= 1
        horizontal_positive = (horizontal_positive << 1) | 1
        positive = ((horizontal_negative << 1) | ~(vertical | horizontal_positive)) & all_bits
        negative = horizontal_positive & vertical
    return distance


def make_barcode_rank_plot(barcode_counts: Mapping[str, int]) -> "Figure":
    """Create a barcode-rank plot from full-barcode counts.

    matplotlib is imported here rather than at module scope because it is the single heaviest
    import in the package -- around half a second -- and this is the only function in the
    module that needs it. Everything else here is small pure-sequence arithmetic that the
    matchers and the chemistry definitions depend on, and those sit on the import path of
    every command, including the ones that draw nothing.
    """
    import matplotlib.pyplot as plt

    counts = sorted((count for count in barcode_counts.values() if count > 0), reverse=True)

    fig, ax = plt.subplots(figsize=(10, 6))

    if counts:
        ranks = list(range(1, len(counts) + 1))
        ax.plot(ranks, counts, color="tab:blue")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(left=1)
    else:
        ax.text(
            0.5,
            0.5,
            "No valid barcodes",
            transform=ax.transAxes,
            ha="center",
            va="center",
        )

    ax.grid(True, which="both", ls="-", alpha=0.2)
    ax.set_xlabel("Barcode rank")
    ax.set_ylabel("Reads per barcode")
    ax.set_title("Barcode Rank Plot")
    fig.tight_layout()

    return fig
