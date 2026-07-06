"""
Rotation-aware, error-tolerant containment between circular repeat patterns,
plus helpers to (a) annotate which patterns are variations of higher-ranked
patterns and (b) select the top-N patterns that are NOT variations of an
already-selected pattern.

Two patterns are treated as the same "family" when, viewing both as circular
sequences, the shorter fits inside the other (any rotation) allowing
sequencing-error mismatches/indels, at >= threshold identity.

Pure standard library — patterns are few and short, so the O(n*m) fuzzy-substring
DP is plenty fast. Shared by main.py's population finder and the standalone
annotate_pattern_variations.py script so the logic stays in one place.
"""

from typing import Dict, List, Optional, Sequence, Tuple

DEFAULT_THRESHOLD = 0.93


# ── core matching ─────────────────────────────────────────────────────────────

def infix_edit_distance(s: str, t: str) -> int:
    """Minimum edit distance of `s` against ANY substring of `t`.

    Free start/end positions in `t` (classic fuzzy substring search): the first
    DP row is all zeros so `s` may begin matching anywhere in `t`, and the answer
    is the minimum over the final row so it may end anywhere.
    """
    m, n = len(s), len(t)
    prev = [0] * (n + 1)
    for i in range(1, m + 1):
        cur = [i] + [0] * n
        si = s[i - 1]
        for j in range(1, n + 1):
            cost = 0 if si == t[j - 1] else 1
            cur[j] = min(prev[j - 1] + cost,  # substitute / match
                         prev[j] + 1,         # delete from s
                         cur[j - 1] + 1)      # insert into s
        prev = cur
    return min(prev)


def containment_identity(short: str, long_: str) -> float:
    """Best identity of circular `short` placed inside circular `long_`.

    Rotation is handled by searching `short` against the doubled `long_`.
    Returns identity in [0, 1] = 1 - edits / len(short). Exact containment is
    fast-pathed before falling back to the DP.
    """
    if not short or len(short) > len(long_):
        return 0.0
    if short in (long_ + long_):          # exact rotation-substring
        return 1.0
    dist = infix_edit_distance(short, long_ + long_)
    return 1.0 - dist / len(short)


def variant_identity(a: str, b: str) -> float:
    """Identity for the 'shorter fits inside longer' relationship between a and b."""
    short, long_ = (a, b) if len(a) <= len(b) else (b, a)
    return containment_identity(short, long_)


def length_ratio(a: str, b: str) -> float:
    """min/max of the two lengths, in (0, 1]. 1.0 = identical length."""
    la, lb = len(a), len(b)
    if la == 0 or lb == 0:
        return 0.0
    return min(la, lb) / max(la, lb)


def is_variant(a: str, b: str, threshold: float = DEFAULT_THRESHOLD,
               min_length_ratio: float = 0.0) -> bool:
    """Same circle if one fits inside the other (>= threshold identity) AND the
    two lengths are within min_length_ratio. With min_length_ratio=0 the length
    guard is off (pure containment). Raise it (e.g. 0.75) to keep a short
    sub-repeat as its own family instead of folding it into a longer pattern."""
    if length_ratio(a, b) < min_length_ratio:
        return False
    return variant_identity(a, b) >= threshold


# ── annotation: mark each pattern's highest-ranked parent ──────────────────────

def assign_variations(seqs: Sequence[str],
                      threshold: float = DEFAULT_THRESHOLD
                      ) -> List[Optional[Tuple[int, float]]]:
    """For each pattern (given best-first), return (parent_index, identity) of the
    HIGHEST-RANKED earlier pattern it is a direct variation of, or None.

    Index-based so callers can map back to ranks/records however they like.
    """
    result: List[Optional[Tuple[int, float]]] = []
    for idx, seq in enumerate(seqs):
        match: Optional[Tuple[int, float]] = None
        for jdx in range(idx):  # earlier == higher-ranked
            ident = variant_identity(seq, seqs[jdx])
            if ident >= threshold:
                match = (jdx, ident)
                break  # first hit is the highest-ranked one
        result.append(match)
    return result


# ── selection: top-N patterns that are not variations of a kept pattern ────────

def select_non_variant(seqs: Sequence[str],
                       threshold: float = DEFAULT_THRESHOLD,
                       limit: Optional[int] = None,
                       min_length_ratio: float = 0.0
                       ) -> Tuple[List[int], Dict[int, int]]:
    """Greedily keep patterns that are not a variation of an already-kept one.

    Walks `seqs` in order (assumed best-first). A pattern is kept as a new
    "primary" if it is not a variant of any previously-kept primary; otherwise it
    is folded into the first such primary.

    Returns (primary_indices, assignment) where:
      - primary_indices: indices kept as primaries, in rank order
      - assignment: maps EVERY index -> the primary index it belongs to
                    (a primary maps to itself)

    The full list is always scanned so variant counts are complete; `limit` only
    bounds how many primaries the caller will typically display, but selection
    continues so later patterns can still be attributed to a primary.
    """
    primaries: List[int] = []
    assignment: Dict[int, int] = {}
    for idx, seq in enumerate(seqs):
        parent = None
        for pidx in primaries:
            if is_variant(seq, seqs[pidx], threshold, min_length_ratio):
                parent = pidx
                break
        if parent is None:
            primaries.append(idx)
            assignment[idx] = idx
        else:
            assignment[idx] = parent
    if limit is not None:
        primaries = primaries[:limit]
    return primaries, assignment


def variant_counts(primary_indices: Sequence[int],
                   assignment: Dict[int, int]) -> Dict[int, int]:
    """Number of variant patterns folded into each primary index."""
    counts = {p: 0 for p in primary_indices}
    primary_set = set(primary_indices)
    for idx, parent in assignment.items():
        if idx != parent and parent in primary_set:
            counts[parent] += 1
    return counts
