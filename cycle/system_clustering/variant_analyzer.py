"""Within-cluster variant classification for transposon systems."""

from __future__ import annotations

import logging

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist

from .config import DEFAULT_FLANKING_EDIT_THRESHOLD

log = logging.getLogger(__name__)


def _edit_distance(a: str, b: str) -> int:
    """Compute Levenshtein edit distance between two strings."""
    if not a:
        return len(b)
    if not b:
        return len(a)

    n, m = len(a), len(b)
    # Optimise: use two-row DP
    prev = list(range(m + 1))
    curr = [0] * (m + 1)
    for i in range(1, n + 1):
        curr[0] = i
        for j in range(1, m + 1):
            cost = 0 if a[i - 1] == b[j - 1] else 1
            curr[j] = min(curr[j - 1] + 1, prev[j] + 1, prev[j - 1] + cost)
        prev, curr = curr, prev
    return prev[m]


class VariantAnalyzer:
    """Classify variants within a cluster of related transposon systems.

    Level 1 (L1): Group by flanking-region similarity.
        Concatenate upstream + downstream flanking sequences, compute pairwise
        edit distances, single-linkage clustering at ``flanking_edit_threshold``.

    Level 2 (L2): Within each L1 group, subdivide by best guide hit sequence
        (exact match on ``seq_noncoding``).
    """

    def __init__(
        self, flanking_edit_threshold: int = DEFAULT_FLANKING_EDIT_THRESHOLD
    ):
        self.flanking_edit_threshold = flanking_edit_threshold

    def classify_cluster(
        self, member_ids: list[str], transposon_metadata: dict[str, dict]
    ) -> dict:
        """Classify members of a single cluster into L1 and L2 variants.

        Returns a dict with:
          - l1_assignments: {is_id: l1_label}
          - l2_assignments: {is_id: l2_label}
          - n_l1_groups: int
          - n_l2_groups: int
          - member_meta: {is_id: summary info}
        """
        n = len(member_ids)

        # Collect member metadata for output
        member_meta = {}
        for mid in member_ids:
            meta = transposon_metadata.get(mid, {})
            guide_hits = meta.get("guide_hits", [])
            best_len = max((h["length"] for h in guide_hits), default=0)
            member_meta[mid] = {
                "sample_id": meta.get("sample_id", ""),
                "organism": meta.get("organism", ""),
                "is_length": meta.get("is_length", 0),
                "n_orfs": len(meta.get("orfs", [])),
                "best_guide_length": best_len,
            }

        # ----- Level 1: flanking similarity -----
        l1_assignments: dict[str, int] = {}

        if n == 1:
            l1_assignments[member_ids[0]] = 0
        else:
            flanking_seqs = []
            for mid in member_ids:
                meta = transposon_metadata.get(mid, {})
                up = meta.get("flanking_upstream", "")
                down = meta.get("flanking_downstream", "")
                flanking_seqs.append(up + down)

            # Pairwise edit distances
            dists = np.zeros(n * (n - 1) // 2, dtype=np.float64)
            idx = 0
            for i in range(n):
                for j in range(i + 1, n):
                    dists[idx] = _edit_distance(flanking_seqs[i], flanking_seqs[j])
                    idx += 1

            if dists.max() == 0:
                # All identical flanking
                for mid in member_ids:
                    l1_assignments[mid] = 0
            else:
                Z = linkage(dists, method="single")
                labels = fcluster(Z, t=self.flanking_edit_threshold, criterion="distance")
                for mid, lab in zip(member_ids, labels):
                    l1_assignments[mid] = int(lab - 1)  # 0-indexed

        # ----- Level 2: guide sequence -----
        l2_assignments: dict[str, str] = {}
        l2_counter = 0

        # Group by L1
        l1_groups: dict[int, list[str]] = {}
        for mid, l1 in l1_assignments.items():
            l1_groups.setdefault(l1, []).append(mid)

        guide_to_l2: dict[tuple, int] = {}
        for l1_label, members in sorted(l1_groups.items()):
            for mid in members:
                meta = transposon_metadata.get(mid, {})
                guide_seq = meta.get("best_guide_seq", "")
                key = (l1_label, guide_seq)
                if key not in guide_to_l2:
                    guide_to_l2[key] = l2_counter
                    l2_counter += 1
                l2_assignments[mid] = f"{l1_label}.{guide_to_l2[key]}"

        n_l1 = len(set(l1_assignments.values()))
        n_l2 = len(set(l2_assignments.values()))

        return {
            "l1_assignments": l1_assignments,
            "l2_assignments": l2_assignments,
            "n_l1_groups": n_l1,
            "n_l2_groups": n_l2,
            "member_meta": member_meta,
        }
