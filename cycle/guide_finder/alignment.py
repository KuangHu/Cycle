"""Short alignment finder for between-sequence comparisons.

Adapted from RNA_guide_editor_finder ShortAlignmentFinder — only the
between-sequence methods needed for guide RNA detection.
"""

from __future__ import annotations

from typing import Optional


_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


class ShortAlignmentFinder:
    """Find short DNA alignments between two sequences.

    Parameters
    ----------
    min_length : int
        Minimum alignment length to report.
    max_mismatches : int
        Maximum number of mismatches allowed.
    check_forward : bool
        Find direct (forward-forward) matches.
    check_revcomp : bool
        Find reverse-complement matches.
    """

    def __init__(
        self,
        min_length: int = 9,
        max_mismatches: int = 1,
        check_forward: bool = True,
        check_revcomp: bool = True,
    ) -> None:
        if min_length < 1:
            raise ValueError("min_length must be at least 1")
        if max_mismatches < 0:
            raise ValueError("max_mismatches cannot be negative")
        if not check_forward and not check_revcomp:
            raise ValueError("At least one of check_forward or check_revcomp must be True")

        self.min_length = min_length
        self.max_mismatches = max_mismatches
        self.check_forward = check_forward
        self.check_revcomp = check_revcomp

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def find_alignments_between(
        self, sequence1: str, sequence2: str
    ) -> list[dict]:
        """Find all short alignments between two different sequences.

        Parameters
        ----------
        sequence1 : str
            First DNA sequence (e.g., flanking region).
        sequence2 : str
            Second DNA sequence (e.g., noncoding region).

        Returns
        -------
        list[dict]
            Each dict contains: pos1, pos2, length, seq1, seq2,
            mismatches, mismatch_positions, orientation.
        """
        sequence1 = sequence1.upper()
        sequence2 = sequence2.upper()

        seeds = self._find_seed_matches_between(sequence1, sequence2)

        extended = []
        for seed in seeds:
            ext = self._extend_match_between(sequence1, sequence2, seed)
            if ext:
                extended.append(ext)

        consolidated = self._consolidate_matches_between(extended)
        consolidated.sort(key=lambda x: (x["pos1"], x["pos2"]))
        return consolidated

    # ------------------------------------------------------------------
    # Seed finding
    # ------------------------------------------------------------------

    def _find_seed_matches_between(
        self, sequence1: str, sequence2: str
    ) -> list[dict]:
        seeds: list[dict] = []
        len1 = len(sequence1)
        len2 = len(sequence2)
        ml = self.min_length

        for i in range(len1 - ml + 1):
            kmer = sequence1[i : i + ml]

            if self.check_forward:
                for j in range(len2 - ml + 1):
                    target = sequence2[j : j + ml]
                    mm, mm_pos = self._count_mismatches(kmer, target)
                    if mm <= self.max_mismatches:
                        seeds.append(
                            {
                                "pos1": i,
                                "pos2": j,
                                "length": ml,
                                "orientation": "forward",
                                "mismatches": mm,
                                "mismatch_positions": mm_pos,
                            }
                        )

            if self.check_revcomp:
                kmer_rc = reverse_complement(kmer)
                for j in range(len2 - ml + 1):
                    target = sequence2[j : j + ml]
                    mm, mm_pos = self._count_mismatches(kmer_rc, target)
                    if mm <= self.max_mismatches:
                        seeds.append(
                            {
                                "pos1": i,
                                "pos2": j,
                                "length": ml,
                                "orientation": "reverse_complement",
                                "mismatches": mm,
                                "mismatch_positions": mm_pos,
                            }
                        )

        return seeds

    # ------------------------------------------------------------------
    # Extension
    # ------------------------------------------------------------------

    def _extend_match_between(
        self, sequence1: str, sequence2: str, seed: dict
    ) -> Optional[dict]:
        pos1 = seed["pos1"]
        pos2 = seed["pos2"]
        length = seed["length"]
        orientation = seed["orientation"]

        # Extend right
        while pos1 + length < len(sequence1) and pos2 + length < len(sequence2):
            test_s1 = sequence1[pos1 : pos1 + length + 1]
            test_s2_raw = sequence2[pos2 : pos2 + length + 1]
            test_s2 = test_s2_raw if orientation == "forward" else reverse_complement(test_s2_raw)

            mm, _ = self._count_mismatches(test_s1, test_s2)
            if mm <= self.max_mismatches:
                length += 1
            else:
                break

        # Extend left
        while pos1 > 0 and pos2 > 0:
            test_s1 = sequence1[pos1 - 1 : pos1 + length]
            test_s2_raw = sequence2[pos2 - 1 : pos2 + length]
            test_s2 = test_s2_raw if orientation == "forward" else reverse_complement(test_s2_raw)

            mm, _ = self._count_mismatches(test_s1, test_s2)
            if mm <= self.max_mismatches:
                pos1 -= 1
                pos2 -= 1
                length += 1
            else:
                break

        if length < self.min_length:
            return None

        seq1 = sequence1[pos1 : pos1 + length]
        seq2_raw = sequence2[pos2 : pos2 + length]
        seq2_cmp = seq2_raw if orientation == "forward" else reverse_complement(seq2_raw)
        mismatches, mismatch_positions = self._count_mismatches(seq1, seq2_cmp)

        return {
            "pos1": pos1,
            "pos2": pos2,
            "length": length,
            "seq1": seq1,
            "seq2": seq2_raw,
            "mismatches": mismatches,
            "mismatch_positions": mismatch_positions,
            "orientation": orientation,
        }

    # ------------------------------------------------------------------
    # Consolidation
    # ------------------------------------------------------------------

    def _consolidate_matches_between(self, matches: list[dict]) -> list[dict]:
        if not matches:
            return []

        matches_sorted = sorted(
            matches, key=lambda x: (-x["length"], x["mismatches"])
        )
        consolidated: list[dict] = []

        for match in matches_sorted:
            if not any(
                self._matches_overlap_between(match, ex) for ex in consolidated
            ):
                consolidated.append(match)

        return consolidated

    def _matches_overlap_between(self, m1: dict, m2: dict) -> bool:
        if m1["orientation"] != m2["orientation"]:
            return False
        s1_ovl = self._regions_overlap(
            m1["pos1"], m1["pos1"] + m1["length"],
            m2["pos1"], m2["pos1"] + m2["length"],
        )
        s2_ovl = self._regions_overlap(
            m1["pos2"], m1["pos2"] + m1["length"],
            m2["pos2"], m2["pos2"] + m2["length"],
        )
        return s1_ovl and s2_ovl

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------

    @staticmethod
    def _regions_overlap(s1: int, e1: int, s2: int, e2: int) -> bool:
        return not (e1 <= s2 or e2 <= s1)

    @staticmethod
    def _count_mismatches(seq1: str, seq2: str) -> tuple[int, list[int]]:
        mm = 0
        positions: list[int] = []
        for i, (a, b) in enumerate(zip(seq1, seq2)):
            if a != b:
                mm += 1
                positions.append(i)
        return mm, positions
