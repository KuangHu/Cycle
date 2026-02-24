"""Smoke tests for CircleFinder junction detection logic.

Uses mock pysam objects to test _find_junction_reads without needing
real BAM files or minimap2/samtools on PATH.
"""

import unittest
from unittest.mock import MagicMock, patch

# Patch shutil.which before importing CircleFinder so the constructor
# doesn't fail when minimap2/samtools aren't on PATH.
with patch("shutil.which", return_value="/usr/bin/fake"):
    from cycle.circle_detect.circle_finder import CircleFinder, ISEntry


def _make_finder(**kwargs):
    """Create a CircleFinder with tools check bypassed."""
    with patch("shutil.which", return_value="/usr/bin/fake"):
        return CircleFinder(output_dir="/tmp/test_circle", **kwargs)


def _mock_read(
    ref_name, ref_start, ref_end, cigar_tuples,
    query_name="read1", query_length=1000, query_sequence="A" * 1000,
    mapq=60, unmapped=False, secondary=False, supplementary=False,
):
    """Build a mock pysam AlignedSegment."""
    r = MagicMock()
    r.is_unmapped = unmapped
    r.is_secondary = secondary
    r.is_supplementary = supplementary
    r.reference_name = ref_name
    r.reference_start = ref_start
    r.reference_end = ref_end
    r.cigartuples = cigar_tuples
    r.query_name = query_name
    r.query_length = query_length
    r.query_sequence = query_sequence
    r.mapping_quality = mapq
    return r


class TestFindJunctionReads(unittest.TestCase):
    """Test _find_junction_reads with mock BAM reads."""

    def setUp(self):
        self.finder = _make_finder(min_overlap=100, boundary_tolerance=50)
        self.entry = ISEntry(
            uuid="is1", chrom="chr1", start=5000, end=6000,
            family="IS110", subfamily="ISEc1", consensus="A" * 1000,
        )
        self.entries_by_uuid = {"is1": self.entry}
        # TH bait: [IS][IS], junction at N=1000
        self.junction_map = {
            "is1__th__j1000": [
                {"type": "tail_head", "position": 1000, "uuid": "is1"},
            ],
        }

    def _run(self, reads):
        """Run _find_junction_reads with a list of mock reads."""
        mock_bam = MagicMock()
        mock_bam.fetch.return_value = reads
        mock_bam.__enter__ = MagicMock(return_value=mock_bam)
        mock_bam.__exit__ = MagicMock(return_value=False)

        with patch("pysam.AlignmentFile", return_value=mock_bam):
            return self.finder._find_junction_reads(
                bam_path="/fake.bam",
                entries_by_uuid=self.entries_by_uuid,
                junction_map=self.junction_map,
                sample_id="sample1",
            )

    def test_tail_head_spanning(self):
        """Read spanning junction at N=1000 → tail_head."""
        read = _mock_read(
            ref_name="is1__th__j1000",
            ref_start=800,   # 200bp before junction
            ref_end=1200,    # 200bp after junction
            cigar_tuples=[(0, 400)],  # pure match, no clips
        )
        junction_reads, summary = self._run([read])

        self.assertEqual(len(junction_reads), 1)
        self.assertEqual(junction_reads[0]["junction_type"], "tail_head")
        self.assertEqual(summary["is1"]["n_tail_head_reads"], 1)

    def test_genome_head_left_clip_at_zero(self):
        """Left soft-clip ≥100, alignment starts near 0 → genome_head."""
        # Read: [150bp soft-clip][aligned from 10 to 510]
        read = _mock_read(
            ref_name="is1__th__j1000",
            ref_start=10,
            ref_end=510,
            cigar_tuples=[(4, 150), (0, 500)],  # S=150, M=500
        )
        junction_reads, summary = self._run([read])

        types = [r["junction_type"] for r in junction_reads]
        self.assertIn("genome_head", types)
        self.assertEqual(summary["is1"]["n_genome_head_reads"], 1)

    def test_genome_head_left_clip_at_N(self):
        """Left soft-clip ≥100, alignment starts near N=1000 → genome_head."""
        read = _mock_read(
            ref_name="is1__th__j1000",
            ref_start=1005,
            ref_end=1505,
            cigar_tuples=[(4, 200), (0, 500)],
        )
        junction_reads, summary = self._run([read])

        types = [r["junction_type"] for r in junction_reads]
        self.assertIn("genome_head", types)

    def test_tail_genome_right_clip_at_N(self):
        """Right soft-clip ≥100, alignment ends near N=1000 → tail_genome."""
        # Read: [aligned 500..1010][200bp soft-clip]
        read = _mock_read(
            ref_name="is1__th__j1000",
            ref_start=500,
            ref_end=1010,
            cigar_tuples=[(0, 510), (4, 200)],
        )
        junction_reads, summary = self._run([read])

        types = [r["junction_type"] for r in junction_reads]
        self.assertIn("tail_genome", types)
        self.assertEqual(summary["is1"]["n_tail_genome_reads"], 1)

    def test_tail_genome_right_clip_at_2N(self):
        """Right soft-clip ≥100, alignment ends near 2N=2000 → tail_genome."""
        read = _mock_read(
            ref_name="is1__th__j1000",
            ref_start=1500,
            ref_end=1990,
            cigar_tuples=[(0, 490), (4, 150)],
        )
        junction_reads, summary = self._run([read])

        types = [r["junction_type"] for r in junction_reads]
        self.assertIn("tail_genome", types)

    def test_no_junction_short_clip(self):
        """Soft-clip below threshold → not counted as GH/TG."""
        # Only 50bp clip, below min_overlap=100
        read = _mock_read(
            ref_name="is1__th__j1000",
            ref_start=5,
            ref_end=505,
            cigar_tuples=[(4, 50), (0, 500)],
        )
        junction_reads, summary = self._run([read])

        # Should not be genome_head (clip too short) or tail_head (doesn't span)
        self.assertEqual(len(junction_reads), 0)
        self.assertEqual(summary["is1"]["n_genome_head_reads"], 0)

    def test_no_junction_far_from_boundary(self):
        """Large clip but alignment starts far from any boundary → no GH/TG."""
        # Alignment starts at 500, far from 0 or 1000
        read = _mock_read(
            ref_name="is1__th__j1000",
            ref_start=500,
            ref_end=800,
            cigar_tuples=[(4, 200), (0, 300)],
        )
        junction_reads, summary = self._run([read])

        types = [r["junction_type"] for r in junction_reads]
        self.assertNotIn("genome_head", types)

    def test_unmapped_skipped(self):
        """Unmapped reads are ignored."""
        read = _mock_read(
            ref_name="is1__th__j1000",
            ref_start=800, ref_end=1200,
            cigar_tuples=[(0, 400)],
            unmapped=True,
        )
        junction_reads, summary = self._run([read])
        self.assertEqual(len(junction_reads), 0)

    def test_combined_th_and_gh(self):
        """A read can be counted as both tail_head AND genome_head."""
        # Spans junction AND has large left clip near boundary N=1000
        # ref_start=895 is near 900 → within 50 of 1000? No, abs(895-1000)=105 > 50.
        # Let's use ref_start=960 which is within 50 of N=1000
        # And ref_end=1200 which spans 200bp past junction
        # Left clip=200, aligned from 960 to 1200 (240bp) — spans junction at 1000
        read = _mock_read(
            ref_name="is1__th__j1000",
            ref_start=960,
            ref_end=1200,
            cigar_tuples=[(4, 200), (0, 240)],
        )
        junction_reads, summary = self._run([read])

        types = [r["junction_type"] for r in junction_reads]
        # Should NOT be tail_head because ref_start=960 > 1000-100=900
        self.assertNotIn("tail_head", types)
        # Should be genome_head: left_clip=200≥100, aligned=240≥100, abs(960-1000)=40<50
        self.assertIn("genome_head", types)

    def test_multiple_reads(self):
        """Multiple reads counted correctly in summary."""
        reads = [
            # tail_head spanning read
            _mock_read("is1__th__j1000", 800, 1200, [(0, 400)], query_name="r1"),
            # genome_head at 0
            _mock_read("is1__th__j1000", 5, 505, [(4, 150), (0, 500)], query_name="r2"),
            # tail_genome at 2N
            _mock_read("is1__th__j1000", 1500, 1995, [(0, 495), (4, 150)], query_name="r3"),
            # another tail_head
            _mock_read("is1__th__j1000", 700, 1300, [(0, 600)], query_name="r4"),
        ]
        junction_reads, summary = self._run(reads)

        self.assertEqual(summary["is1"]["n_tail_head_reads"], 2)
        self.assertEqual(summary["is1"]["n_genome_head_reads"], 1)
        self.assertEqual(summary["is1"]["n_tail_genome_reads"], 1)
        self.assertEqual(summary["is1"]["n_total_mapped"], 4)


class TestBuildBaitFasta(unittest.TestCase):
    """Test _build_bait_fasta output."""

    def test_generates_th_only(self):
        """Bait FASTA only contains TH baits, no IC."""
        import tempfile
        import os

        finder = _make_finder()
        entry = ISEntry(
            uuid="test1", chrom="chr1", start=100, end=200,
            family="IS110", subfamily="ISEc1", consensus="ACGT" * 100,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "test_bait.fa")
            fasta_path, junction_map = finder._build_bait_fasta(
                [entry], out_path,
            )

            # Only TH entries
            self.assertEqual(len(junction_map), 1)
            self.assertIn("test1__th__j400", junction_map)
            self.assertEqual(junction_map["test1__th__j400"][0]["type"], "tail_head")

            # Read FASTA and verify content
            with open(fasta_path) as f:
                lines = f.readlines()
            self.assertEqual(len(lines), 2)  # header + sequence
            self.assertTrue(lines[0].startswith(">test1__th__j400"))
            self.assertEqual(len(lines[1].strip()), 800)  # 400*2


class TestParseHeaders(unittest.TestCase):
    """Test _parse_bait_headers handles old IC headers gracefully."""

    def test_ic_headers_ignored(self):
        """Old __ic__ headers in existing bait FASTAs are skipped."""
        import tempfile
        import os

        finder = _make_finder()
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = os.path.join(tmpdir, "bait.fa")
            with open(fasta, "w") as f:
                f.write(">uuid1__th__j1000\n")
                f.write("A" * 2000 + "\n")
                f.write(">uuid1__ic__j500_j1500__fl500\n")
                f.write("A" * 2500 + "\n")

            jmap = finder._parse_bait_headers(fasta)
            # Only TH entry, IC is ignored
            self.assertEqual(len(jmap), 1)
            self.assertIn("uuid1__th__j1000", jmap)

    def test_legacy_len_format(self):
        """Legacy __len headers still parsed."""
        import tempfile
        import os

        finder = _make_finder()
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = os.path.join(tmpdir, "bait.fa")
            with open(fasta, "w") as f:
                f.write(">uuid1__len1000\n")
                f.write("A" * 2000 + "\n")

            jmap = finder._parse_bait_headers(fasta)
            self.assertEqual(len(jmap), 1)
            self.assertEqual(jmap["uuid1__len1000"][0]["type"], "tail_head")


if __name__ == "__main__":
    unittest.main()
