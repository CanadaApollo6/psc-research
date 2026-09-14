"""Synthetic truth checks for raw-read eligibility, orientation, and FASTQ integrity."""

import io
from pathlib import Path
import random
import sys
import unittest

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import analyze_ubash3a_raw_reads as audit
from run_ubash3a_raw_read_alignment import fastq_records


def read_at_junction(event="normal", left=30, right=30, flag=0, ts="+", mapq=60):
    donor, acceptor = audit.NORMAL if event == "normal" else audit.PLUS29
    read = pysam.AlignedSegment()
    read.query_name = "synthetic"
    read.reference_start = donor-left
    read.flag = flag
    read.mapping_quality = mapq
    read.cigartuples = [(0, left), (3, acceptor-donor-1), (0, right)]
    read.query_sequence = "A"*(left+right)
    read.set_tag("ts", ts)
    return read


class RawReadTruthTests(unittest.TestCase):
    def test_both_query_orientations_infer_positive_genomic_transcript(self):
        for read in [read_at_junction(), read_at_junction(flag=16, ts="-")]:
            self.assertEqual(audit.genomic_transcript_strand(read), "+")
            self.assertTrue(audit.classify(read, False)["qualified"])

    def test_literal_ts_plus_on_reverse_read_is_negative_transcript(self):
        read = read_at_junction(flag=16, ts="+")
        self.assertEqual(audit.genomic_transcript_strand(read), "-")
        self.assertFalse(audit.classify(read, False)["qualified"])

    def test_exact_plus29_event_and_anchor_boundary(self):
        read = read_at_junction("plus29", left=20, right=20)
        info = audit.classify(read, False)
        self.assertEqual(info["event"], "plus29")
        self.assertTrue(info["qualified"])
        read = read_at_junction("plus29", left=19)
        self.assertFalse(audit.classify(read, False)["qualified"])
        self.assertTrue(audit.classify(read, False, 10)["qualified"])

    def test_insertion_next_to_junction_interrupts_anchor(self):
        read = read_at_junction()
        gap = audit.NORMAL[1]-audit.NORMAL[0]-1
        read.cigartuples = [(0, 20), (1, 1), (0, 10), (3, gap), (0, 30)]
        read.query_sequence = "A"*61
        self.assertFalse(audit.classify(read, False)["qualified"])
        self.assertTrue(audit.classify(read, False, 10)["qualified"])

    def test_deletion_cannot_be_counted_as_a_splice(self):
        read = read_at_junction()
        read.cigartuples = [(0, 30), (2, audit.NORMAL[1]-audit.NORMAL[0]-1), (0, 30)]
        self.assertEqual(audit.classify(read, False)["event"], "other")

    def test_other_donor_cannot_be_rounded_to_target(self):
        read = read_at_junction()
        read.reference_start += 1
        self.assertEqual(audit.classify(read, False)["event"], "other")

    def test_chimeric_secondary_and_low_mapq_are_retained_but_unqualified(self):
        for kwargs in [{"flag": 256}, {"flag": 2048}, {"mapq": 19}]:
            self.assertFalse(audit.classify(read_at_junction(**kwargs), False)["qualified"])
        self.assertFalse(audit.classify(read_at_junction(), True)["qualified"])
        read = read_at_junction()
        read.set_tag("SA", "22,100,+,60M,50,0;")
        self.assertFalse(audit.classify(read, False)["qualified"])

    def test_fragment_does_not_match_complete_downstream_chain(self):
        read = read_at_junction("plus29")
        ref = {"event": "plus29", "target_label": "T29", "chain": [audit.PLUS29, (42437580, 42437855)],
               "exons": [(42434832, 42434983), (42437488, 42437580), (42437855, 42440816)]}
        self.assertEqual(audit.reference_matches(audit.classify(read, False), [ref], 20), ([], [], []))

    def test_dynamic_programming_favors_the_correct_distinct_junction_template(self):
        rng = random.Random(814)
        class Genome:
            sequence = "".join(rng.choice("ACGT") for _ in range(3000))
            def fetch(self, chromosome, start, end):
                self_outer.assertEqual(chromosome, "21")
                return self.sequence[start-42434800:end-42434800]
        self_outer = self
        genome = Genome()
        for event, sign in [("normal", -1), ("plus29", 1)]:
            donor, acceptor = audit.NORMAL if event == "normal" else audit.PLUS29
            read = read_at_junction(event, left=donor-42434904+1, right=50)
            read.query_sequence = genome.fetch("21", 42434903, donor) + genome.fetch("21", acceptor-1, acceptor+49)
            scores = audit.local_template_scores(read, genome)
            self.assertGreater(sign*scores["dp_plus29_minus_normal"], 0)

    def test_terminal_run_descriptors(self):
        self.assertEqual(audit.end_run("TTTACAAA", "T", False), 3)
        self.assertEqual(audit.end_run("TTTACAAA", "A"), 3)


class FastqIntegrityTests(unittest.TestCase):
    def test_valid_record_preserves_original_bytes(self):
        source = b"@test comment\nACGTN\n+\nIIIII\n"
        self.assertEqual(list(fastq_records(io.BytesIO(source))), [(source, 5)])

    def test_truncated_and_mismatched_records_fail(self):
        for source in [b"@x\nAC\n+\n", b"@x\nAC\n+\nI\n", b"x\nAC\n+\nII\n"]:
            with self.assertRaises(ValueError):
                list(fastq_records(io.BytesIO(source)))

    def test_wrapped_or_invalid_base_records_fail_explicitly(self):
        for source in [b"@x\nAC\nGT\n+\nIIII\n", b"@x\nAZ\n+\nII\n"]:
            with self.assertRaises(ValueError):
                list(fastq_records(io.BytesIO(source)))


if __name__ == "__main__":
    unittest.main()
