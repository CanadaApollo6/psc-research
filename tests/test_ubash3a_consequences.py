"""Boundary and biological-coordinate tests for conditional splice reconstruction."""

import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("ubash3a_consequences", ROOT / "scripts/analyze_ubash3a_consequences.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def exon(identifier, start, end, strand=1):
    return {"id": identifier, "start": start, "end": end, "strand": strand}


class ReferenceReconstructionTests(unittest.TestCase):
    def test_exon_concatenation_excludes_intronic_bases_and_preserves_map(self):
        cdna, positions, blocks = m.reconstruct_reference("ATGAAA" + "G" * 10 + "CCCTAA", 100, [exon("a", 100, 105), exon("b", 116, 121)])
        self.assertEqual(cdna, "ATGAAACCCTAA")
        self.assertEqual(positions, list(range(100, 106)) + list(range(116, 122)))
        self.assertEqual([b["cdna_end_0based_exclusive"] for b in blocks], [6, 12])

    def test_exon_outside_reference_fails(self):
        with self.assertRaises(ValueError):
            m.reconstruct_reference("ATG", 100, [exon("a", 99, 102)])

    def test_overlap_and_wrong_strand_fail(self):
        with self.assertRaises(ValueError):
            m.ordered_exons({"strand": 1, "Exon": [exon("a", 100, 105), exon("b", 105, 110)]})
        with self.assertRaises(ValueError):
            m.ordered_exons({"strand": -1, "Exon": [exon("a", 100, 105, -1)]})

    def test_cds_convention_includes_stop(self):
        cdna = "CCCATGGCTTAAGGG"
        positions = list(range(100, 100 + len(cdna)))
        start, translated, convention = m.reference_coding_bounds(cdna, positions, {"start": 103, "end": 111, "length": 2}, "ATGGCTTAA", "MA")
        self.assertEqual(start, 3)
        self.assertEqual(translated["stop_end_0based_exclusive"], 12)
        self.assertEqual(convention, "lookup_translation_bounds_equal_API_CDS")

    def test_cds_convention_excludes_stop_in_lookup(self):
        cdna = "CCCATGGCTTAAGGG"
        _, _, convention = m.reference_coding_bounds(cdna, list(range(100, 115)), {"start": 103, "end": 108, "length": 2}, "ATGGCTTAA", "MA")
        self.assertEqual(convention, "API_CDS_adds_stop_after_lookup_translation_bounds")

    def test_disagreeing_cds_or_protein_is_not_silently_trimmed(self):
        for cds, protein in [("ATGGCGTAA", "MA"), ("ATGGCTTAA", "ML")]:
            with self.assertRaises(ValueError):
                m.reference_coding_bounds("ATGGCTTAA", list(range(100, 109)), {"start": 100, "end": 108, "length": 2}, cds, protein)


class TranslationTests(unittest.TestCase):
    def test_all_three_stop_codons_and_first_stop_only(self):
        for stop in ["TAA", "TAG", "TGA"]:
            value = m.translate_from_start("ATGTGG" + stop + "GCC", 0)
            self.assertEqual(value["peptide"], "MW")
            self.assertEqual(value["stop_start_0based"], 6)
            self.assertEqual(value["stop_end_0based_exclusive"], 9)

    def test_incomplete_final_codon_is_not_a_stop(self):
        value = m.translate_from_start("ATGGCTA", 0)
        self.assertEqual(value["peptide"], "MA")
        self.assertTrue(value["no_stop_in_transcript"])
        self.assertIsNone(value["stop_codon"])

    def test_start_and_ambiguous_codons_fail(self):
        for cdna in ["CTGTAA", "ATGNNNTAA"]:
            with self.assertRaises(ValueError):
                m.translate_from_start(cdna, 0)

    def test_coding_retention_changes_frame_after_junction(self):
        value = m.translate_from_start("ATGAAA" + "GTAGT" + "CCCTAA", 0)
        self.assertEqual(value["peptide"], "MKVVP")
        self.assertTrue(value["no_stop_in_transcript"])


class RetentionTests(unittest.TestCase):
    def setUp(self):
        self.cdna, self.positions, self.blocks = m.reconstruct_reference("ATGAAA" + "G" * 40 + "CCCTAA", 100, [exon("a", 100, 105), exon("b", 146, 151)])

    def test_retention_moves_one_boundary_and_downstream_coordinates(self):
        segment = "GTAGT" + "A" * 24
        altered, positions, blocks = m.insert_retained_segment(self.cdna, self.positions, self.blocks, 6, segment, 106)
        self.assertEqual(altered, "ATGAAA" + segment + "CCCTAA")
        self.assertEqual(len(altered) - len(self.cdna), 29)
        self.assertEqual(positions[6:35], list(range(106, 135)))
        self.assertEqual(blocks[0]["genomic_end_1based"], 134)
        self.assertEqual(blocks[1]["genomic_start_1based"], 146)
        self.assertEqual(blocks[1]["cdna_start_0based"], 35)

    def test_nonjunction_insertion_and_wrong_genomic_start_fail(self):
        with self.assertRaises(ValueError):
            m.insert_retained_segment(self.cdna, self.positions, self.blocks, 5, "GTA", 106)
        with self.assertRaises(ValueError):
            m.insert_retained_segment(self.cdna, self.positions, self.blocks, 6, "GTA", 107)

    def test_retention_cannot_consume_entire_intron(self):
        with self.assertRaises(ValueError):
            m.insert_retained_segment(self.cdna, self.positions, self.blocks, 6, "A" * 40, 106)

    def test_prefix_comparison_preserves_shortened_identical_prefix(self):
        self.assertEqual(m.common_prefix_length("MAGA", "MAGATT"), 4)
        self.assertEqual(m.common_prefix_length("MARG", "MAGG"), 2)


class NmdFeatureTests(unittest.TestCase):
    settings = {"long_stop_exon_gt_nt": 400, "start_proximity_coding_nt": 150, "primary_stop_end_to_final_junction_gt_nt": 55, "sensitivity_stop_end_to_final_junction_gt_nt": 50}
    blocks = [{"exon_number": 1, "cdna_start_0based": 0, "cdna_end_0based_exclusive": 200}, {"exon_number": 2, "cdna_start_0based": 200, "cdna_end_0based_exclusive": 400}]

    def feature(self, end, blocks=None):
        return m.nmd_features({"stop_start_0based": end - 3, "stop_end_0based_exclusive": end}, 0, blocks or self.blocks, self.settings)

    def test_55_boundary_and_50_sensitivity_are_explicit(self):
        value = self.feature(145)
        self.assertEqual(value["stop_end_to_final_junction_nt"], 55)
        self.assertFalse(value["nmd_55nt_flag"])
        self.assertTrue(value["nmd_50nt_flag"])
        self.assertTrue(self.feature(144)["nmd_55nt_flag"])
        self.assertFalse(self.feature(150)["nmd_50nt_flag"])

    def test_last_exon_is_not_flagged(self):
        value = self.feature(250)
        self.assertTrue(value["stop_in_final_exon"])
        self.assertFalse(value["nmd_55nt_flag"])
        self.assertIsNone(value["stop_end_to_next_junction_nt"])

    def test_long_exon_and_start_proximity_are_separate_flags(self):
        blocks = [{"exon_number": 1, "cdna_start_0based": 0, "cdna_end_0based_exclusive": 401}, {"exon_number": 2, "cdna_start_0based": 401, "cdna_end_0based_exclusive": 800}]
        value = self.feature(100, blocks)
        self.assertTrue(value["stop_exon_gt_400_nt"])
        self.assertTrue(value["stop_start_within_first_150_coding_nt"])
        self.assertTrue(value["nmd_55nt_flag"])

    def test_no_stop_does_not_create_an_nmd_probability(self):
        value = m.nmd_features({"stop_start_0based": None, "stop_end_0based_exclusive": None}, 0, self.blocks, self.settings)
        self.assertEqual(value["nmd_status"], "unavailable_no_stop")
        self.assertNotIn("nmd_55nt_flag", value)


class ProteinFeatureTests(unittest.TestCase):
    def test_genomic_codon_projection_handles_isoform_offsets(self):
        anchor = {(100, 101, 102): 1, (200, 201, 202): 2, (300, 301, 302): 3}
        rows = m.codon_alignment("MT", [100, 101, 102, 300, 301, 302], 0, anchor, "MAT")
        self.assertEqual([r["uniprot_position_1based"] for r in rows], [1, 3])
        feature = {"type": "Domain", "description": "example", "location": {"start": {"value": 1, "modifier": "EXACT"}, "end": {"value": 3, "modifier": "EXACT"}}}
        summary = m.summarize_feature(feature, rows)
        self.assertEqual(summary["reference_residues_present"], 2)
        self.assertEqual(summary["sequence_status"], "partial_reference_sequence")

    def test_same_position_without_same_genomic_codon_is_unmapped(self):
        rows = m.codon_alignment("M", [101, 102, 103], 0, {(100, 101, 102): 1}, "M")
        self.assertIsNone(rows[0]["uniprot_position_1based"])

    def test_mapped_codon_amino_acid_disagreement_fails(self):
        with self.assertRaises(ValueError):
            m.codon_alignment("L", [100, 101, 102], 0, {(100, 101, 102): 1}, "M")


if __name__ == "__main__":
    unittest.main()
