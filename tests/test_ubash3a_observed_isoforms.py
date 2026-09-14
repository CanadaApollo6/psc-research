"""Synthetic truths for deposited RNA structure identification and selection."""

import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("observed_isoforms", ROOT / "scripts/analyze_ubash3a_observed_isoforms.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class AnnotationAttributesTests(unittest.TestCase):
    def test_gtf_quoted_attributes(self):
        self.assertEqual(m.attributes('gene_id "ENSG1.2"; transcript_id "abc";'), {"gene_id": "ENSG1.2", "transcript_id": "abc"})

    def test_equals_inside_quoted_gtf_value(self):
        self.assertEqual(m.attributes('transcript_id "a=b";')["transcript_id"], "a=b")

    def test_gff3_percent_encoding_and_multiple_parents(self):
        attrs = m.attributes("ID=exon1;Parent=tx1,tx2;Name=RNA%20one")
        self.assertEqual(attrs["Name"], "RNA one")
        self.assertEqual(m.transcript_ids("exon", attrs), ["tx1", "tx2"])

    def test_gff3_transcript_and_gene_are_distinguished(self):
        self.assertEqual(m.transcript_ids("mRNA", {"ID": "tx1", "Parent": "g1"}), ["tx1"])
        self.assertEqual(m.transcript_ids("gene", {"ID": "g1"}), [])

    def test_conflicting_duplicate_ids_fail(self):
        with self.assertRaises(ValueError):
            m.attributes('transcript_id "one"; transcript_id "two";')

    def test_malformed_attributes_fail(self):
        with self.assertRaises(ValueError):
            m.attributes('transcript_id "unfinished')


class JunctionIdentityTests(unittest.TestCase):
    def setUp(self):
        self.plan = {"normal_upstream_exon_end_1based": 120, "retained_upstream_exon_end_1based": 149, "downstream_exon_start_1based": 200}

    def test_exact_plus29_event(self):
        event, index, suffix = m.event_and_suffix(m.junction_chain([(100, 149), (200, 220), (300, 340)]), self.plan)
        self.assertEqual((event, index, suffix), ("plus29", 0, ((220, 300),)))

    def test_normal_splicing_is_separate(self):
        event, _, _ = m.event_and_suffix(m.junction_chain([(100, 120), (200, 220)]), self.plan)
        self.assertEqual(event, "normal")

    def test_nearby_donor_and_shifted_acceptor_do_not_count(self):
        for exons in [[(100, 148), (200, 220)], [(100, 149), (201, 220)]]:
            self.assertEqual(m.event_and_suffix(m.junction_chain(exons), self.plan)[0], "other")

    def test_full_intron_retention_is_not_plus29(self):
        self.assertEqual(m.event_and_suffix(m.junction_chain([(100, 220), (300, 340)]), self.plan)[0], "other")

    def test_exon_skip_changes_downstream_chain(self):
        a = m.event_and_suffix(m.junction_chain([(100, 149), (200, 220), (300, 320), (400, 440)]), self.plan)[2]
        b = m.event_and_suffix(m.junction_chain([(100, 149), (200, 220), (400, 440)]), self.plan)[2]
        self.assertNotEqual(a, b)

    def test_terminal_coordinate_variation_does_not_change_internal_chain(self):
        self.assertEqual(m.junction_chain([(95, 120), (200, 230)]), m.junction_chain([(100, 120), (200, 220)]))

    def test_overlap_fails(self):
        with self.assertRaises(ValueError):
            m.junction_chain([(100, 200), (200, 300)])

    def test_exonic_coordinate_boundaries_and_gap(self):
        exons = [(100, 120), (200, 220)]
        for position in [100, 120, 200, 220]:
            self.assertTrue(m.interval_contains(exons, position))
        for position in [99, 121, 150, 199, 221]:
            self.assertFalse(m.interval_contains(exons, position))


class UnquotedAttributeTests(unittest.TestCase):
    def test_unquoted_support_value_remains_original_string(self):
        self.assertEqual(m.attributes('GENCODE_transcript_support_level 1; score -2.5e-3;'),
                         {'GENCODE_transcript_support_level': '1', 'score': '-2.5e-3'})

    def test_unquoted_multiple_tokens_are_not_silently_accepted(self):
        with self.assertRaises(ValueError):
            m.attributes('malformed two values;')


class LocusSelectionTests(unittest.TestCase):
    def collect(self, text):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "source.gtf"
            path.write_text(text)
            return m.collect_locus(path, 100, 300)

    def test_coordinate_lookup_recovers_anonymous_ids_and_distant_terminal_exon(self):
        text = '21\tsource\texon\t110\t120\t.\t+\t.\tgene_id "PB.1"; transcript_id "PB.1.1";\n21\tsource\texon\t900\t950\t.\t+\t.\tgene_id "PB.1"; transcript_id "PB.1.1";\n'
        models, rows = self.collect(text)
        self.assertEqual(models["PB.1.1"]["exons"], [(110, 120), (900, 950)])
        self.assertEqual(rows, 2)
        self.assertTrue(models["PB.1.1"]["coordinate_selected"])
        self.assertFalse(models["PB.1.1"]["gene_label_selected"])

    def test_gene_label_on_discordant_chromosome_is_retained(self):
        models, _ = self.collect('chr22\tsource\texon\t110\t120\t.\t-\t.\tgene_id "ENSG00000160185.1"; transcript_id "wrong";\n')
        self.assertEqual(models["wrong"]["chromosomes"], {"22"})
        self.assertTrue(models["wrong"]["gene_label_selected"])

    def test_gene_symbol_substring_is_not_an_exact_match(self):
        models, _ = self.collect('chr22\tsource\texon\t110\t120\t.\t+\t.\tgene_name "UBASH3AP1"; transcript_id "other";\n')
        self.assertEqual(models, {})

    def test_source_without_exons_fails_instead_of_reporting_absence(self):
        with self.assertRaises(ValueError):
            self.collect('21\tsource\ttranscript\t110\t120\t.\t+\t.\ttranscript_id "empty";\n')


if __name__ == "__main__":
    unittest.main()
