from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import analyze_ubash3a_proteomics_spectra as analysis
import verify_ubash3a_proteomics_search as verification


class TargetDecoyTests(unittest.TestCase):
    def row(self, scan, score, decoy=False):
        return {"file": "one", "scan": scan, "score": score, "is_decoy": decoy}

    def test_zero_decoys_does_not_mean_zero_error(self):
        result = analysis.q_values([self.row(i, i) for i in range(1, 11)])
        self.assertTrue(all(row["q_value"] == 0.1 for row in result))

    def test_equal_score_group_does_not_depend_on_input_order(self):
        rows = [self.row(1, 1), self.row(2, 1, True), self.row(3, 2), self.row(4, 3)]
        first, second = analysis.q_values(rows), analysis.q_values(list(reversed(rows)))
        self.assertEqual(first, second)
        self.assertEqual(first[0]["q_value"], first[1]["q_value"])
        self.assertAlmostEqual(first[0]["q_value"], 2 / 3)

    def test_no_targets_returns_one(self):
        self.assertEqual(analysis.q_values([self.row(1, 1, True)])[0]["q_value"], 1)

    def test_q_values_are_monotone(self):
        result = analysis.q_values([self.row(1, 1), self.row(2, 2), self.row(3, 3, True), self.row(4, 4)])
        self.assertEqual([r["q_value"] for r in result], [0.5, 0.5, 2 / 3, 2 / 3])

    def test_any_decoy_reference_is_conservative(self):
        self.assertTrue(analysis.has_decoy("P_123,DECOY_P_456"))
        self.assertFalse(analysis.has_decoy("P_123,P_456"))

    def test_multiple_charges_are_one_spectrum_and_ties_favor_decoy(self):
        rows = [
            {"file": "one", "scan": "12", "charge": "2", "num": "1", "e-value": "0.01", "xcorr": "3", "protein": "P_1", "plain_peptide": "PEPTIDER"},
            {"file": "one", "scan": "12", "charge": "3", "num": "1", "e-value": "0.01", "xcorr": "3", "protein": "DECOY_P_1", "plain_peptide": "EDITPEPR"},
        ]
        result = analysis.best_per_scan(rows)
        self.assertEqual(len(result), 1)
        self.assertTrue(result[0]["is_decoy"])
        self.assertEqual(result[0]["charge_queries_for_scan"], 2)

    def test_lower_ranks_do_not_inflate_target_count(self):
        rows = [
            {"file": "one", "scan": "12", "charge": "2", "num": "1", "e-value": "0.01", "xcorr": "3", "protein": "P_1", "plain_peptide": "PEPTIDER"},
            {"file": "one", "scan": "12", "charge": "2", "num": "2", "e-value": "0.02", "xcorr": "2", "protein": "P_2", "plain_peptide": "DIFFERENT"},
        ]
        self.assertEqual(len(analysis.best_per_scan(rows)), 1)

    def test_repeated_rank_one_is_one_conservative_decoy_win(self):
        base = {"file": "one", "scan": "1127", "charge": "2", "num": "1", "e-value": "0.01", "xcorr": "3"}
        rows = [{**base, "protein": "P_1", "plain_peptide": "PEPTIDER", "query_hit_index": 1},
                {**base, "protein": "DECOY_P_2", "plain_peptide": "EDITPEPR", "query_hit_index": 2}]
        result = analysis.best_per_scan(rows)
        self.assertEqual(len(result), 1)
        self.assertTrue(result[0]["is_decoy"])
        self.assertEqual(result[0]["protein"], "DECOY_P_2")
        self.assertEqual(result[0]["rank_one_assignment_count"], 2)
        self.assertEqual(result[0]["reported_top_score_tie_count"], 2)
        self.assertEqual(result, analysis.best_per_scan(list(reversed(rows))))

    def test_repeated_target_rank_does_not_add_a_second_target(self):
        base = {"file": "one", "scan": "1127", "charge": "2", "num": "1", "e-value": "0.01", "xcorr": "3", "protein": "P_1"}
        rows = [{**base, "plain_peptide": "PEPTIDER", "query_hit_index": 1},
                {**base, "plain_peptide": "ANOTHER", "query_hit_index": 2}]
        result = analysis.q_values(analysis.best_per_scan(rows))
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["plain_peptide"], "ANOTHER")
        self.assertEqual(result[0]["q_value"], 1)


class SpectrumIdentityTests(unittest.TestCase):
    def with_xml(self, body, check):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "output.pep.xml"
            path.write_text('<msms_pipeline_analysis xmlns="http://regis-web.systemsbiology.net/pepXML">' + body + '</msms_pipeline_analysis>')
            check(path)

    def test_native_id_is_used_instead_of_exported_scan(self):
        body = '<spectrum_query start_scan="67" end_scan="67" assumed_charge="2" spectrumNativeID="controllerType=0 controllerNumber=1 scan=1127"/>'
        self.with_xml(body, lambda path: self.assertEqual(analysis.query_scan_map(path), {("67", "2"): "1127"}))

    def test_missing_native_id_is_not_silently_replaced(self):
        body = '<spectrum_query start_scan="67" end_scan="67" assumed_charge="2"/>'
        def check(path):
            with self.assertRaises(KeyError):
                analysis.query_scan_map(path)
        self.with_xml(body, check)

    def test_conflicting_exported_ids_for_one_original_scan_fail(self):
        body = ''.join(f'<spectrum_query start_scan="{scan}" end_scan="{scan}" assumed_charge="2" spectrumNativeID="scan=1127"/>' for scan in (67, 68))
        def check(path):
            with self.assertRaises(AssertionError):
                analysis.query_scan_map(path)
        self.with_xml(body, check)

    def test_xml_parser_preserves_dense_rank_and_export_position(self):
        hits = ''.join(f'<search_hit hit_rank="1" peptide="{peptide}" protein="P_{index}" calc_neutral_pep_mass="900.0" num_matched_ions="4" tot_num_ions="12"><search_score name="expect" value="0.01"/><search_score name="xcorr" value="3.0"/></search_hit>' for index, peptide in enumerate(("PEPTIDER", "ANOTHER"), 1))
        body = '<spectrum_query start_scan="67" end_scan="67" assumed_charge="2" precursor_neutral_mass="900.0" spectrumNativeID="scan=1127"><search_result>' + hits + '</search_result></spectrum_query>'
        def check(path):
            rows, _ = verification.parse_pepxml(path, "one")
            self.assertEqual([row["rank"] for row in rows], [1, 1])
            self.assertEqual([row["query_hit_index"] for row in rows], [1, 2])
            self.assertEqual([row["comet_scan"] for row in rows], [67, 67])
            self.assertEqual([row["scan"] for row in rows], [1127, 1127])
        self.with_xml(body, check)

    def test_pin_preserves_decoy_alternative_when_parent_label_is_target(self):
        text = "SpecId\tLabel\tScanNr\tExpMass\tCalcMass\tXcorr\tPeptide\tProteins\nfile_107_2_3\t1\t107\t792.461324\t792.457393\t0.649\tR.KQKSSSK.A\tP_target\tDECOY_P_other\n"
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "output.pin"
            path.write_text(text)
            result = verification.parse_pin(path)[(107, 2, 3)]
            self.assertEqual(result["PIN_reported_label"], 1)
            self.assertTrue(result["hit_decoy"])
            self.assertEqual(result["protein_refs"], ("DECOY_P_other", "P_target"))

    def test_pin_scan_number_must_match_explicit_spectrum_key(self):
        text = "SpecId\tLabel\tScanNr\tExpMass\tCalcMass\tXcorr\tPeptide\tProteins\nfile_107_2_3\t1\t108\t792.461324\t792.457393\t0.649\tR.KQKSSSK.A\tP_target\n"
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "output.pin"
            path.write_text(text)
            with self.assertRaises(AssertionError):
                verification.parse_pin(path)


if __name__ == "__main__":
    unittest.main()
