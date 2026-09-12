"""Unit checks for the ETS2 summary benchmark's conservative joins and signs."""

import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("benchmark_ets2_experiments", ROOT / "scripts/benchmark_ets2_experiments.py")
BENCHMARK = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = BENCHMARK
assert SPEC.loader is not None
SPEC.loader.exec_module(BENCHMARK)


class ETS2ExperimentBenchmarkTests(unittest.TestCase):
    def test_symbol_xref_excludes_both_sides_of_ambiguous_identity(self):
        xref = BENCHMARK.build_one_to_one_xref([
            ("ENSG00000000001", "ONE"),
            ("ENSG00000000002", "AMBIGUOUS"),
            ("ENSG00000000003", "AMBIGUOUS"),
            ("ENSG00000000004", "ALIAS_A"),
            ("ENSG00000000004", "ALIAS_B"),
            ("ENSG00000000005", "USABLE"),
        ])
        self.assertEqual(xref.mapping, {"ALIAS_A": "ENSG00000000004", "ALIAS_B": "ENSG00000000004", "ONE": "ENSG00000000001", "USABLE": "ENSG00000000005"})
        self.assertEqual(xref.ambiguous["AMBIGUOUS"], ("ENSG00000000002", "ENSG00000000003"))
        self.assertEqual(xref.reverse_alias_targets["ENSG00000000004"], ("ALIAS_A", "ALIAS_B"))

    def test_rank_mapping_does_not_double_count_aliases_and_retains_accounting(self):
        xref = BENCHMARK.build_one_to_one_xref([
            ("ENSG00000000001", "A"),
            ("ENSG00000000001", "ALIAS_A"),
            ("ENSG00000000002", "B"),
        ])
        rows = [
            {"source_identifier": "A", "raw_value": "1.0", "value": 1.0, "source_row": 1},
            {"source_identifier": "ALIAS_A", "raw_value": "2.0", "value": 2.0, "source_row": 2},
            {"source_identifier": "B", "raw_value": "0", "value": 0.0, "source_row": 3},
            {"source_identifier": "MISSING", "raw_value": "NA", "value": None, "source_row": 4},
        ]
        mapped, audited, accounting = BENCHMARK.map_rank_rows(rows, xref)
        self.assertEqual(mapped, {
            "ENSG00000000001": rows[0] | {"ensembl_id": "ENSG00000000001", "mapping_status": "mapped_one_to_one"},
            "ENSG00000000002": rows[2] | {"ensembl_id": "ENSG00000000002", "mapping_status": "mapped_one_to_one"},
        })
        self.assertEqual([row["mapping_status"] for row in audited], ["mapped_one_to_one", "duplicate_target_ensembl", "mapped_one_to_one", "missing_xref"])
        self.assertEqual(accounting["input_rows"], 4)
        self.assertEqual(accounting["mapped_unique_ensembl_count"], 2)
        self.assertEqual(accounting["zero_value_rows"], 1)
        self.assertEqual(accounting["na_value_rows"], 1)

    def test_unresolved_meki_orientation_preserves_overlap_but_withholds_direction(self):
        experiments = {
            "g1": {"by_id": {
                "ENSG1": {"t": -2.0, "logFC": -1.0},
                "ENSG2": {"t": 1.0, "logFC": 0.5},
            }},
            "MEKi0.1": {"by_id": {
                "ENSG1": {"value": 3.0},
                "ENSG2": {"value": -4.0},
            }},
        }
        row = BENCHMARK.build_comparison_row("g1", "MEKi0.1", experiments, BENCHMARK.default_mek_orientation())
        self.assertEqual(row["n_shared_unique"], 2)
        self.assertEqual(row["n_t"], 2)
        self.assertIsNone(row["spearman_signed_t"])
        self.assertIsNone(row["direction_agreement_t"])
        self.assertEqual(row["directional_status"], "unavailable_unresolved_orientation")

    def test_unresolved_rank_summary_keeps_original_numeric_median(self):
        source_set = [
            {"ensembl_id": "ENSG1", "is_excluded": False},
            {"ensembl_id": "ENSG2", "is_excluded": False},
        ]
        experiments = {
            experiment: {"by_id": {}}
            for experiment in ("g1", "g2", "chr21", "OE250", "OE500", "MEKi0.1", "MEKi0.5")
        }
        experiments["MEKi0.1"] = {"by_id": {"ENSG1": {"value": -3.0}, "ENSG2": {"value": 5.0}}}
        summaries = BENCHMARK.build_program_summary(
            source_set,
            experiments,
            BENCHMARK.default_mek_orientation(),
            {"source_set_n": 2, "excluded_n": 0},
        )
        summary = next(row for row in summaries if row["experiment"] == "MEKi0.1")
        self.assertEqual(summary["numeric_available_n"], 2)
        self.assertEqual(summary["median_t"], 1.0)
        self.assertIsNone(summary["fraction_negative"])
        self.assertIsNone(summary["fraction_positive"])

    def test_validated_rank_multiplier_is_applied_before_signed_correlation(self):
        experiments = {
            "g1": {"by_id": {
                "ENSG1": {"t": -2.0, "logFC": -1.0},
                "ENSG2": {"t": 1.0, "logFC": 0.5},
            }},
            "MEKi0.1": {"by_id": {
                "ENSG1": {"value": 3.0},
                "ENSG2": {"value": -4.0},
            }},
        }
        orientation = BENCHMARK.default_mek_orientation()
        orientation["MEKi0.1"] = {
            "status": "validated",
            "sign_multiplier": -1,
            "evidence": "test fixture",
            "directional_metrics_available": True,
        }
        row = BENCHMARK.build_comparison_row("g1", "MEKi0.1", experiments, orientation)
        self.assertEqual(row["spearman_signed_t"], 1.0)
        self.assertEqual(row["direction_agreement_t"], 1.0)

    def test_zero_statistics_are_excluded_from_sign_denominator(self):
        summary = BENCHMARK.sign_concordance([0.0, -1.0, 1.0, 0.0], [1.0, -1.0, 0.0, 0.0])
        self.assertEqual(summary["finite_pairs"], 4)
        self.assertEqual(summary["informative_pairs"], 1)
        self.assertEqual(summary["same_sign_n"], 1)
        self.assertEqual(summary["opposite_sign_n"], 0)
        self.assertEqual(summary["both_zero_n"], 1)
        self.assertEqual(summary["one_zero_n"], 2)
        self.assertEqual(summary["fraction"], 1.0)

    def test_constant_signed_vector_has_no_spearman_value(self):
        self.assertIsNone(BENCHMARK.spearman_signed([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]))
        self.assertEqual(BENCHMARK.spearman_signed([1.0, 2.0, 3.0], [3.0, 2.0, 1.0]), -1.0)


if __name__ == "__main__":
    unittest.main()
