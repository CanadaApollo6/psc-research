"""Synthetic-contract tests for the frozen ETS2 atlas program helpers."""

from __future__ import annotations

import io
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import sparse


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.analyze_ets2_program import (  # noqa: E402
    AUTHOR_CELL_TYPES,
    DEFAULT_MIN_PAIRED_DONORS,
    DEFAULT_REPLICATES,
    PRIMARY_PROGRAM_ID,
    aggregate_all_features,
    build_feature_map,
    build_correlations,
    build_paired_contrasts,
    build_program_mappings,
    expression_bins,
    map_source_program,
    matched_seed,
    read_source_column,
    sample_expression_matched_sets,
)


def synthetic_adata(*, normalized_same_as_raw: bool = False) -> SimpleNamespace:
    """Return a tiny AnnData-like object with repeated donor/library cells."""

    # The first four cells are one donor/cell type split across two libraries;
    # the remaining two are a second donor.  Feature 2 is not a target gene,
    # so it checks that the denominator uses every raw feature.
    counts = np.asarray(
        [
            [4, 1, 10, 0],
            [2, 0, 5, 0],
            [3, 1, 0, 0],
            [1, 2, 5, 0],
            [10, 0, 0, 2],
            [10, 0, 0, 2],
        ],
        dtype=np.int64,
    )
    obs = pd.DataFrame(
        {
            "donor_id": ["D1", "D1", "D1", "D1", "D2", "D2"],
            "author_cell_type": ["Monocyte", "Monocyte", "Monocyte", "Monocyte", "Kupffer", "Kupffer"],
            "cell_type": ["myeloid"] * 6,
            "disease": ["PSC"] * 6,
            "assay": ["RNA"] * 6,
            "suspension_type": ["single-cell"] * 6,
            "library_uuid": ["L1", "L1", "L2", "L2", "L3", "L3"],
        },
        index=[f"cell-{i}" for i in range(len(counts))],
    )
    raw_matrix = sparse.csr_matrix(counts)
    normalized = raw_matrix if normalized_same_as_raw else counts.astype(float) + 100.0
    raw_var = pd.DataFrame(
        {"feature_name": ["G1", "G2", "G3", "ETS2"]},
        index=["ENSG00000000001.1", "ENSG00000000002", "ENSG00000000003", "ENSG00000157557"],
    )
    raw = SimpleNamespace(X=raw_matrix, var=raw_var, n_vars=counts.shape[1])
    return SimpleNamespace(
        n_obs=len(counts),
        n_vars=counts.shape[1],
        shape=counts.shape,
        obs=obs,
        obs_names=obs.index,
        raw=raw,
        X=normalized,
    )


class TestAggregation(unittest.TestCase):
    def test_reuses_donor_across_libraries_and_uses_all_raw_features(self) -> None:
        result = aggregate_all_features(synthetic_adata(), "sc", dataset_id="dataset")
        self.assertEqual(len(result.pseudobulk), 2)
        monocyte = result.pseudobulk[result.pseudobulk.author_cell_type.eq("Monocyte")].iloc[0]
        self.assertEqual(int(monocyte.cell_count), 4)
        self.assertEqual(int(monocyte.library_count), 2)
        self.assertEqual(int(monocyte.all_feature_umis), 34)
        self.assertEqual(int(result.count_values[0].sum()), 34)
        # Raw counts sum to 10/34 CPM for G1 in the first group; the large
        # normalized X values must not leak into the denominator.
        self.assertAlmostEqual(float(result.log_values[0, 0]), np.log2(10 / 34 * 1_000_000 + 1))

    def test_normalized_matrix_is_not_accepted_as_raw(self) -> None:
        with self.assertRaises(ValueError):
            aggregate_all_features(synthetic_adata(normalized_same_as_raw=True), "sc", dataset_id="dataset")

    def test_sparse_float_counts_are_cast_before_exact_group_sum(self) -> None:
        adata = synthetic_adata()
        # Float storage is allowed when every stored value is an integer.  The
        # first value is just above float32's exact-integer boundary; summing
        # as float32 before casting would lose the trailing one.
        counts = sparse.csr_matrix(np.asarray([[16_777_216.0, 1.0, 0.0, 0.0]], dtype=np.float32))
        adata.n_obs = 1
        adata.shape = counts.shape
        adata.obs = adata.obs.iloc[:1].copy()
        adata.obs_names = adata.obs.index
        adata.raw = SimpleNamespace(X=counts, var=adata.raw.var.copy(), n_vars=counts.shape[1])
        adata.X = np.zeros(counts.shape, dtype=np.float64)
        result = aggregate_all_features(adata, "sc", dataset_id="dataset")
        self.assertEqual(int(result.pseudobulk.iloc[0].all_feature_umis), 16_777_217)


class TestMapping(unittest.TestCase):
    def test_exact_id_and_name_mapping_keep_ambiguity_visible(self) -> None:
        raw_var = pd.DataFrame(
            {"feature_name": ["SYM", "SYM", "ETS2", "OTHER"]},
            index=["ENSG00000000001.2", "ENSG00000000002", "ENSG00000157557", "ENSG00000000004"],
        )
        features = build_feature_map(raw_var, 4)
        table, mapped, inventory = map_source_program(
            "demo", ["ENSG00000000001", "SYM", "MISSING", "ENSG00000157557"], features, dataset_key="sc"
        )
        statuses = dict(zip(table.source_value, table.mapping_status))
        self.assertEqual(statuses["ENSG00000000001"], "mapped")
        self.assertEqual(statuses["SYM"], "ambiguous_feature_name")
        self.assertEqual(statuses["MISSING"], "missing_feature_name")
        self.assertEqual(statuses["ENSG00000157557"], "excluded_target_gene")
        self.assertEqual(inventory["mapped_unique_genes"], 1)
        self.assertNotIn(2, mapped)

    def test_nonoverlap_uses_mapped_membership_not_source_text(self) -> None:
        raw_var = pd.DataFrame(
            {"feature_name": ["A", "B", "C", "D", "ETS2"]},
            index=[f"ENSG0000000000{i}" for i in range(1, 5)] + ["ENSG00000157557"],
        )
        features = build_feature_map(raw_var, 5)
        source = {
            "ets2_g1_dn": ["ENSG00000000001", "ENSG00000000002", "ENSG00000157557"],
            "ets2_g1_up": ["ENSG00000000003"],
            "ets2_g2_dn": ["ENSG00000000004"],
            "chr21_dn": ["ENSG00000000004"],
            "inflammation": ["A"],
            "interferon_gamma": ["B"],
            "oxidative_stress": ["C"],
            "apoptotic_signaling": ["D"],
        }
        plan = {
            "programs": [
                {"id": "ets2_g1_dn", "role": "primary"},
                {"id": "ets2_g1_up", "role": "companion"},
                {"id": "ets2_g2_dn", "role": "companion"},
                {"id": "chr21_dn", "role": "companion"},
                {"id": "inflammation", "role": "comparison"},
                {"id": "interferon_gamma", "role": "comparison"},
                {"id": "oxidative_stress", "role": "comparison"},
                {"id": "apoptotic_signaling", "role": "comparison"},
                {"id": "ets2_g1_dn_without_comparators", "role": "sensitivity"},
            ]
        }
        membership, inventory, mapped = build_program_mappings(plan, source, {"sc": features})
        derived = inventory[inventory.program_id.eq("ets2_g1_dn_without_comparators")].iloc[0]
        # A and B overlap with the four comparator memberships and are removed;
        # the target set is therefore empty and cannot pass its gate.
        self.assertEqual(int(derived.derived_nonoverlap_mapped_genes), 0)
        self.assertEqual(mapped["sc"]["ets2_g1_dn_without_comparators"], set())
        derived_rows = membership[(membership.program_id == "ets2_g1_dn_without_comparators") & membership.source_value.isin(["ENSG00000000001", "ENSG00000000002"])]
        self.assertTrue(derived_rows.mapping_status.isin(["mapped_comparator_overlap", "mapped_nonoverlap"]).all())


class TestMatching(unittest.TestCase):
    def test_expression_bin_zero_is_low_reference_abundance(self) -> None:
        feature_bin, bins = expression_bins([3.0, 1.0, 2.0, 1.0], ["B", "A", "C", "D"], bins=2)
        self.assertEqual(list(bins[0]), [1, 3])
        self.assertEqual(list(bins[1]), [2, 0])
        np.testing.assert_array_equal(feature_bin, np.asarray([1, 0, 1, 0]))

    def test_matching_excludes_program_union_and_preserves_bin_counts(self) -> None:
        bins = [np.asarray([0, 1, 2]), np.asarray([3, 4, 5]), np.asarray([6, 7, 8])]
        target = [0, 4, 6]
        excluded = {0, 4, 6, 1, 5, 7}
        first = sample_expression_matched_sets(bins, target, excluded, replicates=20, seed=matched_seed(1, "sc", "demo"))
        second = sample_expression_matched_sets(bins, target, excluded, replicates=20, seed=matched_seed(1, "sc", "demo"))
        self.assertTrue(first.available)
        np.testing.assert_array_equal(first.selections, second.selections)
        self.assertEqual(first.target_bin_counts, {0: 1, 1: 1, 2: 1})
        self.assertTrue(all(len(set(row)) == 3 for row in first.selections))
        for row in first.selections:
            self.assertIn(int(row[0]), {2})
            self.assertIn(int(row[1]), {3})
            self.assertIn(int(row[2]), {8})

    def test_insufficient_controls_are_unavailable_without_fallback(self) -> None:
        bins = [np.asarray([0, 1]), np.asarray([2, 3])]
        result = sample_expression_matched_sets(bins, [0, 1], {0, 1, 2}, replicates=DEFAULT_REPLICATES, seed=3)
        self.assertFalse(result.available)
        self.assertEqual(result.reason, "insufficient_control_pool")
        self.assertEqual(result.underflow_bins, (0,))
        self.assertEqual(result.selections.shape, (0, 2))


class TestPairedContrast(unittest.TestCase):
    def _scores(self, donors: int) -> pd.DataFrame:
        rows = []
        for donor_index in range(donors):
            for population, score in (("ActMac", float(donor_index + 1)), ("Kupffer", 0.5)):
                rows.append(
                    {
                        "dataset_key": "sc",
                        "dataset_id": "dataset",
                        "modality": "sc",
                        "suspension_type": "single-cell",
                        "assay": "RNA",
                        "donor_id": f"D{donor_index}",
                        "program_id": "ets2_g1_dn",
                        "threshold": 20,
                        "author_cell_type": population,
                        "program_score": score,
                        "score_status": "evaluable_descriptive",
                    }
                )
        return pd.DataFrame(rows)

    def test_paired_donor_gate_and_leave_one_out_gate(self) -> None:
        plan = {
            "contrasts": [{"left": "ActMac", "right": "Kupffer", "role": "primary"}],
            "minimum_paired_donors": DEFAULT_MIN_PAIRED_DONORS,
        }
        three = build_paired_contrasts(self._scores(3), {}, plan).iloc[0]
        self.assertEqual(three.status, "evaluable_descriptive")
        self.assertEqual(int(three.paired_donors), 3)
        self.assertEqual(int(three.loo_evaluable_count), 0)
        two = build_paired_contrasts(self._scores(2), {}, plan).iloc[0]
        self.assertEqual(two.status, "insufficient_paired_donors")
        self.assertTrue(pd.isna(two.median_delta))

    def test_correlations_keep_missing_author_populations_explicit(self) -> None:
        rows = []
        all_programs = (PRIMARY_PROGRAM_ID,) + ("inflammation", "interferon_gamma", "oxidative_stress", "apoptotic_signaling")
        for donor_index in range(4):
            for population in ("Monocyte", "Kupffer"):
                for program_index, program_id in enumerate(all_programs):
                    rows.append(
                        {
                            "dataset_key": "sn",
                            "dataset_id": "dataset",
                            "modality": "sn",
                            "suspension_type": "single-nucleus",
                            "assay": "RNA",
                            "donor_id": f"D{donor_index}",
                            "author_cell_type": population,
                            "threshold": 20,
                            "program_id": program_id,
                            "program_score": float(donor_index + program_index + 1),
                            "ets2_log2_cpm1": float(donor_index + 1),
                            "score_status": "evaluable_descriptive",
                        }
                    )
        result = build_correlations(pd.DataFrame(rows), {"minimum_correlation_donors": 4})
        self.assertEqual(len(result), 4 * 5)
        missing = result[result.author_cell_type.isin(["ActMac", "LAM-like"])]
        self.assertEqual(len(missing), 2 * 5)
        self.assertTrue((missing.status == "no_observed_context").all())


class TestSourceReading(unittest.TestCase):
    def test_reads_the_named_column_without_inferring_direction(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "source.zip"
            payload = b"A,B\nENSG00000000001,ENSG00000000002\n,ENSG00000000003\n"
            with zipfile.ZipFile(path, "w") as archive:
                archive.writestr("prefix/sets.csv", payload)
            values, metadata = read_source_column(path, "sets.csv", "B")
            self.assertEqual(values, ["ENSG00000000002", "ENSG00000000003"])
            self.assertEqual(metadata["bytes"], len(payload))
            self.assertEqual(len(metadata["sha256"]), 64)


if __name__ == "__main__":
    unittest.main()
