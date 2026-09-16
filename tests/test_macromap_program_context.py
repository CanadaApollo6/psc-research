"""Synthetic and source-pin tests. Never score real MacroMap expression."""
from __future__ import annotations

import copy
import gzip
import hashlib
import json
import math
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts import analyze_macromap_program_context as m


def metadata_fixture(runs=5, lines_per_run=2):
    rows = []
    for protocol in m.PROTOCOLS:
        for run in range(runs):
            for line in range(lines_per_run):
                identity = f"FD_{protocol}_{run}_{line}"
                for time in m.TIMES:
                    for stimulus in m.STIMULI + ("Ctrl", "PIC"):
                        condition = f"{stimulus}_{time}"
                        rows.append({"SampleID": identity + "_" + condition,
                                     "Stimulus_Hours": condition, "RunID": f"{protocol}-R{run}",
                                     "Library_prep": protocol, "HipsciID": f"H-{identity}",
                                     "Date_EB_formation": "same-date", "Sex": "Female"})
    return rows, [r["SampleID"] + "_" + r["RunID"] for r in rows]


def validated_fixture(**kwargs):
    rows, header = metadata_fixture(**kwargs)
    return m.validate_samples(rows, header)


def prediction_fixture():
    rows = []
    for pindex, protocol in enumerate(m.PROTOCOLS):
        for run in range(3):
            # Unequal line counts within runs and protocols test line weighting.
            for donor in range(run + 1 + pindex):
                line = f"{protocol}-{run}-{donor}"
                for time in m.TIMES:
                    stimuli = m.STIMULI if donor == 0 else m.STIMULI[:2]
                    for j, stimulus in enumerate(stimuli):
                        base = 3.0 + run / 10 + j / 100
                        gain = 0.2 + pindex * 0.1 + run * 0.05 + donor * 0.01 + time / 1000
                        rows.append({"protocol": protocol, "run": str(run), "line": line,
                                     "time": time, "stimulus": stimulus, "status": "available",
                                     "baseline_loss": base, "extended_loss": base - gain})
    return rows


class SourceIdentityTests(unittest.TestCase):
    def test_pins_fail_on_mutation_and_absence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "source").write_bytes(b"source")
            pins = {"source": {"bytes": 6, "sha256": hashlib.sha256(b"source").hexdigest()}}
            self.assertEqual(m.verify_pins(root, pins), pins)
            (root / "source").write_bytes(b"Source")
            with self.assertRaises(ValueError):
                m.verify_pins(root, pins)
            (root / "source").unlink()
            with self.assertRaises(FileNotFoundError):
                m.verify_pins(root, pins)

    def test_public_source_pin_invariants(self):
        # Only small preserved scientific files; raw caches need not exist on CI.
        for key in ("config/ets2-program-plan.json", "config/ets2-validation-amended-rules.json",
                    "data/derived/ets2-program-membership.csv", "scripts/analyze_ets2_program.py"):
            self.assertEqual(m.sha256(ROOT / key), m.SOURCE_PINS[key]["sha256"])
        self.assertEqual(m.SOURCE_PINS[m.RAW + "/MacroMap_raw_expression.txt.gz"]["bytes"], 174122135)

    def test_preserved_full_members_not_old_mapped_subset(self):
        programs = m.retained_source_programs(ROOT)
        self.assertEqual(len(programs["ets2_g1_dn"]), 928)
        self.assertIn(m.original.ETS2_STABLE_ID, programs["ets2_g1_dn"])
        self.assertEqual(len(programs["ets2_g1_up"]), 668)
        self.assertEqual(len(programs["inflammation"]), 815)
        self.assertEqual(tuple(programs), m.original.ORIGINAL_PROGRAM_IDS)

    def test_version_suffix_par_y_ambiguity_and_zero_universe(self):
        universe = [{"Geneid": "ENSG000001.1"}, {"Geneid": "ENSG000002.4"}, {"Geneid": "ENSG000003.2"}]
        annotation = ["chr1 1 10 + ENSG000001.1 SAME", "chr1 20 30 + ENSG000002.4 SAME",
                      "chr1 40 50 - ENSG000003.2 UNIQUE", "chrY 1 10 + ENSG000001.1_PAR_Y SAME"]
        features, absent = m.reference_features(universe, annotation)
        self.assertEqual(absent, ["ENSG000001.1_PAR_Y"])
        self.assertEqual(len(features), 3)
        audit, indices, inventory = m.original.map_source_program(
            "toy", ["SAME", "UNIQUE", "ENSG000001.7"], features, dataset_key="synthetic")
        self.assertEqual(indices, {0, 2})
        self.assertEqual(inventory["ambiguous_source_genes"], 1)
        self.assertFalse(inventory["mapping_gate_pass"])
        self.assertEqual(audit.iloc[0].mapping_status, "ambiguous_feature_name")
        with self.assertRaises(ValueError):
            m.reference_features(universe + [{"Geneid": "ENSG000001.2"}],
                                 annotation + ["chr1 1 10 + ENSG000001.2 OTHER"])

    def test_mapping_keeps_parent_gate_after_overlap_removal(self):
        import pandas as pd
        ids = [f"ENSG{i:011d}" for i in range(40)]
        features = m.original.build_feature_map(pd.DataFrame({"feature_name": ids}, index=ids), 40)
        plan = json.loads((ROOT / "config/ets2-program-plan.json").read_text())
        source = {p: ids[:20] for p in m.original.ORIGINAL_PROGRAM_IDS}
        for p in m.COMPARATORS:
            source[p] = ids[10:30]
        _, inventory, mapped = m.original.build_program_mappings(plan, source, {"toy": features})
        nonoverlap = inventory[inventory.program_id.eq("ets2_g1_dn_without_comparators")].iloc[0]
        self.assertTrue(nonoverlap.mapping_gate_pass)
        self.assertEqual(nonoverlap.mapping_fraction, 0.5)
        self.assertEqual(nonoverlap.parent_mapping_fraction, 1.0)
        self.assertEqual(len(mapped["toy"]["ets2_g1_dn_without_comparators"]), 10)


class PairAndFoldTests(unittest.TestCase):
    def test_exact_header_join_is_order_independent(self):
        rows, header = metadata_fixture()
        samples = m.validate_samples(rows[::-1], header[::-1])
        self.assertEqual([r["sample_id"] for r in samples], header[::-1])
        self.assertEqual([r["sample_index"] for r in samples], list(range(len(header))))
        with self.assertRaises(ValueError):
            m.validate_samples(rows, header[:-1])

    def test_csv_pairs_restore_boolean_and_time_index_types(self):
        pairs = m.pair_registry(validated_fixture())
        pairs[0]["identity_mapped"] = False
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "pairs.csv"
            m.write_csv(path, pairs)
            restored = m.read_pairs(path)
            self.assertEqual(restored, pairs)
            self.assertIs(restored[0]["identity_mapped"], False)
            self.assertIsInstance(restored[0]["time"], int)
        with self.assertRaises(ValueError):
            m._bool_field("false-or-missing")

    def test_pairs_require_same_time_and_pic_is_unavailable(self):
        samples = validated_fixture()
        line = samples[0]["line"]
        samples = [r for r in samples if not (r["line"] == line and r["condition"] == "Ctrl_6")]
        pairs = m.pair_registry(samples)
        chosen = [r for r in pairs if r["line"] == line]
        self.assertTrue(all(r["status"] == "unavailable_mock_identity" for r in chosen if r["stimulus"] == "PIC"))
        self.assertTrue(all(r["status"] == "missing_time_matched_control" for r in chosen if r["time"] == 6 and r["stimulus"] != "PIC"))
        self.assertTrue(all(r["status"] == "available" for r in chosen if r["time"] == 24 and r["stimulus"] != "PIC"))
        self.assertNotIn(line, {r["line"] for r in m.prediction_cohort(pairs)})

    def test_nomatch_labels_remain_separate(self):
        rows, header = metadata_fixture()
        chosen = {rows[0]["HipsciID"], rows[22]["HipsciID"]}
        for row in rows:
            if row["HipsciID"] in chosen:
                row["HipsciID"] = "NOMATCH"
        samples = m.validate_samples(rows, header)
        unmapped = {r["line"] for r in samples if not r["identity_mapped"]}
        self.assertEqual(len(unmapped), 2)
        primary = m.pair_registry(samples)
        self.assertTrue(all(r["status"] == "excluded_NOMATCH_identity" for r in primary if r["line"] in unmapped))
        sensitivity = m.pair_registry(samples, include_nomatch=True)
        self.assertEqual(len({r["line"] for r in sensitivity if r["status"] == "available"}), 20)

    def test_static_culture_run_protocol_duplicates_fail_closed(self):
        rows, header = metadata_fixture()
        for field in ("Date_EB_formation", "RunID", "Library_prep", "HipsciID"):
            altered = copy.deepcopy(rows)
            altered[0][field] = "different"
            altered_header = [r["SampleID"] + "_" + r["RunID"] for r in altered]
            with self.assertRaises(ValueError):
                m.validate_samples(altered, altered_header)
        with self.assertRaises(ValueError):
            m.validate_samples(rows + [rows[0]], header + [header[0]])
        altered = copy.deepcopy(rows)
        first, second = altered[0]["HipsciID"], altered[22]["HipsciID"]
        for row in altered:
            if row["HipsciID"] == second:
                row["HipsciID"] = first
        with self.assertRaises(ValueError):
            m.validate_samples(altered, header)

    def test_fixed_run_folds_no_line_or_protocol_leak(self):
        samples = validated_fixture()
        assignment = m.build_folds(samples)
        self.assertEqual(assignment, m.build_folds(samples[::-1]))
        self.assertEqual(set(assignment.values()), set(range(5)))
        pairs = m.prediction_cohort(m.pair_registry(samples))
        for p in m.PROTOCOLS:
            for t in m.TIMES:
                for f in range(5):
                    train, test, gate = m.split_pairs(pairs, assignment, p, t, f)
                    self.assertEqual(gate["status"], "available")
                    self.assertFalse({r["run"] for r in train} & {r["run"] for r in test})
                    self.assertFalse({r["line"] for r in train} & {r["line"] for r in test})
                    self.assertEqual(gate["missing_test_stimuli"], [])
        with self.assertRaises(ValueError):
            m.build_folds(validated_fixture(runs=4))

    def test_fold_class_gate_never_remakes_folds(self):
        samples = validated_fixture()
        assignment = m.build_folds(samples)
        pairs = m.prediction_cohort(m.pair_registry(samples))
        p = m.PROTOCOLS[0]
        keep = [r for r in pairs if r["stimulus"] != "IFNG" or assignment[(r["protocol"], r["run"])] == 0]
        _, _, gate = m.split_pairs(keep, assignment, p, 6, 0)
        self.assertEqual(gate["status"], "insufficient_training_class_lines")
        self.assertEqual(assignment, m.build_folds(samples))


class TrainingOnlyAndScoreTests(unittest.TestCase):
    def test_training_controls_only_reference_shared_controls_once(self):
        samples = validated_fixture()
        pairs = m.pair_registry(samples)
        assignment = m.build_folds(samples)
        values = np.arange(3 * len(samples), dtype=float).reshape(3, len(samples))
        p, time, fold = m.PROTOCOLS[0], 6, 0
        reference, audit = m.training_control_reference(values, pairs, assignment, p, time, fold)
        self.assertEqual(audit["reference_lines"], 8)
        self.assertEqual(audit["reference_runs"], 4)
        assert_allclose(reference, values[:, audit["control_indices"]].mean(axis=1))
        corrupted = values.copy()
        other = np.setdiff1d(np.arange(len(samples)), audit["control_indices"])
        corrupted[:, other] = np.nan  # includes all treatments and held-out controls
        again, audit2 = m.training_control_reference(corrupted, pairs, assignment, p, time, fold)
        assert_array_equal(reference, again)
        self.assertEqual(audit, audit2)

    def test_matched_weight_algebra_exact_random_sets(self):
        reference = np.linspace(0, 3, 100)
        ids = [f"ENSG{i:011d}" for i in range(100)]
        target = [3, 25, 45, 70, 98]
        excluded = target + [2, 22]
        info = m.matched_weights(reference, ids, {"toy": target}, excluded, "P", 6, 0, replicates=30, bins=5)["toy"]
        self.assertTrue(info["available"])
        self.assertEqual(info["background_weights"][excluded].sum(), 0)
        _, bins = m.original.expression_bins(reference, ids, bins=5)
        sets = m.original.sample_expression_matched_sets(bins, target, excluded, replicates=30, seed=info["seed"])
        x = np.random.default_rng(6).uniform(size=(100, 7))
        scalar_raw = np.array([x[target, j].mean() for j in range(7)])
        scalar_background = np.array([np.mean([x[indices, j].mean() for indices in sets.selections]) for j in range(7)])
        scored = m.score_with_weights(x, info, list(range(7)))
        assert_allclose(scored["raw"], scalar_raw, atol=1e-14)
        assert_allclose(scored["background"], scalar_background, atol=1e-14)
        assert_allclose(scored["adjusted"], scalar_raw - scalar_background, atol=1e-14)
        again = m.matched_weights(reference, ids, {"toy": target}, excluded, "P", 6, 0, replicates=30, bins=5)["toy"]
        self.assertEqual(info["selection_sha256"], again["selection_sha256"])

    def test_matching_requires_target_exclusion_and_positive_draws(self):
        ids = [f"ENSG{i:011d}" for i in range(20)]
        with self.assertRaises(ValueError):
            m.matched_weights(np.zeros(20), ids, {"toy": [0]}, [], "P", 6, 0, bins=2)
        with self.assertRaises(ValueError):
            m.matched_weights(np.zeros(20), ids, {"toy": [0]}, [0], "P", 6, 0, bins=2, replicates=0)

    def test_matching_underflow_unavailable_not_zero(self):
        ids = [f"ENSG{i:011d}" for i in range(6)]
        info = m.matched_weights(np.zeros(6), ids, {"toy": [0, 1]}, list(range(6)), "P", 6, 0,
                                 replicates=3, bins=2)["toy"]
        self.assertFalse(info["available"])
        self.assertEqual(info["reason"], "insufficient_control_pool")
        scored = m.score_with_weights(np.zeros((6, 1)), info, [0])
        self.assertIsNone(scored["adjusted"])

    def test_equal_line_weights_scaler_and_no_test_fitting(self):
        weights = m.line_weights(["a", "a", "b"])
        assert_allclose(weights, [0.25, 0.25, 0.5])
        x = np.array([[0.0, 1], [2.0, 1], [4.0, 1]])
        scaler = m.fit_scaler(x, ["a", "a", "b"])
        assert_allclose(scaler.mean, [2.5, 1])
        assert_allclose(scaler.scale, [math.sqrt(2.75), 1])
        assert_array_equal(scaler.constant, [False, True])
        before = scaler.mean.copy()
        scaler.transform(np.array([[1e12, 90]]))
        assert_array_equal(before, scaler.mean)


class RidgeTests(unittest.TestCase):
    def test_weighted_objective_gradient_finite_difference(self):
        rng = np.random.default_rng(16)
        x = rng.normal(size=(7, 3))
        y = np.array([0, 1, 2, 0, 1, 2, 2])
        w = m.line_weights(["a", "a", "b", "c", "d", "e", "e"])
        theta = rng.normal(size=12)
        objective, grad = m.ridge_objective(theta, x, y, w, 3)
        numerical = np.zeros(12)
        for i in range(12):
            plus, minus = theta.copy(), theta.copy()
            plus[i] += 1e-5
            minus[i] -= 1e-5
            numerical[i] = (m.ridge_objective(plus, x, y, w, 3)[0] -
                            m.ridge_objective(minus, x, y, w, 3)[0]) / 2e-5
        assert_allclose(grad, numerical, atol=2e-10)
        self.assertTrue(math.isfinite(objective))

    def test_binary_ridge_matches_independent_scalar_root(self):
        x = np.array([[-1.0], [1.0]])
        fit = m.fit_multinomial_ridge(x, [0, 1], ["a", "b"], n_classes=2)
        self.assertEqual(fit.status, "available")
        # log(1+exp(-b)) + lambda*b^2/4 for symmetric class slopes +/-b/2.
        slope = brentq(lambda b: -1 / (1 + math.exp(b)) + m.RIDGE_LAMBDA * b / 2, 0, 20)
        self.assertAlmostEqual(fit.coefficients[1, 1] - fit.coefficients[1, 0], slope, places=6)
        assert_allclose(fit.coefficients[0], 0, atol=1e-8)
        assert_allclose(np.exp(fit.log_proba(x)).sum(axis=1), 1, atol=1e-14)

    def test_intercept_fit_equal_line_empirical_priors_and_class_symmetry(self):
        x = np.zeros((4, 2))
        fit = m.fit_multinomial_ridge(x, [0, 0, 0, 1], ["a", "a", "a", "b"], n_classes=2)
        self.assertEqual(fit.status, "available")
        assert_allclose(np.exp(fit.log_proba(x)), np.full((4, 2), 0.5), atol=1e-8)
        fit = m.fit_multinomial_ridge(x, [0, 0, 0, 1], ["a", "b", "c", "d"], n_classes=2)
        assert_allclose(np.exp(fit.log_proba(x))[0], [0.75, 0.25], atol=1e-6)
        swapped = m.fit_multinomial_ridge(x, [1, 1, 1, 0], ["a", "b", "c", "d"], n_classes=2)
        assert_allclose(fit.log_proba(x), swapped.log_proba(x)[:, ::-1], atol=1e-7)

    def test_missing_stimulus_and_nonconvergence_return_no_model(self):
        missing = m.fit_multinomial_ridge(np.zeros((3, 2)), [0, 0, 1], ["a", "b", "c"], n_classes=3)
        self.assertEqual(missing.status, "missing_training_stimulus")
        with self.assertRaises(ValueError):
            missing.log_proba(np.zeros((1, 2)))
        failed = m.fit_multinomial_ridge(np.array([[-1.0], [1.0]]), [0, 1], ["a", "b"], n_classes=2, max_iter=0)
        self.assertEqual(failed.status, "optimizer_failed")
        self.assertIsNone(failed.coefficients)
        invalid = m.fit_multinomial_ridge(np.full((2, 1), np.nan), [0, 1], ["a", "b"], n_classes=2)
        self.assertEqual(invalid.status, "invalid_training_data")

    def test_nine_class_paired_fit_ties_missing_values_and_leakage(self):
        train, test = [], []
        for group, bucket in (("train", train), ("test", test)):
            for i in range(3):
                for stimulus in m.STIMULI:
                    bucket.append({"line": f"{group}-{i}", "run": f"{group}-r{i}", "time": 6,
                                   "protocol": m.PROTOCOLS[0], "stimulus": stimulus, "status": "available"})
        tx = np.zeros((27, 5))
        ex = np.zeros((27, 5))
        result = m.heldout_increment(tx, ex, train, test)
        self.assertEqual(result["status"], "available")
        for row in result["predictions"]:
            self.assertAlmostEqual(row["baseline_loss"], math.log(9))
            self.assertEqual(row["baseline_class"], "CIL")
            self.assertAlmostEqual(row["gain"], 0)
        ex[0, 4] = np.nan
        self.assertEqual(m.heldout_increment(tx, ex, train, test)["status"], "missing_predictor")
        with self.assertRaises(ValueError):
            m.heldout_increment(tx, tx, train, train)

    def test_nonfinite_prediction_is_not_clipped_or_silently_retained(self):
        fit = m.RidgeFit("available", np.array([[0., 0.], [1e308, -1e308]]), 2)
        with self.assertRaises(ValueError):
            fit.log_proba(np.array([[1e308]]))

    def test_log_softmax_stays_finite_for_extreme_heldout_features(self):
        fit = m.fit_multinomial_ridge(np.array([[-1.], [1.]]), [0, 1], ["a", "b"], n_classes=2)
        result = fit.log_proba(np.array([[1e6], [-1e6]]))
        self.assertTrue(np.isfinite(result).all())
        self.assertTrue(np.all(result <= 0))


class ClusteredMetricTests(unittest.TestCase):
    def test_loss_aggregation_equal_lines_times_once_and_protocol_line_shares(self):
        rows = prediction_fixture()
        result = m.aggregate_prediction(rows, expected_pairs=rows)
        self.assertEqual(result["status"], "available")
        lines = m.line_time_losses(rows)
        self.assertEqual(result["line_counts"], {m.PROTOCOLS[0]: 6, m.PROTOCOLS[1]: 9})
        expected = np.mean([r["gain"] for r in lines])  # both times for every line
        self.assertAlmostEqual(result["primary_gain"], expected)
        naive_row_mean = np.mean([r["baseline_loss"] - r["extended_loss"] for r in rows])
        self.assertNotAlmostEqual(result["primary_gain"], naive_row_mean, places=5)

    def test_missing_failed_or_duplicate_predictions_not_dropped(self):
        rows = prediction_fixture()
        result = m.aggregate_prediction(rows[:-1], expected_pairs=rows)
        self.assertEqual(result["status"], "incomplete_predictions")
        corrupt = copy.deepcopy(rows)
        corrupt[0]["baseline_loss"] = np.nan
        with self.assertRaises(ValueError):
            m.aggregate_prediction(corrupt, expected_pairs=rows)
        with self.assertRaises(ValueError):
            m.line_time_losses(rows + [rows[0]])
        missing_class = [r for r in rows if r["stimulus"] != "IFNG"]
        self.assertEqual(m.aggregate_prediction(missing_class, expected_pairs=missing_class)["status"], "missing_evaluation_stimulus")
        missing_time = [r for r in rows if not (r["line"] == rows[0]["line"] and r["time"] == 24)]
        self.assertEqual(m.aggregate_prediction(missing_time, expected_pairs=missing_time)["status"], "unpaired_times")

    def test_run_bootstrap_matches_explicit_cluster_replication(self):
        rows = prediction_fixture()
        result = m.cluster_bootstrap_gain(rows, expected_pairs=rows, replicates=25, seed=4)
        self.assertEqual(result["interval_status"], "fixed_prediction_conditional_only")
        rng = np.random.Generator(np.random.PCG64(4))
        manual = []
        for _ in range(25):
            multiplicity = {}
            for p in m.PROTOCOLS:
                runs = sorted({r["run"] for r in rows if r["protocol"] == p})
                selected = rng.integers(len(runs), size=len(runs))
                counts = Counter(selected)
                multiplicity.update({(p, run): counts[i] for i, run in enumerate(runs)})
            manual.append(m.aggregate_prediction(rows, expected_pairs=rows, multiplicity=multiplicity,
                                                protocol_weights=result["protocol_weights"])["primary_gain"])
        assert_allclose(result["interval"], np.quantile(manual, [0.025, 0.975]), atol=1e-14)
        # Row order must not determine random cluster ordering or weights.
        shuffled = m.cluster_bootstrap_gain(rows[::-1], expected_pairs=rows[::-1], replicates=25, seed=4)
        assert_allclose(result["interval"], shuffled["interval"], atol=1e-14)

    def test_too_few_runs_not_replaced_by_many_treatment_rows(self):
        rows = prediction_fixture()
        for row in rows:
            row["run"] = "one"
        result = m.cluster_bootstrap_gain(rows, expected_pairs=rows, replicates=10)
        self.assertEqual(result["interval_status"], "insufficient_runs")
        self.assertIsNone(result["interval"])

    def test_response_summary_preserves_old_gates_and_uses_clusters(self):
        result = m.clustered_response_summary([1, 1, 1, 9], ["a", "b", "c", "d"], ["A", "A", "A", "B"], replicates=10)
        self.assertEqual(result["median"], 1)
        self.assertEqual(result["mean"], 3)
        self.assertEqual(result["interval_status"], "insufficient_runs")
        self.assertEqual(result["leave_one_line_out_median_range"], [1, 1])
        with self.assertRaises(ValueError):
            m.clustered_response_summary([1, 2, 3], ["a", "a", "b"], ["A", "A", "B"])
        invalid = m.clustered_response_summary([1, np.nan, 3], ["a", "b", "c"], ["A", "B", "C"])
        self.assertEqual(invalid["status"], "missing_response")
        result = m.clustered_response_summary([1, 2, 3, 9], ["a", "b", "c", "d"], ["A", "A", "B", "C"], replicates=20)
        self.assertEqual(result["interval_status"], "pointwise_fixed_score_conditional_only")
        self.assertEqual(len(result["median_interval"]), 2)


class OutputGuardTests(unittest.TestCase):
    def test_private_anchor_symlink_and_reused_run_directory_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "work").mkdir()
            (root / "reports").mkdir()
            (root / "work/macromap-design").symlink_to(root / "reports", target_is_directory=True)
            with self.assertRaises(ValueError):
                m.fresh_private_run(root, root / "work/macromap-design/private-design")
            self.assertFalse((root / "reports/private-design").exists())
            (root / "work/macromap-design").unlink()
            run = m.fresh_private_run(root, root / "work/macromap-design/fresh")
            target = root / "reports/keep.csv"
            target.write_text("unchanged")
            (run / "sample-map.csv").symlink_to(target)
            with self.assertRaises(FileExistsError):
                m.prepare(root, run)
            self.assertEqual(target.read_text(), "unchanged")
            with self.assertRaises(ValueError):
                m.fresh_private_run(root, root / "work/macromap-design/../../reports/elsewhere")

    def test_literal_parent_steps_cannot_hide_aliases_from_any_file_helper(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            private, public = directory / "private", directory / "reports"
            private.mkdir()
            (public / "nested").mkdir(parents=True)
            (private / "link").symlink_to(public / "nested", target_is_directory=True)
            unsafe = private / "link/../private.csv"
            source = private / "source"
            source.write_text("source")
            operations = [lambda: m.exclusive_text(unsafe, "private"),
                          lambda: m.write_csv(unsafe, [{"secret": "private"}]),
                          lambda: m.save_npz_exclusive(unsafe, {"data": np.array([1])}),
                          lambda: m.publish_without_replacement(source, unsafe),
                          lambda: m.sha256(unsafe), lambda: m.read_csv(unsafe)]
            for operation in operations:
                with self.assertRaises(ValueError):
                    operation()
                self.assertFalse((public / "private.csv").exists())
                self.assertTrue(source.exists())

    def test_source_change_during_normalization_does_not_publish_cache(self):
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source, output = directory / "toy.gz", directory / "normalized"
            CountCacheTests().write_counts(source, [[1, 2]])
            original_parser = m.exact_count_row
            def mutate_after_read(text, columns):
                values = original_parser(text, columns)
                # A timestamp/inode-signature change is enough to block this
                # cache, even if a concurrent source writer restores bytes.
                stat = source.stat()
                import os
                os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1000))
                return values
            with patch.object(m, "exact_count_row", side_effect=mutate_after_read):
                with self.assertRaises(ValueError):
                    m.normalize_count_cache(source, output, ["ENSG00000000000.1"], ["s1", "s2"])
            self.assertFalse(output.exists())
            self.assertFalse(output.with_name(output.name + ".partial").exists())
            self.assertFalse(output.with_name(output.name + ".lock").exists())

    def test_source_symlinks_are_rejected_even_if_bytes_match(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "source.gz"
            CountCacheTests().write_counts(source, [[1, 2]])
            alias = directory / "alias.gz"
            alias.symlink_to(source)
            with self.assertRaises(ValueError):
                m.sha256(alias)
            with self.assertRaises(ValueError):
                m.normalize_count_cache(alias, directory / "out", ["ENSG00000000000.1"], ["s1", "s2"])
            self.assertFalse((directory / "out").exists())

    def test_exclusive_writes_and_publication_reject_artifact_aliases(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            target = directory / "sentinel"
            target.write_text("unchanged")
            link = directory / "alias"
            link.symlink_to(target)
            with self.assertRaises(ValueError):
                m.write_csv(link, [{"x": 1}])
            with self.assertRaises(ValueError):
                m.exclusive_text(link, "private")
            source = directory / "source"
            source.write_text("new")
            with self.assertRaises(ValueError):
                m.publish_without_replacement(source, link)
            with self.assertRaises(FileExistsError):
                m.publish_without_replacement(source, target)
            self.assertEqual(target.read_text(), "unchanged")
            self.assertTrue(source.exists())

    def test_cache_namespace_reservation_blocks_concurrent_writer(self):
        from concurrent.futures import ThreadPoolExecutor
        from threading import Event
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source, output = directory / "toy.gz", directory / "normalized"
            CountCacheTests().write_counts(source, [[1, 2]])
            entered, release = Event(), Event()
            original_parser = m.exact_count_row
            def pause_after_reservation(text, columns):
                entered.set()
                if not release.wait(timeout=5):
                    raise RuntimeError("Test synchronization timed out")
                return original_parser(text, columns)
            with patch.object(m, "exact_count_row", side_effect=pause_after_reservation):
                with ThreadPoolExecutor(max_workers=1) as executor:
                    first = executor.submit(m.normalize_count_cache, source, output,
                                            ["ENSG00000000000.1"], ["s1", "s2"])
                    self.assertTrue(entered.wait(timeout=5))
                    try:
                        with self.assertRaises(FileExistsError):
                            m.normalize_count_cache(source, output, ["ENSG00000000000.1"], ["s1", "s2"])
                    finally:
                        release.set()
                    audit = first.result(timeout=5)
            self.assertEqual(m.sha256(output), audit["cache_sha256"])
            self.assertTrue(output.with_name(output.name + ".lock").is_file())
            with self.assertRaises(FileExistsError):
                m.normalize_count_cache(source, output, ["ENSG00000000000.1"], ["s1", "s2"])
            self.assertEqual(m.sha256(output), audit["cache_sha256"])

    def test_count_cache_dangling_names_and_parent_alias_are_never_followed(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "toy.gz"
            CountCacheTests().write_counts(source, [[1, 2]])
            for suffix in ("", ".partial", ".json", ".json.partial", ".lock"):
                output = directory / ("output" + str(len(suffix)))
                target = directory / ("must-not-create" + str(len(suffix)))
                alias = output.with_name(output.name + suffix)
                alias.symlink_to(target)
                with self.subTest(suffix=suffix), self.assertRaises(ValueError):
                    m.normalize_count_cache(source, output, ["ENSG00000000000.1"], ["s1", "s2"])
                self.assertFalse(target.exists())
                self.assertTrue(alias.is_symlink())
                alias.unlink()
            external = directory / "external"
            external.mkdir()
            parent_alias = directory / "parent-alias"
            parent_alias.symlink_to(external, target_is_directory=True)
            with self.assertRaises(ValueError):
                m.normalize_count_cache(source, parent_alias / "output", ["ENSG00000000000.1"], ["s1", "s2"])
            self.assertFalse((external / "output").exists())


class CountCacheTests(unittest.TestCase):
    def write_counts(self, path, rows):
        text = "\t".join(m.FEATURE_FIELDS + ("s1", "s2")) + "\n"
        for i, values in enumerate(rows):
            text += f"1\tENSG{i:011d}.1\t1\t10\t+\t10\t" + "\t".join(map(str, values)) + "\n"
        with gzip.open(path, "wt") as fh:
            fh.write(text)

    def test_released_all_features_denominator_measured_zeros_chunk_invariance(self):
        counts = np.array([[2, 4], [18, 36], [0, 0]])
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source = directory / "toy.gz"
            self.write_counts(source, counts)
            ids = [f"ENSG{i:011d}.1" for i in range(3)]
            results = []
            for chunk in (1, 3):
                output = directory / f"normalized-{chunk}"
                audit = m.normalize_count_cache(source, output, ids, ["s1", "s2"], chunk_rows=chunk)
                self.assertTrue(audit["complete"])
                array = np.memmap(output, dtype="<f8", mode="r", shape=(3, 2))
                results.append(np.array(array))
                del array
            assert_array_equal(results[0], results[1])
            expected = np.log2(counts / counts.sum(axis=0) * 1e6 + 1)
            assert_allclose(results[0], expected, atol=1e-14)
            assert_array_equal(results[0][-1], 0)

    def test_exact_count_tokens_reject_fractional_rounding_and_underflow(self):
        invalid = ("9.9999999999999999", "1.00000000000000001", "-1e-400", "1e-400",
                   "9007199254740991.00000000000001", "9007199254740992", "1_000", "NaN", "Inf")
        for token in invalid:
            with self.subTest(token=token), self.assertRaises(ValueError):
                m.exact_count_row(token + "\t1", 2)
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            for i, token in enumerate(invalid[:4]):
                source, output = directory / f"precision-{i}.gz", directory / f"precision-{i}.dat"
                self.write_counts(source, [[token, "1"]])
                with self.subTest(token=token), self.assertRaises(ValueError):
                    m.normalize_count_cache(source, output, ["ENSG00000000000.1"], ["s1", "s2"])
                self.assertFalse(output.exists())
                self.assertFalse(output.with_name(output.name + ".partial").exists())

    def test_exact_integral_decimal_forms_preserve_source_quantity(self):
        valid = ("0001", "+10.000", "10.", "1e1", "100e-1", ".0", "-0e-400", "9007199254740991")
        expected = [1, 10, 10, 10, 10, 0, 0, 9007199254740991]
        assert_array_equal(m.exact_count_row("\t".join(valid), len(valid)), expected)
        assert_array_equal(m.exact_count_row("0\t+10.00\t0002\t3.", 4), [0, 10, 2, 3])

    def test_fractional_negative_missing_and_zero_library_fail_no_output(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            for i, values in enumerate(([[0.5, 1]], [[-1, 1]], [["NA", 1]], [[0, 1]], [[1]], [[1, 2, 3]])):
                source, output = directory / f"bad-{i}.gz", directory / f"bad-{i}.dat"
                self.write_counts(source, values)
                with self.assertRaises(ValueError):
                    m.normalize_count_cache(source, output, ["ENSG00000000000.1"], ["s1", "s2"])
                self.assertFalse(output.exists())
                self.assertFalse(output.with_name(output.name + ".partial").exists())


if __name__ == "__main__":
    unittest.main()
