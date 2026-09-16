"""Synthetic frozen-runner and private-audit tests; no real source scores."""
from __future__ import annotations
import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

ROOT = Path(__file__).resolve().parents[1]
for directory in (ROOT, ROOT / "tests"):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
from scripts import analyze_macromap_program_context as m
from test_macromap_program_context import validated_fixture

class ExecutionRunnerTests(unittest.TestCase):
    def test_execution_requires_separate_authorized_hash_and_unchanged_contract(self):
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "draft.json"
            path.write_text(json.dumps({"analysis_id": "macromap-fixed-program-context-v1", "status": "DRAFT_NOT_AUTHORIZED"}))
            with self.assertRaises(ValueError):
                m.run_frozen_analysis(path, expected_plan_sha256="wrong", root=root)
            with self.assertRaises(ValueError):
                m.run_frozen_analysis(path, expected_plan_sha256=m.sha256(path), root=root)
            self.assertFalse((root / "work").exists())
            plan = {"analysis_id": "macromap-fixed-program-context-v1", "status": "frozen", "execute_real_expression": True,
                    "frozen_at_utc": "2026-09-16T00:00:00+00:00", "script_sha256": m.sha256(Path(m.__file__)),
                    "runner_script_sha256": m.sha256(Path(m.__file__).with_name("run_macromap_program_context.py")),
                    "runtime": m.runtime_versions(), "method_contract": m.execution_contract(),
                    "prepared_directory": "work/macromap-design/prepared", "prepared_summary_sha256": "mock-for-unit-test",
                    "output_directory": "work/macromap-design/execution"}
            design = {"summary": {"preparation_script_sha256": plan["script_sha256"]}}
            with patch.object(m, "load_prepared_design", return_value=design), patch.object(m, "verify_pins"):
                self.assertIs(m.validate_execution_plan(plan, root=root), design)
                changed = copy.deepcopy(plan)
                changed["method_contract"]["ridge_lambda"] = 1
                with self.assertRaises(ValueError):
                    m.validate_execution_plan(changed, root=root)
                changed = {**plan, "output_directory": "reports/private-data"}
                with self.assertRaises(ValueError):
                    m.validate_execution_plan(changed, root=root)

    def test_runner_checks_normalized_source_pin_before_any_matching_or_fit(self):
        import gzip
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / m.RAW / "MacroMap_raw_expression.txt.gz"
            source.parent.mkdir(parents=True)
            with gzip.open(source, "wt") as fh:
                fh.write("\t".join(m.FEATURE_FIELDS + ("s1", "s2")) + "\n1\tENSG00000000000.1\t1\t10\t+\t10\t1\t2\n")
            plan = {"output_directory": "work/macromap-design/frozen-toy", "prepared_summary_sha256": "toy",
                    "runner_script_sha256": m.sha256(Path(m.__file__).with_name("run_macromap_program_context.py"))}
            path = root / "plan.json"
            path.write_text(json.dumps(plan))
            design = {"feature_ids": ["ENSG00000000000.1"], "sample_ids": ["s1", "s2"], "samples": [{}, {}]}
            with patch.object(m, "validate_execution_plan", return_value=design), \
                 patch.object(m, "SOURCE_PINS", {m.RAW + "/MacroMap_raw_expression.txt.gz": {"sha256": "0" * 64}}), \
                 patch.object(m, "execute_fixed_context") as engine:
                with self.assertRaisesRegex(ValueError, "Normalized raw source differs"):
                    m.run_frozen_analysis(path, expected_plan_sha256=m.sha256(path), root=root)
            engine.assert_not_called()
            output = root / plan["output_directory"]
            self.assertTrue((output / "execution-failed.json").is_file())
            self.assertFalse((output / "execution-complete.json").exists())

    def test_blockwise_program_means_detection_and_pair_direction(self):
        x = np.array([[0., 2., 4.], [1., 3., 5.], [2., 4., 6.], [0., 0., 8.]])
        targets = np.array([[0.5, 0.5, 0, 0], [0, 0, 0.5, 0.5]])
        backgrounds = targets[::-1]
        expected = x.T @ targets.T
        for chunk in (1, 4):
            raw, background, detected = m.score_program_blockwise(x, [0, 1, 2], targets, backgrounds, chunk_rows=chunk)
            assert_allclose(raw, expected)
            assert_allclose(background, x.T @ backgrounds.T)
            assert_allclose(detected, (x > 0).T @ targets.T)
            delta = m._pair_matrix([{"treatment_index": 2, "control_index": 0}], {0: 0, 1: 1, 2: 2}, raw - background)
            assert_allclose(delta[0], (raw - background)[2] - (raw - background)[0])

    def test_private_fit_bundle_reconstructs_scaler_probabilities_and_losses(self):
        train, test = [], []
        for group, bucket in (("train", train), ("test", test)):
            for i in range(3):
                for j, stimulus in enumerate(m.STIMULI):
                    bucket.append({"line": f"{group}-{i}", "run": f"{group}-r{i}", "time": 6,
                                   "protocol": m.PROTOCOLS[0], "stimulus": stimulus, "status": "available",
                                   "treatment_id": f"{group}-{i}-{j}", "control_id": f"{group}-{i}-control",
                                   "treatment_index": i * 10 + j + 1, "control_index": i * 10})
        rng = np.random.default_rng(17)
        tx, ex = rng.normal(size=(27, 5)), rng.normal(size=(27, 5))
        tx[:, 0] = 2  # constant training feature, still audited
        result = m.heldout_increment(tx, ex, train, test)
        self.assertEqual(result["status"], "available")
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            m.save_fit_bundle(directory, result, tx, ex, train, test, predictor_program="ets2_g1_dn", fold_gate={"status": "available"})
            with np.load(directory / "fit-arrays.npz", allow_pickle=False) as audit:
                assert_array_equal(audit["train_x_raw"], tx)
                assert_array_equal(audit["train_x"], tx)
                assert_array_equal(audit["train_y"], audit["train_class_index"])
                assert_array_equal(audit["test_lines"], audit["test_line"])
                self.assertEqual(audit["test_pair_ids"][0], m.canonical_pair_key(test[0]))
                assert_allclose(audit["baseline_logp"], audit["baseline_log_probabilities"])
                assert_allclose(audit["test_x_standardized"], (ex - audit["scaler_mean"]) / audit["scaler_scale"])
                self.assertTrue(audit["scaler_constant"][0])
                beta = audit["baseline_coefficients"]
                logits = audit["test_x_standardized"][:, :4] @ beta[1:] + beta[0]
                centered = logits - logits.max(axis=1, keepdims=True)
                logp = centered - np.log(np.exp(centered).sum(axis=1, keepdims=True))
                assert_allclose(audit["baseline_log_probabilities"], logp, atol=1e-13)
                loss = -logp[np.arange(27), audit["test_class_index"]]
                assert_allclose(audit["baseline_loss"], loss)
                self.assertFalse(set(audit["train_line"]) & set(audit["test_line"]))
            with self.assertRaises(ValueError):
                m.save_npz_exclusive(directory / "bad.npz", {"object": np.array([{}], dtype=object)})

    def test_pooled_response_is_descriptive_without_weighted_median_interval(self):
        result = m.clustered_response_summary([1, 2, 10], ["a", "b", "c"], ["r1", "r2", "r3"], compute_interval=False)
        self.assertEqual(result["median"], 2)
        self.assertEqual(result["interval_status"], "descriptive_only_no_pooled_interval")
        self.assertIsNone(result["median_interval"])
        self.assertFalse(m.execution_contract()["pooled_response_median_interval"])

    def test_pooled_response_no_interval_and_failed_values_not_dropped(self):
        result = m.pooled_response_summary([1., 2., 10.], ["a", "b", "c"], ["1", "2", "3"],
                                          [m.PROTOCOLS[0], m.PROTOCOLS[0], m.PROTOCOLS[1]])
        self.assertAlmostEqual(result["mean"], 13 / 3)
        self.assertEqual(result["median"], 2)
        self.assertEqual(result["mean_interval_status"], "not_computed_prespecified")
        self.assertEqual(result["median_interval_status"], "not_computed_prespecified")
        self.assertIsNone(result["mean_interval"])
        self.assertIsNone(result["median_interval"])
        self.assertFalse(m.execution_contract()["pooled_response_mean_interval"])
        failed = m.pooled_response_summary([1., 2., np.nan], ["a", "b", "c"], ["1", "2", "3"],
                                          [m.PROTOCOLS[0], m.PROTOCOLS[0], m.PROTOCOLS[1]])
        self.assertEqual(failed["status"], "missing_response")
        self.assertEqual(failed["original_contrast_protocol_line_shares"], {m.PROTOCOLS[0]: 2 / 3, m.PROTOCOLS[1]: 1 / 3})

    def test_failed_matching_keeps_all_fit_manifest_rows_and_blocks_primary(self):
        from unittest.mock import patch
        samples = validated_fixture()
        indices = {p: list(range(i * 10, (i + 1) * 10)) for i, p in enumerate(m.original.ORIGINAL_PROGRAM_IDS)}
        indices["ets2_g1_dn_without_comparators"] = indices["ets2_g1_dn"]
        design = {"samples": samples, "sample_ids": [r["sample_id"] for r in samples],
                  "feature_ids": [f"ENSG{i:011d}.1" for i in range(100)], "program_indices": indices,
                  "pool_exclusion_indices": list(range(100)), "mapping_gates": {p: True for p in m.PROGRAMS},
                  "folds": m.build_folds(samples)}
        original_response = m.clustered_response_summary
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            with patch.object(m, "fit_multinomial_ridge") as fitter, \
                 patch.object(m, "clustered_response_summary", side_effect=lambda *a, **k: original_response(*a, **k, replicates=2)):
                result = m.execute_fixed_context(np.ones((100, len(samples))), design, output)
            fitter.assert_not_called()
            manifest = json.loads((output / "fit-manifest.json").read_text())
            self.assertEqual(len(manifest), 80)
            self.assertTrue(all(r["status"] == "missing_predictor" for r in manifest))
            for scope in result:
                for variant in result[scope]:
                    self.assertEqual(result[scope][variant]["status"], "incomplete_predictions")
                    self.assertIsNone(result[scope][variant]["primary_gain"])
            with np.load(output / manifest[0]["fit_npz"], allow_pickle=False) as arrays:
                self.assertTrue(np.isnan(arrays["test_x"]).all())
                self.assertNotIn("baseline_logp", arrays.files)
                self.assertEqual(len(arrays["test_pair_ids"]), 18)
            with np.load(output / manifest[0]["stratum_npz"], allow_pickle=False) as arrays:
                self.assertFalse(arrays["matching_available"].any())
                self.assertTrue(np.isnan(arrays["background_weights"]).all())

    def test_end_to_end_synthetic_engine_and_private_reference_bundle(self):
        from unittest.mock import patch
        # Reduced draw counts ONLY inside this synthetic test keep unittest fast.
        # The execution freeze contract and real runner always require1000/5000.
        samples = validated_fixture()
        indices = {p: list(range(i * 10, (i + 1) * 10)) for i, p in enumerate(m.original.ORIGINAL_PROGRAM_IDS)}
        indices["ets2_g1_dn_without_comparators"] = indices["ets2_g1_dn"]
        design = {"samples": samples, "sample_ids": [r["sample_id"] for r in samples],
                  "feature_ids": [f"ENSG{i:011d}.1" for i in range(1000)],
                  "program_indices": indices, "pool_exclusion_indices": list(range(80)),
                  "mapping_gates": {p: True for p in m.PROGRAMS}, "folds": m.build_folds(samples)}
        expression = np.random.default_rng(12).uniform(size=(1000, len(samples)))
        matched, response, bootstrap = m.matched_weights, m.clustered_response_summary, m.cluster_bootstrap_gain
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            with patch.object(m, "matched_weights", side_effect=lambda *a, **k: matched(*a, **k, replicates=5)), \
                 patch.object(m, "clustered_response_summary", side_effect=lambda *a, **k: response(*a, **k, replicates=5)), \
                 patch.object(m, "cluster_bootstrap_gain", side_effect=lambda *a, **k: bootstrap(*a, **k, replicates=5)):
                result = m.execute_fixed_context(expression, design, output)
            manifest = json.loads((output / "fit-manifest.json").read_text())
            self.assertEqual(len(manifest), 80)
            self.assertTrue(all(row["status"] == "available" for row in manifest))
            self.assertTrue(all(m.sha256(output / row["fit_npz"]) == row["fit_npz_sha256"] for row in manifest))
            self.assertEqual(set(result), {"mapped_primary", "all_source_lines_sensitivity"})
            for scope in result:
                for variant in result[scope]:
                    self.assertEqual(result[scope][variant]["status"], "available")
                    self.assertEqual(result[scope][variant]["line_counts"], {p: 10 for p in m.PROTOCOLS})
            context = output / "mapped_primary/Mod_SmartSeq2-t6-fold0"
            with np.load(context / "score-reference-arrays.npz", allow_pickle=False) as audit:
                refs = audit["training_control_indices"]
                assert_allclose(audit["reference"], expression[:, refs].mean(axis=1))
                for index in refs:
                    sample = samples[int(index)]
                    self.assertEqual(sample["condition"], "Ctrl_6")
                    self.assertNotEqual(design["folds"][(sample["protocol"], sample["run"])], 0)
                assert_allclose(audit["sample_raw_scores"], expression[:, audit["selected_sample_indices"]].T @ audit["target_weights"].T, atol=1e-13)
                assert_allclose(audit["background_weights"][:, :80], 0)
            panel = json.loads((output / "mapped_primary/response-panel.json").read_text())
            self.assertTrue(all(r["status"] == "unavailable_mock_identity" for r in panel if r["stimulus"] == "PIC"))
            self.assertTrue(all(r["median_interval_status"] == "not_computed_prespecified" for r in panel
                                if r["protocol"] == "pooled_descriptive" and r["stimulus"] != "PIC"))


if __name__ == "__main__":
    unittest.main()
