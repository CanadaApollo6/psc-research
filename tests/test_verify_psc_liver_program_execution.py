"""Artificial native-writer fixtures only. Never read real count/score payloads."""
import copy
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest import mock
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj
E = module("execution_verifier_test_module", ROOT / "scripts/verify_psc_liver_program_execution.py")
# Author is used ONLY to emit artificial native fixtures. The adapter never imports it.
A = module("artificial_native_writer", ROOT / "scripts/analyze_psc_liver_programs.py")
V = E.load_math(ROOT)

def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, sort_keys=True, indent=2, allow_nan=False) + "\n")

def pin(root, path):
    return {"path": str(path.relative_to(root)), "sha256": E.digest(path), "bytes": path.stat().st_size}

def fixture(root, full=False, missing=False, scale=1):
    """Synthetic values, fake identities, temporary root and injected native auth."""
    group_sizes = {"PSC": 17, "ASC": 17, "AIH": 30} if full else {"PSC": 3, "ASC": 3, "AIH": 3}
    remaining = dict(group_sizes)
    groups = []
    while any(remaining.values()):
        for group in ("AIH", "PSC", "ASC"):
            if remaining[group]:
                groups.append(group)
                remaining[group] -= 1
    counts = np.array([[1 + ((i * i + 7 * i * j + 13 * j) % 53) for j in range(len(groups))] for i in range(24)] +
                      [[0] * len(groups) for _ in range(10)], dtype=np.float64)
    counts *= scale
    if full:
        counts = np.vstack((counts, np.ones((41096 - len(counts), len(groups)), dtype=np.float64)))
    labels = ["ARTIFICIAL_GENE_" + str(i) for i in range(len(counts))]
    units = [{"unit_id": "ARTIFICIAL_UNIT_" + str(i), "group": group, "column_index": i,
              "geo_accession": "ARTIFICIAL_ACCESSION_" + str(i), "unit_type": "source_patient_sample", "verified_person_id": None}
             for i, group in enumerate(groups)]
    specs = []
    for tier in (E.CORE, E.STRICT):
        for number in range(10):
            indices = list(range(24, 34)) if number == 9 else [(number + 2 * k + (tier == E.STRICT)) % 24 for k in range(10)]
            available = not (missing and number == 8)
            specs.append({"interface": tier, "program_id": "artificial_program_" + str(number),
                          "score_kind": "unsigned_equal_gene", "coverage": {"score_status": "available" if available else "unavailable_coverage"},
                          "components": [{"name": "all", "indices": indices, "source_members": 10}]})
    external = {"synthetic_fixture": True, "no_real_source_identifiers": True}
    plan = A.draft_plan(root)
    plan.update(state="frozen", expected_compiled_programs_sha256=E.canonical(specs))
    plan["root_authorization"] = {"approver_role": "root", "real_execution_authorized": True,
                                  "reference": "ARTIFICIAL_TEST_ONLY_NOT_REAL_PERMISSION", "approved_at": "synthetic_fixture"}
    plan["method"]["bootstrap"]["replicates"] = 10000 if full else 17
    raw = root / "fixture/artificial-original.tsv.gz"
    cache = root / "fixture/artificial-cache.npy"
    raw.parent.mkdir(parents=True)
    with gzip.open(raw, "wt", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow([""] + [u["unit_id"] for u in units])
        for label, values in zip(labels, counts):
            writer.writerow([label] + [str(int(value)) for value in values])
    np.save(cache, counts)
    feature_path = root / "fixture/artificial-features.json"
    unit_path = root / "fixture/artificial-units.json"
    write_json(feature_path, [{"label_as_submitted": label, "source_row_index_0based": i} for i, label in enumerate(labels)])
    write_json(unit_path, [{"matrix_column_id_as_submitted": u["unit_id"], "matrix_column_index_0based": u["column_index"],
                           "condition_as_submitted": u["group"], "geo_accession": u["geo_accession"], "separate_published_patient_identifier": None} for u in units])
    plan["analysis_id"] = "ARTIFICIAL_NATIVE_SCHEMA_FIXTURE"
    plan["dataset"] = {"accession": "ARTIFICIAL_NOT_STUDY_DATA", "shape": list(counts.shape), "groups": group_sizes,
                       "raw_matrix_pin": pin(root, raw), "matrix_pin": pin(root, cache)}
    review_path = root / "fixture/artificial-accepted-review.json"
    write_json(review_path, {"status": "pass_synthetic_and_guard_review", "real_execution_reviewed": False,
                            "compiled_programs_sha256": E.canonical(specs), "implementation_sha256": {}})
    review_pin = pin(root, review_path)
    plan["independent_review"] = {"status": "accepted_synthetic_and_guard_review", "reference": "ARTIFICIAL_REVIEW", "artifact_pin": review_pin}
    plan_path = root / "config/psc-liver-program-freeze.json"
    write_json(plan_path, plan)
    plan_sha = E.digest(plan_path)
    run = root / E.RUNS_REL / "artificial-native-run"
    # This invokes the REAL native writer with purely generated inputs. It is not
    # a positive authorization test on the real repository or its RNA payloads.
    with mock.patch.object(A, "authorize_execution", return_value=(plan, specs, units, external, run, cache, {})), mock.patch.object(A, "SHAPE", counts.shape):
        A.execute(plan_path, plan_sha, "artificial_test_driver_only", run, root)
    contract = {"synthetic_fixture": True, "shape": list(counts.shape), "groups": group_sizes, "method": plan["method"],
        "contrasts": plan["contrast_order"], "fixed_plan_fields": {k: v for k, v in plan.items() if k not in
            {"state", "root_authorization", "independent_review", "expected_compiled_programs_sha256"}},
        "required_plan_fields": sorted(plan), "plan_paths": [str(plan_path.relative_to(root))],
        "compiled_sha256": E.canonical(specs), "source_metadata_pins": {str(p.relative_to(root)): pin(root, p) for p in (feature_path, unit_path)},
        "raw_pin": pin(root, raw), "cache_pin": pin(root, cache), "feature_audit_path": str(feature_path.relative_to(root)),
        "unit_source_path": str(unit_path.relative_to(root)), "external_audit_sha256": E.digest(run / "private-external-identifier-audit.json"),
        "accepted_review_pin": review_pin, "immutable_pins": {}, "runtime": E.runtime()}
    return {"root": root, "contract": contract, "plan": plan, "plan_path": plan_path, "plan_sha": plan_sha,
            "run": run, "completion": run / E.COMPLETE, "completion_sha": E.digest(run / E.COMPLETE), "counts": counts}

def audit(f, output="verification"):
    return E.run_authenticated_audit(f["root"], f["contract"], f["plan_path"], f["plan_sha"],
        f["completion"], f["completion_sha"], f["root"] / E.OUTPUT_REL / output, E.REAL_ACK, V)

def repin_artifacts(f):
    receipt = E.json_read(f["completion"])
    for row in receipt["artifacts"]:
        path = f["run"] / row["path"]
        row.update(bytes=path.stat().st_size, sha256=E.digest(path))
    write_json(f["completion"], receipt)
    f["completion_sha"] = E.digest(f["completion"])

class NativeFixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base = tempfile.TemporaryDirectory()
        cls.original = fixture(Path(cls.base.name))
    @classmethod
    def tearDownClass(cls):
        cls.base.cleanup()
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name) / "root"
        shutil.copytree(self.original["root"], self.root)
        self.f = copy.deepcopy(self.original)
        old = self.original["root"]
        for key in ("root", "plan_path", "run", "completion"):
            self.f[key] = self.root / self.original[key].relative_to(old)
    def tearDown(self):
        self.temp.cleanup()
    def mutate(self, name, operation, repin=True):
        path = self.f["run"] / name
        value = E.json_read(path)
        operation(value)
        write_json(path, value)
        if repin:
            repin_artifacts(self.f)
    def refused(self):
        with self.assertRaises(E.AuditError):
            audit(self.f)
        self.assertFalse((self.root / E.OUTPUT_REL / "verification").exists())
    def test_complete_native_writer_fixture_and_private_aggregate_receipt(self):
        result = audit(self.f)
        self.assertFalse(result["real_data_audited"])
        self.assertEqual(result["checks"]["endpoint_rows_checked"], 60)
        self.assertEqual(result["checks"]["score_values_compared"], 180)
        self.assertEqual(result["checks"]["private_omissions_checked"], 540)
        self.assertEqual(result["checks"]["bootstrap_differences_compared"], 1020)
        path = self.root / E.OUTPUT_REL / "verification/audit-receipt.json"
        self.assertEqual(path.parent.stat().st_mode & 0o777, 0o700)
        self.assertEqual(path.stat().st_mode & 0o777, 0o600)
        text = path.read_text()
        for token in ("ARTIFICIAL_UNIT_", "ARTIFICIAL_GENE_", "ARTIFICIAL_ACCESSION_", '"values":', '"p_raw":'):
            self.assertNotIn(token, text)
    def test_deterministic_receipt_replay(self):
        self.assertEqual(audit(self.f, "first"), audit(self.f, "second"))
    def test_no_ack_no_expression_or_output_access(self):
        with mock.patch.object(E, "parse_original_gzip", side_effect=AssertionError("payload")) as read:
            with self.assertRaises(E.AuditError):
                E.run_authenticated_audit(self.root, self.f["contract"], None, None, None, None, None, None, V)
            read.assert_not_called()
    def test_stale_plan_hash(self):
        self.f["plan_sha"] = "0" * 64
        self.refused()
    def test_stale_completion_hash(self):
        self.f["completion_sha"] = "0" * 64
        self.refused()
    def test_unapproved_plan_rejected(self):
        value = self.f["plan"]
        value["root_authorization"]["real_execution_authorized"] = False
        write_json(self.f["plan_path"], value)
        self.f["plan_sha"] = E.digest(self.f["plan_path"])
        self.refused()
    def test_numeric_boolean_approval_refused(self):
        value = self.f["plan"]
        value["root_authorization"]["real_execution_authorized"] = 1
        write_json(self.f["plan_path"], value)
        self.f["plan_sha"] = E.digest(self.f["plan_path"])
        self.refused()
    def test_unfrozen_plan_refused(self):
        value = self.f["plan"]
        value["state"] = "draft_not_frozen"
        write_json(self.f["plan_path"], value)
        self.f["plan_sha"] = E.digest(self.f["plan_path"])
        self.refused()
    def test_method_changes_refused(self):
        value = copy.deepcopy(self.f["plan"])
        value["method"]["bootstrap"]["seed"] += 1
        write_json(self.f["plan_path"], value)
        self.f["plan_sha"] = E.digest(self.f["plan_path"])
        self.refused()
    def test_incomplete_run_receipt(self):
        receipt = E.json_read(self.f["completion"])
        receipt["status"] = "failed_not_complete"
        write_json(self.f["completion"], receipt)
        self.f["completion_sha"] = E.digest(self.f["completion"])
        self.refused()
    def test_extra_output(self):
        (self.f["run"] / "EXTRA.json").write_text("{}")
        self.refused()
    def test_running_sentinel_refused(self):
        (self.f["run"] / "RUNNING.json").write_text("{}")
        self.refused()
    def test_failed_sentinel_refused(self):
        (self.f["run"] / "FAILED.json").write_text("{}")
        self.refused()
    def test_missing_output(self):
        (self.f["run"] / "private-scores.json").unlink()
        self.refused()
    def test_unpinned_artifact_mutation(self):
        self.mutate("private-scores.json", lambda x: x[0]["values"].__setitem__(0, 99.0), repin=False)
        self.refused()
    def test_hash_repinned_score_error_still_fails(self):
        self.mutate("private-scores.json", lambda x: x[0]["values"].__setitem__(0, 99.0))
        self.refused()
    def test_library_total_error(self):
        self.mutate("private-library-totals.json", lambda x: x.__setitem__(0, x[0] + 1))
        self.refused()
    def test_missing_status_row(self):
        self.mutate("aggregate-public-summary.json", lambda x: x["endpoint_rows"].pop())
        self.refused()
    def test_extra_public_patient_identifier(self):
        self.mutate("aggregate-public-summary.json", lambda x: x.__setitem__("unit_id", "ARTIFICIAL_PRIVATE_ID"))
        self.refused()
    def test_strict_discovery_q_refused(self):
        self.mutate("aggregate-public-summary.json", lambda x: x["endpoint_rows"][30].__setitem__("q_bh", 0.01))
        self.refused()
    def test_missing_p_not_observed_one(self):
        self.mutate("aggregate-public-summary.json", lambda x: x["endpoint_rows"][27].__setitem__("p_raw", 1.0))
        self.refused()
    def test_missing_public_q_not_placeholder_one(self):
        self.mutate("aggregate-public-summary.json", lambda x: x["endpoint_rows"][27].__setitem__("q_bh", 1.0))
        self.refused()
    def test_wrong_private_correction_placeholder(self):
        self.mutate("private-correction-vector.json", lambda x: x["p_vector_for_correction"].__setitem__(27, 0.0))
        self.refused()
    def test_family_subset_refused(self):
        self.mutate("private-correction-vector.json", lambda x: x["ordered_keys"].pop())
        self.refused()
    def test_wrong_public_omission_denominator(self):
        self.mutate("aggregate-public-summary.json", lambda x: x["endpoint_rows"][0]["omission"].__setitem__("relevant_unit_denominator", 9))
        self.refused()
    def test_outside_contrast_omission_not_recomputed_or_changed(self):
        def change(rows):
            row = next(r for r in rows if not r["in_contrast"])
            row["effect_delta"] = 1.0
        self.mutate("private-all-unit-omissions.json", change)
        self.refused()
    def test_full_omission_statistics_checked(self):
        self.mutate("private-all-unit-omissions.json", lambda x: x[0]["statistics"].__setitem__("degrees_of_freedom", 1.0))
        self.refused()
    def test_source_units_changed(self):
        self.mutate("private-source-units.json", lambda x: x[0].__setitem__("group", "ASC"))
        self.refused()
    def test_compiled_mapping_cannot_self_attest(self):
        self.mutate("private-compiled-programs.json", lambda x: x[0]["components"][0]["indices"].__setitem__(0, 23))
        self.refused()
    def test_original_frozen_plan_bytes_must_match(self):
        path = self.f["run"] / "private-frozen-plan-original.json"
        path.write_text(path.read_text() + " ")
        repin_artifacts(self.f)
        self.refused()
    def test_changed_external_identity_audit(self):
        self.mutate("private-external-identifier-audit.json", lambda x: x.__setitem__("changed", True))
        self.refused()
    def test_output_reuse(self):
        output = self.root / E.OUTPUT_REL / "verification"
        output.mkdir(parents=True)
        with self.assertRaises(E.AuditError):
            audit(self.f)
        self.assertEqual(list(output.iterdir()), [])
    def test_run_directory_privacy(self):
        self.f["run"].chmod(0o755)
        self.refused()
    def test_output_artifact_symlink(self):
        path = self.f["run"] / "private-scores.json"
        target = self.root / "fixture/scores-copy.json"
        path.rename(target)
        path.symlink_to(target)
        self.refused()
    def test_output_artifact_hardlink(self):
        path = self.f["run"] / "private-scores.json"
        (self.root / "fixture/hardlink.json").hardlink_to(path)
        self.refused()
    def test_plan_hardlink_before_body_read(self):
        self.f["plan_path"].unlink()
        self.f["plan_path"].hardlink_to(self.root / self.f["contract"]["raw_pin"]["path"])
        with mock.patch.object(E, "json_read", wraps=E.json_read) as reader:
            self.refused()
            self.assertNotIn(self.f["plan_path"], [call.args[0] for call in reader.call_args_list])
    def test_unexpected_npz_member(self):
        path = self.f["run"] / "private-bootstrap-draws.npz"
        with zipfile.ZipFile(path, "a") as stream:
            stream.writestr("unexpected.npy", b"invalid")
        repin_artifacts(self.f)
        self.refused()
    def test_npz_duplicate_member(self):
        path = self.f["run"] / "private-bootstrap-draws.npz"
        with zipfile.ZipFile(path, "a") as stream:
            stream.writestr("PSC.npy", stream.read("PSC.npy"))
        repin_artifacts(self.f)
        self.refused()
    def test_wrong_rng_index(self):
        path = self.f["run"] / "private-bootstrap-draws.npz"
        with np.load(path, allow_pickle=False) as stream:
            values = {key: stream[key] for key in stream.files}
        values["PSC"][0, 0] = (values["PSC"][0, 0] + 1) % 9
        np.savez_compressed(path, **values)
        repin_artifacts(self.f)
        self.refused()
    def test_bootstrap_difference_error(self):
        path = self.f["run"] / "private-bootstrap-differences.npz"
        with np.load(path, allow_pickle=False) as stream:
            values = {key: stream[key] for key in stream.files}
        values[next(iter(values))][0] += 1
        np.savez_compressed(path, **values)
        repin_artifacts(self.f)
        self.refused()
    def test_cache_count_error_even_when_re_pinned(self):
        cache = self.root / self.f["contract"]["cache_pin"]["path"]
        values = np.load(cache, allow_pickle=False)
        values[0, 0] += 1
        np.save(cache, values)
        self.f["contract"]["cache_pin"] = pin(self.root, cache)
        self.refused()

    def test_tiny_nonzero_outside_delta_refused_exactly(self):
        def change(rows):
            next(r for r in rows if not r["in_contrast"])["effect_delta"] = 1e-13
        self.mutate("private-all-unit-omissions.json", change)
        self.refused()
    def test_outside_baseline_stat_copy_is_exact(self):
        def change(rows):
            next(r for r in rows if not r["in_contrast"])["statistics"]["effect"] += 1e-13
        self.mutate("private-all-unit-omissions.json", change)
        self.refused()
    def test_available_probability_zero_refused(self):
        self.mutate("aggregate-public-summary.json", lambda x: x["endpoint_rows"][0].__setitem__("p_raw", 0.0))
        self.refused()
    def test_available_probability_above_one_refused(self):
        self.mutate("aggregate-public-summary.json", lambda x: x["endpoint_rows"][0].__setitem__("p_raw", 1.0 + 1e-11))
        self.refused()
    def test_finite_positive_variance_cannot_be_zero(self):
        self.mutate("aggregate-public-summary.json", lambda x: x["endpoint_rows"][0].__setitem__("variance_a", 0.0))
        self.refused()
    def test_missing_private_placeholder_is_exact_one(self):
        self.mutate("private-correction-vector.json", lambda x: x["p_vector_for_correction"].__setitem__(27, 1.0 - 1e-13))
        self.refused()
    def test_tiny_negative_unsigned_score_refused(self):
        self.mutate("private-scores.json", lambda x: x[9]["values"].__setitem__(0, -1e-13))
        self.refused()

    def test_fixed_confidence_perturbation_is_not_tolerated(self):
        value = copy.deepcopy(self.f["plan"])
        value["method"]["confidence"] += 1e-11
        write_json(self.f["plan_path"], value)
        self.f["plan_sha"] = E.digest(self.f["plan_path"])
        self.refused()
    def test_fractional_tiny_library_total_change_refused(self):
        self.mutate("private-library-totals.json", lambda x: x.__setitem__(0, x[0] + 1e-9))
        self.refused()
    def test_degenerate_interval_must_be_exactly_degenerate(self):
        self.mutate("aggregate-public-summary.json", lambda x: x["endpoint_rows"][27]["bootstrap"].__setitem__("ci_high", 1e-13))
        self.refused()

    def test_re_pinned_raw_schema_and_integrality_errors(self):
        raw = self.root / self.f["contract"]["raw_pin"]["path"]
        with gzip.open(raw, "rt") as stream:
            original = list(csv.reader(stream, delimiter="\t"))
        changes = [
            lambda rows: rows[0].__setitem__(1, "WRONG_ARTIFICIAL_HEADER"),
            lambda rows: rows[1].__setitem__(0, "WRONG_ARTIFICIAL_GENE"),
            lambda rows: rows[1].__setitem__(1, "1.5"),
            lambda rows: rows[1].__setitem__(1, "-1"),
            lambda rows: rows[1].__setitem__(1, str(2**53 + 1)),
            lambda rows: rows[1].pop(),
            lambda rows: rows.pop(),
            lambda rows: rows.append(rows[-1])]
        for change in changes:
            with self.subTest(case=changes.index(change)):
                rows = copy.deepcopy(original)
                change(rows)
                with gzip.open(raw, "wt", newline="") as stream:
                    csv.writer(stream, delimiter="\t", lineterminator="\n").writerows(rows)
                self.f["contract"]["raw_pin"] = pin(self.root, raw)
                self.refused()
    def test_re_pinned_gzip_crc_error(self):
        raw = self.root / self.f["contract"]["raw_pin"]["path"]
        data = bytearray(raw.read_bytes())
        data[-8] ^= 1
        raw.write_bytes(data)
        self.f["contract"]["raw_pin"] = pin(self.root, raw)
        self.refused()
    def test_object_cache_cannot_load_pickle(self):
        cache = self.root / self.f["contract"]["cache_pin"]["path"]
        np.save(cache, self.f["counts"].astype(object))
        self.f["contract"]["cache_pin"] = pin(self.root, cache)
        self.refused()
    def test_unsafe_public_output_refused_before_count_read(self):
        with mock.patch.object(E, "parse_original_gzip", side_effect=AssertionError("count_read")) as reader:
            with self.assertRaises(E.AuditError):
                E.run_authenticated_audit(self.root, self.f["contract"], self.f["plan_path"], self.f["plan_sha"],
                    self.f["completion"], self.f["completion_sha"], self.root / "reports/unsafe", E.REAL_ACK, V)
            reader.assert_not_called()
    def test_run_collection_root_cannot_be_the_run(self):
        collection = self.root / E.RUNS_REL
        run = self.f["run"]
        for path in run.iterdir():
            path.rename(collection / path.name)
        run.rmdir()
        collection.chmod(0o700)
        self.f["run"] = collection
        self.f["completion"] = collection / E.COMPLETE
        self.refused()

class StandaloneBoundaryTests(unittest.TestCase):
    def test_large_integer_library_total_one_count_error_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            f = fixture(Path(temp), scale=20000000)
            path = f["run"] / "private-library-totals.json"
            totals = E.json_read(path)
            original = totals[0]
            totals[0] += 1
            self.assertTrue(E.math.isclose(totals[0], original, rel_tol=1e-10, abs_tol=1e-12))
            write_json(path, totals)
            repin_artifacts(f)
            with self.assertRaises(E.AuditError):
                audit(f)
    def test_nonrepresentable_sd_becomes_explicit_unavailable_variance(self):
        small = E.expected_welch(V, [0.0] * 16 + [5e-324], [0.0] * 17)
        self.assertEqual(small["welch_status"], "unavailable_numeric_variance_underflow")
        self.assertEqual(small["zero_variance_groups"], ["b"])
        large = E.expected_welch(V, [-1.7e308, 1.7e308], [0.0, 0.0])
        self.assertEqual(large["welch_status"], "unavailable_numeric_variance")
        self.assertEqual(large["effect"], 0.0)
    def test_small_probabilities_require_relative_not_absolute_agreement(self):
        with self.assertRaises(E.AuditError):
            E.Compare().equal({"p_raw": 1e-13}, {"p_raw": 1e-20})
        with self.assertRaises(E.AuditError):
            E.Compare().equal({"variance_a": 1e-250}, {"variance_a": 1e-300})
    def test_source_contract_covers_every_frozen_producer_dependency(self):
        contract = E.production_contract(ROOT)
        pins = contract["source_metadata_pins"]
        self.assertTrue(set(A.CLOSED_VERSIONED_PINS).issubset(pins))
        self.assertTrue({p["path"] for p in A.EXTERNAL_RAW_SOURCE_PINS.values()}.issubset(pins))
        self.assertIn(A.EXTERNAL_SOURCE_REGISTRY[A.ECM_PROGRAM]["schema_handoff_pin"]["path"], pins)

    def test_decimal_integer_parser(self):
        self.assertEqual([E.count_integer(x) for x in ("0", "1.0", "1e2", str(2**53))], [0, 1, 100, 2**53])
    def test_bad_numeric_counts(self):
        for value in ("NaN", "Infinity", "-1", "0.1", str(2**53 + 1), "", "1_000", " 1 "):
            with self.subTest(value=value), self.assertRaises(E.AuditError):
                E.count_integer(value)
    def test_exact_constant_welch_and_one_zero_arm(self):
        self.assertEqual(E.expected_welch(V, [2.1] * 30, [2.1] * 17)["effect"], 0)
        self.assertIsNone(E.expected_welch(V, [2.1] * 30, [2.1] * 17)["p_raw"])
        self.assertEqual(E.expected_welch(V, [0, 0], [0, 2])["degrees_of_freedom"], 1)
    def test_variance_underflow_is_not_exact_zero(self):
        result = E.expected_welch(V, [0, 1e-200], [1e-200, 2e-200])
        self.assertEqual(result["welch_status"], "unavailable_numeric_variance_underflow")
        self.assertEqual(result["zero_variance_groups"], [])
        self.assertIsNotNone(result["effect"])
    def test_tail_underflow_is_unavailable_not_zero_p(self):
        result = E.expected_welch(V, [0, 1e-150] * 15, [1.0] * 17)
        self.assertEqual(result["welch_status"], "unavailable_numeric_tail_underflow")
        self.assertIsNone(result["p_raw"])
        self.assertIsNone(result["degrees_of_freedom"])
    def test_missing_score_inventory_is_explicit(self):
        with tempfile.TemporaryDirectory() as temp:
            result = audit(fixture(Path(temp), missing=True))
            self.assertEqual(result["checks"]["unavailable_score_vectors_checked"], 2)
            self.assertEqual(result["checks"]["unavailable_core_p_values_checked"], 6)
            self.assertEqual(result["checks"]["bootstrap_differences_compared"], 918)
    def test_full_approved_size_and_seed_fixture(self):
        with tempfile.TemporaryDirectory() as temp:
            result = audit(fixture(Path(temp), full=True))
            c = result["checks"]
            self.assertEqual(c["raw_source_cells_and_cache_values_compared"], 41096 * 64)
            self.assertEqual((c["source_columns_checked"], c["full_library_totals_compared"], c["score_values_compared"],
                              c["shared_rng_indices_compared"], c["bootstrap_differences_compared"], c["private_omissions_checked"]),
                             (64, 64, 1280, 640000, 600000, 3840))
    def test_no_author_import_or_synthetic_bundle_relabel(self):
        text = (ROOT / E.SELF_REL).read_text()
        self.assertNotIn("analyze_psc_liver_programs", text)
        self.assertNotIn("verify_numeric_bundle(", text)
    def test_default_metadata_only_no_arrays(self):
        with mock.patch.object(E.np, "load", side_effect=AssertionError("array read")), mock.patch.object(E, "parse_original_gzip", side_effect=AssertionError("raw read")), mock.patch("sys.stdout", new_callable=io.StringIO) as output:
            self.assertEqual(E.main([]), 0)
            self.assertFalse(json.loads(output.getvalue())["real_expression_accessed"])
    def test_real_only_hash_arguments_not_ignored(self):
        with mock.patch("sys.stdout", new_callable=io.StringIO) as output:
            self.assertEqual(E.main(["--expected-plan-sha256", "0" * 64]), 2)
            self.assertEqual(json.loads(output.getvalue())["reason"], "real_arguments_without_mode")
    def test_real_mode_needs_root_ack(self):
        with mock.patch("sys.stdout", new_callable=io.StringIO) as output:
            self.assertEqual(E.main(["--verify-real"]), 2)
            self.assertEqual(json.loads(output.getvalue())["reason"], "explicit_root_real_audit_ack_required")
    def test_errors_never_export_values_or_tracebacks(self):
        with mock.patch.object(E, "production_contract", side_effect=ValueError("PRIVATE_UNIT_AND_VALUE_12345")), mock.patch("sys.stdout", new_callable=io.StringIO) as output:
            self.assertEqual(E.main([]), 2)
            self.assertNotIn("PRIVATE_UNIT", output.getvalue())
            self.assertNotIn("12345", output.getvalue())
            self.assertNotIn("Traceback", output.getvalue())
    def test_duplicate_json_keys_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "value.json"
            path.write_text('{"key":1,"key":2}')
            with self.assertRaises(E.AuditError):
                E.json_read(path)
    def test_nan_json_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "value.json"
            path.write_text('{"key":NaN}')
            with self.assertRaises(E.AuditError):
                E.json_read(path)

if __name__ == "__main__":
    unittest.main()
