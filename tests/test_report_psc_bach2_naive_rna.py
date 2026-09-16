"""Public aggregate only: schema, disclosure, deterministic presentation and file gates."""
import copy
import csv
import importlib.util
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location(
    "bach2_reporter", Path(__file__).resolve().parents[1] / "scripts/report_psc_bach2_naive_rna.py")
r = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(r)


class PublicValidationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.handoff, cls.blobs, cls.objects = r.load_sources()
        cls.schema = cls.objects["public_schema"]
        cls.expected = cls.objects["public_aggregate"]["provenance"]

    def setUp(self):
        self.public = copy.deepcopy(self.objects["public_aggregate"])

    def validate(self):
        r.validate_public(self.public, self.schema, self.expected)

    def test_actual_approved_schema_and_pins(self):
        self.validate()
        self.assertEqual(r.digest(self.blobs["public_aggregate"]), r.AGGREGATE_SHA)

    def test_recursive_private_canaries_rejected_at_every_object(self):
        def objects(value, path=()):
            if isinstance(value, dict):
                yield path
                for key, nested in value.items():
                    yield from objects(nested, path + (key,))
        for path in objects(self.public):
            with self.subTest(path=path):
                changed = copy.deepcopy(self.public)
                target = changed
                for key in path:
                    target = target[key]
                target["donor_private_canary"] = {"measurements": ["PRIVATE_CANARY"]}
                with self.assertRaises(r.ReportingError):
                    r.validate_public(changed, self.schema, self.expected)

    def test_decimal_strings_reject_nonfinite_nonnumeric_wrong_types(self):
        bad = ("NaN", "Infinity", "-Infinity", "sNaN", "", "0x1", " 1", "1 ",
               "1_0", "PRIVATE_CANARY", "1e999999", None, True, 1, 0.1, [], {})
        for value in bad:
            with self.subTest(value=value):
                self.public["groups"]["SNP"]["mean"] = value
                with self.assertRaises(r.ReportingError):
                    self.validate()

    def test_every_decimal_leaf_rejects_nested_private_payload(self):
        def leaves(schema, path=()):
            if isinstance(schema, dict):
                for key, nested in schema.items():
                    yield from leaves(nested, path + (key,))
            elif "decimal" in schema:
                yield path
        for path in leaves(self.schema):
            with self.subTest(path=path):
                changed = copy.deepcopy(self.public)
                target = changed
                for key in path[:-1]:
                    target = target[key]
                target[path[-1]] = {"donor": "PRIVATE_CANARY"}
                with self.assertRaises(r.ReportingError):
                    r.validate_public(changed, self.schema, self.expected)

    def test_exact_grid_and_integer_types(self):
        for group in ("SNP", "noSNP"):
            for n in (3, 5, 4.0, True, "4"):
                with self.subTest(group=group, n=n):
                    changed = copy.deepcopy(self.public)
                    changed["groups"][group]["n"] = n
                    with self.assertRaises(r.ReportingError):
                        r.validate_public(changed, self.schema, self.expected)
        self.public["groups"]["other"] = self.public["groups"].pop("noSNP")
        with self.assertRaises(r.ReportingError):
            self.validate()

    def test_ci_shape_order_missingness_and_finite_values(self):
        for ci in ([], ["0"], ["0", "1", "2"], ["1", "0"], ["0", "1"],
                   ["-1", "NaN"], ["-1", 1], None, {"donor": "PRIVATE_CANARY"}):
            with self.subTest(ci=ci):
                self.public["contrast"]["CI95"] = ci
                with self.assertRaises(r.ReportingError):
                    self.validate()

    def test_all_unseen_statuses_fail_closed(self):
        for status in ("unavailable_zero_SE", "unavailable_nonfinite_inference",
                       "CI_available_P_unavailable_numeric_tail", "passed", None):
            with self.subTest(status=status):
                self.public["contrast"]["inference_status"] = status
                with self.assertRaises(r.ReportingError):
                    self.validate()

    def test_inference_ranges_and_missing_fields(self):
        for key, values in {"SE": ["0", "-1"], "df": ["2", "7", None],
                            "P_two_sided": ["0", "-0.1", "1.1", None]}.items():
            for value in values:
                with self.subTest(key=key, value=value):
                    changed = copy.deepcopy(self.public)
                    changed["contrast"][key] = value
                    with self.assertRaises(r.ReportingError):
                        r.validate_public(changed, self.schema, self.expected)

    def test_loo_complete_denominator_counts_extrema_signs(self):
        for key, value in (("count", 7), ("count", "8"), ("negative_count", 7),
                           ("zero_count", -1), ("strict_sign_reversal_count", 1),
                           ("full_point_sign", 1), ("maximum", "0"), ("minimum", "1")):
            with self.subTest(key=key, value=value):
                changed = copy.deepcopy(self.public)
                changed["leave_one_donor_out"][key] = value
                with self.assertRaises(r.ReportingError):
                    r.validate_public(changed, self.schema, self.expected)

    def test_qc_fixed_gate_and_bounds(self):
        for key, value in (("RNA_rows", 36639), ("selected_cells", 0), ("donors", True),
                           ("source_and_membership_gates", {"donor": "passed"}),
                           ("BACH2_positive_selected_cells", -1),
                           ("BACH2_positive_selected_cells", 16660),
                           ("zero_BACH2_donors_retained", 9)):
            with self.subTest(key=key, value=value):
                changed = copy.deepcopy(self.public)
                changed["QC"][key] = value
                with self.assertRaises(r.ReportingError):
                    r.validate_public(changed, self.schema, self.expected)

    def test_provenance_drift_each_leaf(self):
        for key in self.expected:
            with self.subTest(key=key):
                changed = copy.deepcopy(self.public)
                changed["provenance"][key] = "0" * 64
                with self.assertRaises(r.ReportingError):
                    r.validate_public(changed, self.schema, self.expected)

    def test_duplicate_keys_and_json_constants(self):
        for blob in (b'{"x":1,"x":2}', b'{"x":{"y":1,"y":2}}', b'{"x":NaN}', b'{"x":Infinity}'):
            with self.subTest(blob=blob), self.assertRaises(r.ReportingError):
                r.parse_json(blob)

    def test_no_unknown_descriptor_or_arbitrary_passed_subtrees(self):
        with self.assertRaises(r.ReportingError):
            r.validate_tree({"private": "passed"}, "literal:passed")
        with self.assertRaises(r.ReportingError):
            r.validate_tree("passed", "anything")


class PresentationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payloads = r.prepare_payloads()
        cls.handoff, cls.blobs, cls.objects = r.load_sources()
        cls.real_root = r.ROOT
        (r.ROOT / r.WORK).mkdir(parents=True, exist_ok=True)

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="test-", dir=self.real_root / r.WORK)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.root_patch = patch.object(r, "ROOT", self.root)
        self.root_patch.start()
        self.addCleanup(self.root_patch.stop)
        (self.root / ".gitignore").write_text("work/\n")
        self.rel = r.WORK + "/fixture"

    def create(self):
        with patch.object(r, "prepare_payloads", return_value=self.payloads):
            return r.preview(self.rel)

    def test_all_exported_strings_and_rows_exact(self):
        public = self.objects["public_aggregate"]
        self.assertEqual(self.payloads[r.DATA + "/public-aggregate.json"], self.blobs["public_aggregate"])
        endpoint = list(csv.DictReader(io.StringIO(self.payloads[r.DATA + "/endpoint.csv"].decode())))
        self.assertEqual(len(endpoint), 1)
        self.assertEqual(tuple(endpoint[0]), r.ENDPOINT_COLUMNS)
        for key, value in public["contrast"].items():
            if key != "CI95":
                self.assertEqual(endpoint[0][key], value)
        self.assertEqual([endpoint[0][k] for k in ("CI95_lower", "CI95_upper")], public["contrast"]["CI95"])
        groups = list(csv.DictReader(io.StringIO(self.payloads[r.DATA + "/groups.csv"].decode())))
        self.assertEqual([g["group"] for g in groups], ["SNP", "noSNP"])
        for row in groups:
            self.assertEqual(tuple(row), r.GROUP_COLUMNS)
            for key, value in public["groups"][row["group"]].items():
                self.assertEqual(row[key], str(value))
        for row in endpoint + groups:
            for key in r.COMMON:
                self.assertEqual(row[key], str(public[key]))
            for key in r.PROVENANCE:
                self.assertEqual(row["provenance." + key], public["provenance"][key])
        influence = r.parse_json(self.payloads[r.DATA + "/influence-qc.json"])
        self.assertEqual(influence, {k: public[k] for k in r.COMMON + ("QC", "leave_one_donor_out", "provenance")})
        self.assertEqual(self.payloads[r.DATA + "/aggregate-verification.json"], self.blobs["independent_verification"])
        self.assertEqual(self.payloads[r.DATA + "/mechanical-replay-verification.json"], self.blobs["replay_verification"])

    def test_exact_full_presentation_replay_and_stable_figure_metadata(self):
        with patch.object(r, "ROOT", self.real_root):
            self.assertEqual(r.prepare_payloads(), self.payloads)
        svg = self.payloads[r.FIGURE + ".svg"]
        self.assertNotIn(b"<dc:date>", svg)
        self.assertNotIn(b"<text", svg)
        self.assertTrue(self.payloads[r.FIGURE + ".png"].startswith(b"\x89PNG\r\n\x1a\n"))

    def test_reads_only_seven_fixed_public_and_metadata_inputs(self):
        read = r.read_file
        visited = []
        def recording(relative, maximum=1_000_000):
            visited.append(relative)
            return read(relative, maximum)
        with patch.object(r, "ROOT", self.real_root), patch.object(r, "read_file", side_effect=recording):
            r.load_sources()
        self.assertEqual(set(visited), {r.HANDOFF_PATH} | {self.handoff[k]["path"] for k in r.SOURCE_KEYS})
        self.assertEqual(len(visited), 7)

    def test_actual_hash_drift_rejects_even_numerically_valid_p1_substitution(self):
        original = r.read_file
        source_path = self.handoff["public_aggregate"]["path"]
        changed = copy.deepcopy(self.objects["public_aggregate"])
        changed["contrast"]["P_two_sided"] = "1.00000000000000000"
        altered = r.json_bytes(changed)
        def replacement(relative, maximum=1_000_000):
            return altered if relative == source_path else original(relative, maximum)
        with patch.object(r, "ROOT", self.real_root), patch.object(r, "read_file", side_effect=replacement):
            with self.assertRaises(r.ReportingError):
                r.load_sources()

    def test_fresh_bundle_complete_then_exact_verify(self):
        result = self.create()
        with patch.object(r, "prepare_payloads", return_value=self.payloads):
            self.assertEqual(r.verify(self.rel, result["completion_sha256"])["status"], "EXACT_PRESENTATION_REPLAY_PASS")
        self.assertEqual(result["files"], len(r.INVENTORY))

    def test_existing_directory_alias_traversal_and_source_collision_rejected(self):
        self.create()
        bad = (self.rel, self.rel + "/../other", self.rel + "/", "./" + self.rel,
               "/" + self.rel, self.rel.replace("/fixture", "//fixture"),
               "work/psc-bach2-naive-rna-plan/runs/primary-v1", r.DATA)
        for relative in bad:
            with self.subTest(relative=relative), self.assertRaises(r.ReportingError):
                r.preview(relative)

    def test_symlink_and_hardlink_inputs_rejected(self):
        original = self.root / "original.json"
        original.write_text("{}")
        (self.root / "alias.json").symlink_to(original)
        with self.assertRaises(r.ReportingError):
            r.read_file("alias.json")
        os.link(original, self.root / "hard.json")
        for name in ("hard.json", "original.json"):
            with self.assertRaises(r.ReportingError):
                r.read_file(name)

    def test_symlink_output_ancestor_rejected(self):
        (self.root / "work").symlink_to(self.real_root / r.WORK, target_is_directory=True)
        with self.assertRaises(r.ReportingError):
            r.preview(self.rel)

    def test_failure_before_or_during_completion_never_leaves_complete_bundle(self):
        for stage in ("write_new", "write_completion"):
            with self.subTest(stage=stage):
                with patch.object(r, "prepare_payloads", return_value=self.payloads), \
                     patch.object(r, stage, side_effect=OSError("synthetic write failure")):
                    with self.assertRaises(OSError):
                        r.preview(self.rel)
                self.assertFalse((self.root / self.rel).exists())

    def test_invalid_input_writes_no_preview(self):
        with patch.object(r, "prepare_payloads", side_effect=r.ReportingError("invalid public input")):
            with self.assertRaises(r.ReportingError):
                r.preview(self.rel)
        self.assertFalse((self.root / self.rel).exists())

    def test_wrong_completion_hash_rejected(self):
        self.create()
        with patch.object(r, "prepare_payloads", return_value=self.payloads):
            with self.assertRaises(r.ReportingError):
                r.verify(self.rel, "0" * 64)

    def test_extra_file_directory_symlink_and_missing_file_rejected(self):
        result = self.create()
        extra = self.root / self.rel / "PRIVATE_CANARY.json"
        extra.write_text("{}")
        with self.assertRaises(r.ReportingError):
            r.verify(self.rel, result["completion_sha256"])
        extra.unlink()
        extra.mkdir()
        with self.assertRaises(r.ReportingError):
            r.verify(self.rel, result["completion_sha256"])
        extra.rmdir()
        extra.symlink_to(self.root / self.rel / r.COMPLETION)
        with self.assertRaises(r.ReportingError):
            r.verify(self.rel, result["completion_sha256"])
        extra.unlink()
        (self.root / self.rel / r.DATA / "groups.csv").unlink()
        with self.assertRaises(r.ReportingError):
            r.verify(self.rel, result["completion_sha256"])

    def test_nested_export_canary_and_payload_tampering_rejected(self):
        result = self.create()
        output = self.root / self.rel / r.DATA / "influence-qc.json"
        changed = r.parse_json(output.read_bytes())
        changed["QC"]["private"] = "PRIVATE_CANARY"
        output.write_bytes(r.json_bytes(changed))
        with patch.object(r, "prepare_payloads", return_value=self.payloads):
            with self.assertRaises(r.ReportingError):
                r.verify(self.rel, result["completion_sha256"])

    def test_publish_cannot_start_without_exact_separate_root_approval(self):
        # No publication is executed: missing and wrong approvals fail before destination creation.
        result = self.create()
        approval = "work/psc-bach2-naive-rna-freeze/test-root-approval.json"
        with patch.object(r, "verified_payloads", return_value=self.payloads):
            with self.assertRaises(r.ReportingError):
                r.publish(self.rel, result["completion_sha256"], approval, "0" * 64)
            path = self.root / approval
            path.parent.mkdir(parents=True)
            body = r.approval_template(self.rel, result["completion_sha256"])
            body["schema_version"] = True  # bool must not masquerade as integer one.
            path.write_bytes(r.json_bytes(body))
            with self.assertRaises(r.ReportingError):
                r.publish(self.rel, result["completion_sha256"], approval, r.digest(path.read_bytes()))
        self.assertFalse((self.root / r.DATA).exists())


if __name__ == "__main__":
    unittest.main()
