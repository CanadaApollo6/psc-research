"""Reporting-only tests. No native analysis imports, private payloads or reruns."""
from __future__ import annotations

from copy import deepcopy
import csv
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).absolute().parent.parent
SPEC = importlib.util.spec_from_file_location("psc_liver_reporter", ROOT / "scripts/report_psc_liver_programs.py")
reporter = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(reporter)


class ReportingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sha = reporter.digest((ROOT / reporter.SCRIPT).read_bytes())
        cls.context = reporter.load_context(ROOT, cls.sha)
        cls.aggregate = cls.context["aggregate"]
        cls.documents = {key: reporter.parse_json(reporter.authenticated(ROOT, pin))
                         for key, pin in reporter.INPUT_PINS.items()
                         if key not in ("protocol", "verification_design")}
        cls.bundle = reporter.build_bundle(cls.context)

    def changed(self, change):
        aggregate = deepcopy(self.aggregate)
        change(aggregate)
        return aggregate

    def refused(self, aggregate):
        with self.assertRaises(reporter.ReportingError):
            reporter.project_tables(aggregate)

    def sandbox(self):
        temporary = tempfile.TemporaryDirectory(prefix="test-", dir=ROOT / reporter.WORK)
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        paths = [pin["path"] for pin in reporter.INPUT_PINS.values()] + [reporter.SCRIPT, reporter.TESTS, reporter.DESIGN]
        for relative in paths:
            path = root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes((ROOT / relative).read_bytes())
        (root / reporter.WORK).mkdir(parents=True, exist_ok=True)
        (root / "data/derived").mkdir(parents=True, exist_ok=True)
        (root / "reports/figures").mkdir(parents=True, exist_ok=True)
        return root, reporter.load_context(root, self.sha)

    def saved_preview(self, root, context, name="preview"):
        relative = reporter.WORK + "/" + name
        with patch.object(reporter, "build_bundle", return_value=self.bundle):
            result = reporter.write_preview(context, relative)
        return relative, result["completion_sha256"]

    def test_exact_grid_and_snapshot_accepted(self):
        reporter.validate_aggregate(self.aggregate)
        self.assertEqual(len(self.aggregate["endpoint_rows"]), 60)
        expected = [(tier, program, contrast[0]) for tier in reporter.TIERS
                    for program in reporter.PROGRAMS for contrast in reporter.CONTRASTS]
        observed = [(r["interface"], r["program_id"], r["contrast_id"]) for r in self.aggregate["endpoint_rows"]]
        self.assertEqual(expected, observed)

    def test_missing_extra_duplicate_reordered_grid_refused(self):
        for kind in ("missing", "extra", "duplicate", "reordered"):
            with self.subTest(kind=kind):
                value = deepcopy(self.aggregate)
                rows = value["endpoint_rows"]
                if kind == "missing": rows.pop()
                if kind == "extra": rows.append(deepcopy(rows[0]))
                if kind == "duplicate": rows[1] = deepcopy(rows[0])
                if kind == "reordered": rows[0], rows[1] = rows[1], rows[0]
                self.refused(value)

    def test_source_orders_groups_and_n_are_exact(self):
        for field, value in (("interface", "Strict"), ("program_id", "selected_program"),
                             ("contrast_id", "AIH_minus_PSC"), ("group_a", "ASC"),
                             ("n_a", True), ("n_b", 29)):
            with self.subTest(field=field):
                self.refused(self.changed(lambda x: x["endpoint_rows"][0].__setitem__(field, value)))
        self.refused(self.changed(lambda x: x["endpoint_order"].reverse()))

    def test_core_family_and_strict_no_q(self):
        core_fields = {"family_size": 29, "family_id": "selected_subset", "discovery_eligible": False,
                       "inferential_role": "identity_sensitivity_only", "p_was_substituted": True,
                       "p_for_correction": 1.0, "q_bh": None, "q_by": None}
        strict_fields = {"family_size": 30, "family_id": "Core_full_panel_three_contrasts",
                         "discovery_eligible": True, "q_bh": 0.0, "q_by": 1.0,
                         "p_for_correction": 1.0}
        for index, changes in ((0, core_fields), (30, strict_fields)):
            for field, value in changes.items():
                with self.subTest(index=index, field=field):
                    self.refused(self.changed(lambda x: x["endpoint_rows"][index].__setitem__(field, value)))
        for field, value in (("core_family_size", 29), ("strict_discoveries_permitted", True),
                             ("unavailable_core_p_values", 1)):
            self.refused(self.changed(lambda x: x.__setitem__(field, value)))

    def test_nulls_empty_fields_and_values_preserved(self):
        tables = reporter.project_tables(self.aggregate)
        endpoints = list(csv.DictReader(io.StringIO(tables[reporter.DATA + "/endpoints.csv"].decode())))
        self.assertEqual(len(endpoints), 60)
        for index, row in enumerate(endpoints):
            self.assertNotIn("p_for_correction", row)
            self.assertEqual(row["p_raw"], str(self.aggregate["endpoint_rows"][index]["p_raw"]))
            self.assertEqual(row["effect"], str(self.aggregate["endpoint_rows"][index]["effect"]))
            self.assertEqual(row["numeric_inference_failure"], "")
            self.assertEqual(row["zero_variance_groups"], "[]")
            if index >= 30:
                self.assertEqual(row["q_bh"], "")
                self.assertEqual(row["q_by"], "")
                self.assertEqual(row["family_size"], "")
        self.assertEqual(endpoints[0]["q_by"], "1.0")
        self.assertEqual(endpoints[0]["parent_gate"], "")
        self.assertEqual(endpoints[0]["target_literal_symbol"], "")
        for value in (None, 0, False, "", "1", [0], {"private_id": "CANARY"}):
            self.refused(self.changed(lambda x: x["endpoint_rows"][0].__setitem__("p_raw", value)))

    def test_complete_export_row_counts_and_projection_without_calculation(self):
        tables = reporter.project_tables(self.aggregate)
        for name, expected in reporter.ROW_COUNTS.items():
            rows = list(csv.DictReader(io.StringIO(tables[reporter.DATA + "/" + name].decode())))
            self.assertEqual(len(rows), expected)
        groups = list(csv.DictReader(io.StringIO(tables[reporter.DATA + "/group-summaries.csv"].decode())))
        first = self.aggregate["endpoint_rows"][0]
        self.assertEqual(groups[0]["mean"], str(first["mean_a"]))
        self.assertEqual(groups[0]["sample_variance_ddof1"], str(first["variance_a"]))
        self.assertEqual([row["group"] for row in groups[:3]], ["PSC", "ASC", "AIH"])

    def test_repeated_group_summary_must_agree(self):
        self.refused(self.changed(lambda x: x["endpoint_rows"][1].__setitem__("mean_a", 1.0)))
        self.refused(self.changed(lambda x: x["endpoint_rows"][1].__setitem__("variance_a", .2)))

    def test_nonoverlap_and_ecm_gate_special_cases(self):
        for index, mapped in ((24, 695), (54, 673)):
            row = self.aggregate["endpoint_rows"][index]
            self.assertFalse(row["coverage"]["components"][0]["coverage_gate"])
            self.assertTrue(row["coverage"]["parent_gate"])
            self.assertEqual(row["coverage"]["mapped_unique_total"], mapped)
            reporter.validate_aggregate(self.aggregate)
        self.assertEqual(self.aggregate["endpoint_rows"][27]["coverage"]["components"][0]["original_source_members"], 321)
        self.refused(self.changed(lambda x: x["endpoint_rows"][27]["coverage"]["components"][0].__setitem__("original_source_members", 340)))
        self.refused(self.changed(lambda x: x["endpoint_rows"][24]["coverage"].__setitem__("parent_gate", False)))

    def test_each_omission_denominator_and_count(self):
        for index in range(60):
            expected = [47, 34, 47][index % 3]
            row = self.aggregate["endpoint_rows"][index]
            self.assertEqual(row["omission"]["relevant_unit_denominator"], expected)
            for key in ("in_contrast_rows", "defined_effect_rows", "relevant_unit_denominator"):
                with self.subTest(index=index, key=key):
                    self.refused(self.changed(lambda x: x["endpoint_rows"][index]["omission"].__setitem__(key, 64)))
        for key, value in (("strict_sign_reversal_rows", 48), ("strict_sign_reversal_rows", -1),
                           ("max_absolute_effect_change", -1), ("status", "unavailable"),
                           ("effect_min", 999.0)):
            self.refused(self.changed(lambda x: x["endpoint_rows"][0]["omission"].__setitem__(key, value)))

    def test_nested_private_canaries_at_every_object_class(self):
        paths = [(), ("method",), ("method", "bootstrap"), ("method", "coverage"),
                 ("method", "multiplicity"), ("method", "omissions"), ("software",),
                 ("group_sizes",), ("public_omission_denominators",), ("external_candidates", 0)]
        paths += [("endpoint_rows", index) + tail for index in (0, 30, 24, 27)
                  for tail in ((), ("coverage",), ("coverage", "components", 0),
                               ("bootstrap",), ("omission",))]
        paths.append(("endpoint_rows", 27, "target_gene_membership"))
        for path in paths:
            with self.subTest(path=path):
                value = deepcopy(self.aggregate)
                target = value
                for key in path:
                    target = target[key]
                target["private_source_unit"] = "PRIVATE_CANARY_ID"
                self.refused(value)

    def test_nested_canaries_in_allowlisted_fields(self):
        paths = [("endpoint_rows", 0, "coverage", "reasons"),
                 ("endpoint_rows", 0, "zero_variance_groups"),
                 ("interpretation_limits",), ("endpoint_order",),
                 ("method", "bootstrap", "group_order"), ("method", "multiplicity", "methods")]
        for path in paths:
            value = deepcopy(self.aggregate)
            target = value
            for key in path[:-1]: target = target[key]
            target[path[-1]] = ["PRIVATE_CANARY_ID"]
            self.refused(value)
        for key in ("welch_status", "effect_status", "interface", "score_kind", "numeric_inference_failure"):
            self.refused(self.changed(lambda x: x["endpoint_rows"][0].__setitem__(key, "PRIVATE_CANARY_ID")))
        self.refused(self.changed(lambda x: x["endpoint_rows"][0]["bootstrap"].__setitem__("ci_low", {"unit": "CANARY"})))
        self.refused(self.changed(lambda x: x["endpoint_rows"][27]["target_gene_membership"].__setitem__("literal_symbol", "CANARY")))
        self.refused(self.changed(lambda x: x["external_candidates"][0].__setitem__("original_source_members", 0)))

    def test_type_nonfinite_domain_and_nested_interval_rejection(self):
        for key in reporter.NUMERIC_FIELDS:
            for value in (float("nan"), float("inf"), -float("inf"), True, "PRIVATE_CANARY", [1], {}, None, 10**400):
                with self.subTest(key=key, type=type(value).__name__):
                    self.refused(self.changed(lambda x: x["endpoint_rows"][0].__setitem__(key, value)))
        for key, value in (("standard_error", 0), ("degrees_of_freedom", 0),
                           ("variance_a", -1), ("mean_a", -1), ("p_raw", 1.1), ("log_p_raw", .01),
                           ("ci_low", 100), ("q_bh", .01), ("q_by", 1.1)):
            self.refused(self.changed(lambda x: x["endpoint_rows"][0].__setitem__(key, value)))
        for key, value in (("ci_low", None), ("ci_high", float("inf")), ("ci_low", 1),
                           ("valid_replicates", 9999), ("requested_replicates", True),
                           ("p_value_computed", True), ("status", "degenerate")):
            self.refused(self.changed(lambda x: x["endpoint_rows"][0]["bootstrap"].__setitem__(key, value)))

    def test_json_duplicate_nonfinite_and_overflow_refused(self):
        for text in ('{"a":1,"a":2}', '{"a":NaN}', '{"a":Infinity}', '{"a":-Infinity}'):
            with self.assertRaises(reporter.ReportingError): reporter.parse_json(text.encode())
        value = reporter.parse_json(reporter.json_bytes(self.aggregate).replace(b'0.11249742883016767', b'1e999', 1))
        self.refused(value)

    def test_allowed_input_hash_mismatches_and_reporter_hash(self):
        root, _ = self.sandbox()
        for key, pin in reporter.INPUT_PINS.items():
            path = root / pin["path"]
            original = path.read_bytes()
            path.write_bytes(original + b" ")
            with self.subTest(key=key), self.assertRaises(reporter.ReportingError):
                reporter.load_context(root, self.sha)
            path.write_bytes(original)
        with self.assertRaises(reporter.ReportingError): reporter.load_context(ROOT, "0" * 64)

    def test_receipt_crossbindings_and_recursive_private_canaries(self):
        for key, field, bad in (("independent_verification", "execution_receipt_sha256", "0"*64),
                                ("mechanical_replay", "actual_hashes_checked", 22),
                                ("mechanical_replay", "primary_completion_sha256", "0"*64),
                                ("handoff", "actual_source_identity_token_hits_in_aggregate", 1),
                                ("plan", "contrast_order", [])):
            documents = deepcopy(self.documents)
            documents[key][field] = bad
            with self.assertRaises(reporter.ReportingError): reporter.validate_receipts(documents)
        for subtree in ("checks", "agreement", "runtime"):
            documents = deepcopy(self.documents)
            documents["independent_verification"][subtree]["private_id"] = "CANARY"
            with self.assertRaises(reporter.ReportingError): reporter.validate_receipts(documents)
            context = dict(self.context, audit=documents["independent_verification"])
            with self.assertRaises(reporter.ReportingError): reporter.verification_summary(context)
        documents = deepcopy(self.documents)
        documents["independent_verification"]["agreement"]["relative_only_fields"].append("CANARY")
        with self.assertRaises(reporter.ReportingError): reporter.validate_receipts(documents)

    def test_unsafe_paths_reuse_source_collision_and_symlinks(self):
        root, context = self.sandbox()
        for relative in (reporter.WORK + "/../escape", reporter.WORK + "//new", "./" + reporter.WORK + "/new",
                         reporter.WORK + "/new/", "/tmp/absolute", reporter.DATA, reporter.INPUT_PINS["aggregate"]["path"]):
            with self.subTest(relative=relative), self.assertRaises(reporter.ReportingError):
                reporter.fresh_preview_path(root, relative)
        path = root / reporter.WORK / "exists"
        path.mkdir()
        with self.assertRaises(reporter.ReportingError): reporter.fresh_preview_path(root, reporter.WORK + "/exists")
        link = root / reporter.WORK / "alias"
        link.symlink_to(path, target_is_directory=True)
        with self.assertRaises(reporter.ReportingError): reporter.fresh_preview_path(root, reporter.WORK + "/alias")
        ancestor = root / "aliasdir"
        ancestor.symlink_to(root / reporter.WORK, target_is_directory=True)
        with self.assertRaises(reporter.ReportingError): reporter.safe_path(root, "aliasdir/new")
        source = root / reporter.INPUT_PINS["aggregate"]["path"]
        source.unlink()
        source.symlink_to(ROOT / reporter.INPUT_PINS["aggregate"]["path"])
        with self.assertRaises(reporter.ReportingError): reporter.load_context(root, self.sha)

    def test_hardlinks_fifo_and_exclusive_write_refuse(self):
        root, _ = self.sandbox()
        original = root / "file"
        original.write_bytes(b"original")
        alias = root / "alias"
        os.link(original, alias)
        with self.assertRaises(reporter.ReportingError): reporter.read_regular(root, "file")
        with self.assertRaises(FileExistsError): reporter.exclusive_write(original, b"replacement")
        self.assertEqual(original.read_bytes(), b"original")
        pipe = root / "fifo"
        os.mkfifo(pipe)
        with self.assertRaises(reporter.ReportingError): reporter.read_regular(root, "fifo")

    def test_prevalidation_refusal_leaves_no_result(self):
        root, context = self.sandbox()
        context["aggregate"] = deepcopy(self.aggregate)
        context["aggregate"]["endpoint_rows"][0]["coverage"]["reasons"] = ["CANARY"]
        destination = reporter.WORK + "/refused"
        with self.assertRaises(reporter.ReportingError): reporter.write_preview(context, destination)
        self.assertFalse((root / destination).exists())
        self.assertEqual(list((root / reporter.WORK).glob(".stage-*")), [])

    def test_deterministic_bundle_and_figure_replay(self):
        second = reporter.build_bundle(self.context)
        self.assertEqual(self.bundle, second)
        self.assertEqual(set(second), set(reporter.ALL_FILES))
        svg = second[reporter.FIGURE + ".svg"]
        self.assertNotIn(b"<dc:date>", svg)
        self.assertTrue(second[reporter.FIGURE + ".png"].startswith(b"\x89PNG\r\n\x1a\n"))
        self.assertNotIn(b"p_for_correction", second[reporter.DATA + "/endpoints.csv"])

    def test_actual_preview_inventory_completion_and_replay(self):
        root, context = self.sandbox()
        preview, completion = self.saved_preview(root, context)
        saved = reporter.verify_bundle(root, preview, completion, context)
        self.assertEqual(saved, self.bundle)
        with patch.object(reporter, "build_bundle", return_value=self.bundle):
            result = reporter.write_preview(context, reporter.WORK + "/replay", saved)
        self.assertEqual(result["status"], "all11_reporting_files_byte_identical")
        self.assertEqual(result["completion_sha256"], completion)
        with self.assertRaises(reporter.ReportingError): reporter.write_preview(context, preview)

    def test_incomplete_extra_tampered_and_symlink_previews_refused(self):
        root, context = self.sandbox()
        preview, completion = self.saved_preview(root, context)
        directory = root / preview
        extra = directory / "EXTRA"
        extra.write_bytes(b"x")
        with self.assertRaises(reporter.ReportingError): reporter.verify_bundle(root, preview, completion, context)
        extra.unlink()
        target = directory / reporter.FIGURE / "wrong"  # Not created; only path arithmetic.
        png = directory / (reporter.FIGURE + ".png")
        original = png.read_bytes()
        png.write_bytes(original + b"x")
        with self.assertRaises(reporter.ReportingError): reporter.verify_bundle(root, preview, completion, context)
        png.write_bytes(original)
        png.unlink()
        with self.assertRaises(reporter.ReportingError): reporter.verify_bundle(root, preview, completion, context)
        png.symlink_to(ROOT / reporter.SCRIPT)
        with self.assertRaises(reporter.ReportingError): reporter.verify_bundle(root, preview, completion, context)

    def test_extra_manifest_canary_even_when_resealed_refused(self):
        root, context = self.sandbox()
        preview, _ = self.saved_preview(root, context)
        directory = root / preview
        manifest = reporter.parse_json((directory / reporter.MANIFEST).read_bytes())
        manifest["private_canary"] = "CANARY"
        (directory / reporter.MANIFEST).write_bytes(reporter.json_bytes(manifest))
        completion = reporter.parse_json((directory / reporter.COMPLETION).read_bytes())
        completion["manifest_sha256"] = reporter.digest((directory / reporter.MANIFEST).read_bytes())
        completion["inventory"] = reporter.inventory({name: (directory / name).read_bytes()
                                                       for name in reporter.ALL_FILES if name != reporter.COMPLETION})
        encoded = reporter.json_bytes(completion)
        (directory / reporter.COMPLETION).write_bytes(encoded)
        with self.assertRaises(reporter.ReportingError):
            reporter.verify_bundle(root, preview, reporter.digest(encoded), context)

    def test_failed_preview_write_is_cleaned(self):
        root, context = self.sandbox()
        destination = reporter.WORK + "/failure"
        with patch.object(reporter, "build_bundle", return_value=self.bundle), \
             patch.object(reporter.os, "fsync", side_effect=OSError("synthetic fsync failure")):
            with self.assertRaises(OSError): reporter.write_preview(context, destination)
        self.assertFalse((root / destination).exists())
        self.assertEqual(list((root / reporter.WORK).glob(".stage-*")), [])

    def test_sandbox_publication_exact_bytes_ack_and_reuse(self):
        root, context = self.sandbox()
        preview, completion = self.saved_preview(root, context)
        with self.assertRaises(reporter.ReportingError): reporter.publish(context, preview, completion, "wrong")
        self.assertFalse((root / reporter.DATA).exists())
        result = reporter.publish(context, preview, completion, reporter.PUBLISH_ACK)
        self.assertEqual(result["status"], "published_exact_authenticated_preview")
        for name, data in self.bundle.items(): self.assertEqual((root / name).read_bytes(), data)
        with self.assertRaises(reporter.ReportingError): reporter.publish(context, preview, completion, reporter.PUBLISH_ACK)

    def test_each_public_write_failure_removes_owned_partials(self):
        # Ten payload writes plus the completion temporary write. Figure failures
        # must not leave files outside the new data directory.
        for failure_index in range(11):
            with self.subTest(failure_index=failure_index):
                root, context = self.sandbox()
                calls = [0]
                def fail_at(_fd):
                    index = calls[0]
                    calls[0] += 1
                    if index == failure_index:
                        raise OSError("synthetic fsync failure")
                with patch.object(reporter, "verify_bundle", return_value=self.bundle), \
                     patch.object(reporter.os, "fsync", side_effect=fail_at):
                    with self.assertRaises(OSError):
                        reporter.publish(context, reporter.WORK + "/unused", "unused", reporter.PUBLISH_ACK)
                self.assertFalse((root / reporter.DATA).exists())
                self.assertFalse((root / (reporter.FIGURE + ".png")).exists())
                self.assertFalse((root / (reporter.FIGURE + ".svg")).exists())

    def test_existing_figure_never_overwritten(self):
        root, context = self.sandbox()
        figure = root / (reporter.FIGURE + ".png")
        figure.write_bytes(b"keep")
        with patch.object(reporter, "verify_bundle", return_value=self.bundle):
            with self.assertRaises(reporter.ReportingError):
                reporter.publish(context, reporter.WORK + "/unused", "unused", reporter.PUBLISH_ACK)
        self.assertEqual(figure.read_bytes(), b"keep")
        self.assertFalse((root / reporter.DATA).exists())

    def test_atomic_completion_never_overwrites(self):
        root, _ = self.sandbox()
        destination = root / "completion"
        destination.write_bytes(b"original")
        temporary = root / "temporary"
        with self.assertRaises(FileExistsError): reporter.install_completion(temporary, destination, b"new")
        self.assertEqual(destination.read_bytes(), b"original")
        self.assertFalse(temporary.exists())


    def test_postlink_cleanup_failure_never_leaves_success_marker(self):
        root, context = self.sandbox()
        original_unlink = Path.unlink
        failed = []
        def fail_once(path, *args, **kwargs):
            if path.name == ".completion.tmp" and not failed:
                failed.append(True)
                raise OSError("synthetic unlink failure")
            return original_unlink(path, *args, **kwargs)
        with patch.object(reporter, "verify_bundle", return_value=self.bundle), \
             patch.object(Path, "unlink", new=fail_once):
            with self.assertRaises(OSError):
                reporter.publish(context, reporter.WORK + "/unused", "unused", reporter.PUBLISH_ACK)
        self.assertTrue(failed)
        self.assertFalse((root / reporter.COMPLETION).exists())
        self.assertFalse((root / reporter.DATA).exists())
        self.assertFalse((root / (reporter.FIGURE + ".png")).exists())
        self.assertFalse((root / (reporter.FIGURE + ".svg")).exists())

    def test_initial_write_fstat_failure_closes_fd_and_removes_new_file(self):
        root, _ = self.sandbox()
        path = root / "newfile"
        original_open = os.open
        descriptors = []
        def record_open(*args, **kwargs):
            fd = original_open(*args, **kwargs)
            descriptors.append(fd)
            return fd
        with patch.object(reporter.os, "open", side_effect=record_open), \
             patch.object(reporter.os, "fstat", side_effect=OSError("synthetic fstat failure")):
            with self.assertRaises(OSError): reporter.exclusive_write(path, b"new")
        self.assertFalse(path.exists())
        for fd in descriptors:
            with self.assertRaises(OSError): os.fstat(fd)

    def test_read_ignores_atime_only_change(self):
        root, _ = self.sandbox()
        path = root / "file"
        path.write_bytes(b"content")
        original_fstat = os.fstat
        calls = []
        def changed_atime(fd):
            original = original_fstat(fd)
            class StatView:
                def __getattr__(self, name):
                    if name == "st_atime": return original.st_atime + len(calls)
                    return getattr(original, name)
            calls.append(fd)
            return StatView()
        with patch.object(reporter.os, "fstat", side_effect=changed_atime):
            self.assertEqual(reporter.read_regular(root, "file"), b"content")


if __name__ == "__main__":
    unittest.main()
