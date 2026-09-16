"""Independent annotation/selection and synthetic endpoint verification only."""
import csv
from fractions import Fraction
import importlib.util
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("independent_bach2_naive", ROOT / "scripts/verify_psc_bach2_naive_rna.py")
V = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)


class IndependentMath(unittest.TestCase):
    def test_exact_constant_not_roundoff_variance(self):
        mean, var = V.exact_moments([0.1] * 4)
        self.assertEqual(mean, Fraction(0.1))
        self.assertEqual(var, 0)

    def test_genuine_tiny_variation_not_clamped(self):
        a = [0.1, 0.1, 0.1, math.nextafter(0.1, math.inf)]
        self.assertGreater(V.exact_moments(a)[1], 0)
        out = V.independent_welch(a, [0.1] * 4)
        self.assertEqual(out["inference_status"], "available")
        self.assertEqual(out["df"], 3)

    def test_sub_double_variance_has_finite_scaled_se(self):
        out = V.independent_welch([0.0, 1e-200, 2e-200, 3e-200], [1e-200] * 4)
        self.assertEqual(out["inference_status"], "available")
        self.assertGreater(out["se"], 0)
        self.assertEqual(out["df"], 3)

    def test_both_constant_point_only_no_forced_p(self):
        for a, b in [(0.1, 0.1), (0.1, 0.2), (0.0, 1.0)]:
            out = V.independent_welch([a] * 4, [b] * 4)
            self.assertEqual(out["delta"], float(Fraction(a) - Fraction(b)))
            self.assertEqual(out["inference_status"], "point_only")
            self.assertIsNone(out["ci95"])
            self.assertIsNone(out["p_two_sided"])

    def test_ordinary_welch_actual_df(self):
        out = V.independent_welch([0, 1, 2, 3], [1, 2, 3, 4])
        self.assertEqual(out["delta"], -1.0)
        self.assertAlmostEqual(out["se"], math.sqrt(5/6), places=15)
        self.assertEqual(out["df"], 6.0)
        critical = (out["ci95"][1] - out["delta"]) / out["se"]
        self.assertAlmostEqual(critical, 2.44691185114497, places=11)
        self.assertAlmostEqual(out["p_two_sided"], 0.3153335962012298, places=12)

    def test_one_constant_arm_valid(self):
        out = V.independent_welch([0.1] * 4, [0, 1, 2, 5])
        self.assertEqual(out["inference_status"], "available")
        self.assertEqual(out["df"], 3.0)
        self.assertGreater(out["p_two_sided"], 0)

    def test_group_swap_sign_and_same_inference(self):
        a, b = [0, 1, 2, 5], [0.2, 0.4, 0.6, 0.8]
        one, two = V.independent_welch(a,b), V.independent_welch(b,a)
        self.assertEqual(one["delta"], -two["delta"])
        self.assertEqual(one["p_two_sided"], two["p_two_sided"])
        self.assertEqual(one["df"], two["df"])
        self.assertEqual(one["ci95"], [-two["ci95"][1], -two["ci95"][0]])

    def test_missing_nonfinite_donor_rejected(self):
        for a in ([0, 1, 2], [0, 1, 2, math.nan], [0, 1, 2, math.inf]):
            with self.assertRaises(ValueError):
                V.independent_welch(a, [0] * 4)

    def test_source_units_to_y_zero_and_fixed_pseudocount(self):
        self.assertEqual(V.independent_y(0, 9), 0)
        self.assertEqual(V.independent_y(1, 1000000), 1)
        self.assertEqual(V.independent_y(3, 1000000), 2)
        self.assertAlmostEqual(V.independent_y(1, 1), math.log2(1000001), places=14)

    def test_invalid_full_rna_denominator_or_totals(self):
        for args in [(0, 0), (1, -1), (-1, 10), (11, 10), (1.0, 10), (True, 10)]:
            with self.assertRaises(ValueError):
                V.independent_y(*args)

    def test_eight_loo_points_retain_both_groups(self):
        out = V.independent_loo([0, 1, 2, 3], [4, 5, 6, 7])
        self.assertEqual(len(out), 8)
        self.assertEqual(out, [-3.5, -23/6, -25/6, -4.5, -4.5, -25/6, -23/6, -3.5])
        self.assertTrue(all(isinstance(value,float) for value in out))


class IndependentSelection(unittest.TestCase):
    def test_strip_only_exact_matching_sample_prefix(self):
        self.assertEqual(V.cell_key("sample01", "sample01_ACGT-1"), ("sample01", "ACGT-1"))
        self.assertEqual(V.cell_key("sample01", "sample02_ACGT-1"), ("sample01", "sample02_ACGT-1"))
        self.assertEqual(V.cell_key("sample01", "ACGT-1"), ("sample01", "ACGT-1"))

    def make_metadata(self, directory, alter=None):
        fields = ["sample_name", "sex", "condition", "seurat_clusters", "paper_clusters",
                  "paper_clusters_short", "UMAP_1", "UMAP_2", "barcode"]
        rows = [["sample01", "male", "SNP", "0", "C0: CD4+ TN RTE", "C0", "0", "0", "sample01_shared-1"],
                ["sample02", "male", "noSNP", "1", "C1: CD4+ TN mature", "C1", "0", "0", "sample02_shared-1"],
                ["sample01", "male", "SNP", "3", "C3: CD8+ TN", "C3", "0", "0", "sample01_other-1"]]
        if alter:
            alter(rows)
        path = Path(directory) / "meta.tsv"
        with path.open("w",newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(fields)
            writer.writerows(rows)
        mapping = {"sample01": {"group":"SNP", "donor":"one"}, "sample02": {"group":"noSNP", "donor":"two"}}
        barcodes = {"sample01": ["shared-1", "other-1", "source-only-1"], "sample02": ["shared-1"]}
        return path, mapping, barcodes

    def test_repeated_bare_barcode_distinct_full_keys(self):
        with tempfile.TemporaryDirectory() as tmp:
            keys, selection, totals, states = V.metadata_selection(*self.make_metadata(tmp))
        self.assertEqual(len(keys), 3)
        self.assertEqual(len(selection), 2)
        self.assertEqual({r["source_column_1based"] for r in selection}, {1})
        self.assertEqual({r["sample"] for r in selection}, {"sample01", "sample02"})
        self.assertFalse(any(r["barcode"] == "source-only-1" for r in selection))

    def test_cross_sample_prefix_not_removed(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = self.make_metadata(tmp, lambda rows: rows[0].__setitem__(8,"sample02_shared-1"))
            with self.assertRaisesRegex(ValueError, "own sample"):
                V.metadata_selection(*args)

    def test_changed_crosswalk_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = self.make_metadata(tmp, lambda rows: rows[0].__setitem__(5, "C1"))
            with self.assertRaisesRegex(ValueError, "crosswalk"):
                V.metadata_selection(*args)

    def test_duplicate_full_key_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            args = self.make_metadata(tmp, lambda rows: rows.append(rows[0].copy()))
            with self.assertRaisesRegex(ValueError, "duplicate"):
                V.metadata_selection(*args)

    def test_wrong_genotype_and_missing_retained_not_ignored(self):
        for column,value in [(2,"noSNP"),(8,"sample01_missing-1")]:
            with tempfile.TemporaryDirectory() as tmp:
                args = self.make_metadata(tmp, lambda rows: rows[0].__setitem__(column,value))
                with self.assertRaises(ValueError):
                    V.metadata_selection(*args)


class ProductionAdapterSynthetic(unittest.TestCase):
    """Import author functions only as systems under test, never as an oracle."""
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location("production_bach2_test_adapter", ROOT / "scripts/analyze_psc_bach2_naive_rna.py")
        cls.A = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.A)

    def compare_welch(self, a, b):
        from decimal import Decimal
        got = self.A.welch({"SNP": [Decimal.from_float(float(x)) for x in a],
                            "noSNP": [Decimal.from_float(float(x)) for x in b]})
        expected = V.independent_welch(a, b)
        self.assertAlmostEqual(float(got["delta"]), expected["delta"], delta=max(1e-14,abs(expected["delta"])*2e-14))
        if expected["inference_status"] == "point_only":
            self.assertIsNone(got["CI95"])
            self.assertIsNone(got["P"])
        else:
            self.assertAlmostEqual(float(got["SE"]), expected["se"], delta=abs(expected["se"])*2e-14)
            self.assertAlmostEqual(float(got["df"]), expected["df"], places=13)
            self.assertAlmostEqual(float(got["P"]), expected["p_two_sided"], delta=5e-12)
            # Student-t ppf implementations use finite numerical approximation.
            for x,y in zip(got["CI95"],expected["ci95"]):
                self.assertAlmostEqual(float(x), y, delta=max(1e-14,abs(y)*5e-9,expected["se"]*5e-9))
        return got

    def test_independent_welch_parity_finite_grid(self):
        cases = [([0,1,2,3],[1,2,3,4]), ([0,2,6,10],[.1,.2,.5,.9]),
                 ([.1]*4,[0,1,2,5]), ([0]*4,[0,0,0,1]), ([.1]*4,[.1]*4),
                 ([.1]*4,[.2]*4), ([0]*4,[1]*4),
                 ([.1,.1,.1,math.nextafter(.1, math.inf)],[.1]*4),
                 ([0,1e-200,2e-200,3e-200],[1e-200]*4)]
        for a,b in cases:
            with self.subTest(a=a,b=b):
                self.compare_welch(a,b)

    def test_unequal_donor_arm_rejected(self):
        from decimal import Decimal
        with self.assertRaises(self.A.GateError):
            self.A.welch({"SNP":[Decimal(0)]*3,"noSNP":[Decimal(0)]*4})

    def test_nonfinite_score_rejected(self):
        from decimal import Decimal
        for bad in [Decimal("NaN"),Decimal("Infinity")]:
            with self.assertRaises(self.A.GateError):
                self.A.welch({"SNP":[Decimal(0)]*3+[bad],"noSNP":[Decimal(0)]*4})

    def test_transform_matches_independent_decimal(self):
        for b,r in [(0,19),(1,1000000),(3,1000000),(2,19),(17,17),(9007199254740993,10**20),(1,10**100)]:
            # Oracle precision 90 is enough except the deliberately extreme ratio,
            # which is checked as strictly positive, without a clamp to zero.
            if r == 10**100:
                self.assertGreater(self.A.score_units(b,r),0)
            else:
                self.assertAlmostEqual(float(self.A.score_units(b,r)),V.independent_y(b,r),places=14)
        for b,r in [(0,0),(1,-1),(-1,2),(3,2),(1.0,3),(True,3)]:
            with self.assertRaises(self.A.GateError):
                self.A.score_units(b,r)

    def sparse(self, entries=None):
        import hashlib
        if entries is None:
            entries = [(1,1,10),(2,1,2),(4,1,9999),(2,2,10000),(3,3,7),(4,3,888),(2,4,8)]
        body = ("%%MatrixMarket matrix coordinate integer general\n"
                '%metadata_json: {"software_version": "cellranger-6.1.1", "format_version": 2}\n'
                f"4 4 {len(entries)}\n" + "".join(f"{r} {c} {v}\n" for r,c,v in entries)).encode()
        pin = {"rows":4,"columns":4,"bytes":len(body),"sha256":hashlib.sha256(body).hexdigest()}
        return body,pin

    def test_sparse_same_selected_cells_all_rna_and_no_adt(self):
        import io
        body,pin = self.sparse()
        got = self.A.stream_selected_matrix(io.BytesIO(body),pin,[1,3],[1,2,3],target_row=2)
        self.assertEqual(got,{"B":2,"R":19,"B_positive_cells":1,"selected_cells":2})

    def test_zero_target_and_all_zero_selected_donor_are_retained_until_R_gate(self):
        import io
        body,pin = self.sparse([(1,1,3),(4,1,900),(3,3,7),(4,3,500)])
        got = self.A.stream_selected_matrix(io.BytesIO(body),pin,[1,3],[1,2,3],target_row=2)
        self.assertEqual((got["B"],got["R"],got["B_positive_cells"]),(0,10,0))
        self.assertEqual(self.A.score_units(got["B"],got["R"]),0)
        body,pin = self.sparse([(4,1,900),(2,2,3),(4,3,500)])
        got = self.A.stream_selected_matrix(io.BytesIO(body),pin,[1,3],[1,2,3],target_row=2)
        with self.assertRaises(self.A.GateError):
            self.A.score_units(got["B"],got["R"])

    def test_stream_integer_precision_duplicate_and_malformed(self):
        import io
        body,pin = self.sparse([(2,1,"9007199254740993")])
        got = self.A.stream_selected_matrix(io.BytesIO(body),pin,[1],[1,2,3],target_row=2)
        self.assertEqual(got["B"],9007199254740993)
        for entries in [[(2,1,1),(2,1,2)],[(2,1,1),(1,1,2)],[(2,1,"0.5")],[(2,1,-1)],[(2,1,0)],[(2,5,1)]]:
            with self.subTest(entries=entries):
                body,pin = self.sparse(entries)
                with self.assertRaises(self.A.GateError):
                    self.A.stream_selected_matrix(io.BytesIO(body),pin,[1,3],[1,2,3],target_row=2)

    def test_stream_pin_and_axis_fail_closed(self):
        import io
        body,pin = self.sparse()
        for altered in [{**pin,"sha256":"0"*64},{**pin,"bytes":len(body)+1},{**pin,"rows":5}]:
            with self.assertRaises(self.A.GateError):
                self.A.stream_selected_matrix(io.BytesIO(body),altered,[1,3],[1,2,3],target_row=2)
        for columns,rna in [([1,1],[1,2,3]),([1,3],[1,3]),([0],[1,2,3])]:
            with self.assertRaises(self.A.GateError):
                self.A.stream_selected_matrix(io.BytesIO(body),pin,columns,rna,target_row=2)

    def make_donors(self):
        selection = {"donors":[]}
        measurements = {}
        for i,b in enumerate([0,1,3,7,2,4,9,12],1):
            sample = f"sample{i:02d}"
            selection["donors"].append({"sample":sample,"donor":f"synthetic-{i}","group":"SNP" if i<=4 else "noSNP",
                "selected":[{"column_1based":1},{"column_1based":2}],"state_counts":{"C0":1,"C1":1}})
            measurements[sample] = {"B":b,"R":1000000,"B_positive_cells":int(b>0),"selected_cells":2}
        return measurements,selection

    def test_donor_pipeline_independent_scores_and_all8_loo(self):
        measurements,selection = self.make_donors()
        got = self.A.analyze_donors(measurements,selection)
        expected = [V.independent_y(measurements[f"sample{i:02d}"]["B"],1000000) for i in range(1,9)]
        for row,value in zip(got["donors"],expected):
            self.assertAlmostEqual(float(row["Y"]),value,places=14)
        loo = V.independent_loo(expected[:4],expected[4:])
        self.assertEqual(len(got["LOO_private"]),8)
        for row,value in zip(got["LOO_private"],loo):
            self.assertEqual(set(row),{"omitted_sample","delta"})
            self.assertAlmostEqual(float(row["delta"]),value,places=14)
        self.assertEqual(got["zero_B_donors"],1)
        self.assertEqual(got["LOO_summary"]["count"],8)

    def test_missing_donor_changed_membership_or_invalid_R_blocks_all(self):
        for change in ("missing","cells","R","group"):
            measurements,selection = self.make_donors()
            if change == "missing": del measurements["sample01"]
            if change == "cells": measurements["sample01"]["selected_cells"] = 1
            if change == "R": measurements["sample01"]["R"] = 0
            if change == "group": selection["donors"][0]["group"] = "noSNP"
            with self.subTest(change=change),self.assertRaises(self.A.GateError):
                self.A.analyze_donors(measurements,selection)

    def fake_plan(self):
        return {"selection_sha256":"a"*64,"source_manifest_sha256":"b"*64,
                "implementation":[{"sha256":"c"*64}],"runtime":{"synthetic":True},
                "public_schema_sha256":"d"*64}

    def test_public_projection_drops_all_unapproved_keys(self):
        measurements,selection = self.make_donors()
        private = self.A.analyze_donors(measurements,selection)
        private["extra"] = {"sample01":"CANARY_DONOR_SECRET"}
        private["inference"]["secret"] = "CANARY_DONOR_SECRET"
        private["inference"]["groups"]["SNP"]["donor"] = "CANARY_DONOR_SECRET"
        private["LOO_summary"]["donors"] = "CANARY_DONOR_SECRET"
        public = self.A.public_projection(private,self.fake_plan(),"e"*64)
        self.assertNotIn("CANARY_DONOR_SECRET",json.dumps(public))
        self.assertNotIn("sample01",json.dumps(public))
        self.assertEqual(set(public),set(self.A.PUBLIC_SCHEMA))

    def test_public_projection_rejects_nested_private_objects_in_scalar_slots(self):
        measurements,selection = self.make_donors()
        for slot in ("CI95","delta"):
            private = self.A.analyze_donors(measurements,selection)
            private["inference"][slot] = [{"sample01":"CANARY_DONOR_SECRET"}]
            with self.subTest(slot=slot),self.assertRaises(self.A.GateError):
                self.A.public_projection(private,self.fake_plan(),"e"*64)

    def test_fresh_output_rejects_reuse_traversal_alias_link_and_source(self):
        with tempfile.TemporaryDirectory() as tmp, mock.patch.object(self.A,"ROOT",Path(tmp)):
            base = Path(tmp)/self.A.BASE
            (base/"runs").mkdir(parents=True)
            fresh = self.A.BASE+"/runs/synthetic-one"
            created = self.A.fresh_directory(fresh,create=True)
            self.assertTrue(created.is_dir())
            for bad in [fresh,self.A.BASE+"/runs/../escape",self.A.BASE+"/runs/a/b",self.A.BASE+"/runs/./other",
                        self.A.BASE+"/runs/other/", "/tmp/absolute", "data/raw/source", self.A.BASE+"/runs/a\\b"]:
                with self.subTest(bad=bad),self.assertRaises(self.A.GateError):
                    self.A.fresh_directory(bad,create=True)
            (base/"runs"/"linked").symlink_to(created,target_is_directory=True)
            with self.assertRaises(self.A.GateError):
                self.A.fresh_directory(self.A.BASE+"/runs/linked")
            (base/"runs").rename(base/"old-runs")
            (base/"runs").symlink_to(base/"old-runs",target_is_directory=True)
            with self.assertRaises(self.A.GateError):
                self.A.fresh_directory(self.A.BASE+"/runs/parentlink")

    def test_root_authorization_all_bindings_and_external_sha_required(self):
        plan = self.fake_plan()
        output = self.A.BASE+"/runs/synthetic-auth"
        auth = {"schema_version":1,"status":self.A.AUTH_STATUS,"expression_authorized":True,"authorizer":"root",
                "one_fixed_extraction_only":True,"output_directory":output,"plan_sha256":self.A.digest(plan),
                "selection_sha256":plan["selection_sha256"],"source_manifest_sha256":plan["source_manifest_sha256"],
                "implementation_sha256":plan["implementation"][0]["sha256"],"runtime_sha256":self.A.digest(plan["runtime"]),
                "public_schema_sha256":plan["public_schema_sha256"],"candidate_sha256":self.A.CANDIDATE_SHA}
        with tempfile.TemporaryDirectory() as tmp,mock.patch.object(self.A,"ROOT",Path(tmp)):
            raw_sha = self.A.digest(auth)
            self.A.check_authorization(auth,raw_sha,raw_sha,plan,output)
            for key in auth:
                changed = {**auth,key:None}
                with self.subTest(key=key),self.assertRaises(self.A.GateError):
                    self.A.check_authorization(changed,self.A.digest(changed),self.A.digest(changed),plan,output)
            with self.assertRaises(self.A.GateError):
                self.A.check_authorization(auth,"0"*64,raw_sha,plan,output)
            with self.assertRaises(self.A.GateError):
                self.A.check_authorization(auth,None,raw_sha,plan,output)

            for key,value in [("schema_version",True),("expression_authorized",1),("one_fixed_extraction_only",1)]:
                changed={**auth,key:value}
                with self.subTest(boolean_alias=key),self.assertRaises(self.A.GateError):
                    self.A.check_authorization(changed,self.A.digest(changed),self.A.digest(changed),plan,output)

    def test_no_authorization_cannot_reach_real_matrix_or_validation(self):
        with mock.patch.object(self.A,"validate_plan",side_effect=AssertionError("must not validate without authorization")):
            with self.assertRaises(self.A.GateError):
                self.A.execute("not-used",None,None,None)
            with self.assertRaises(self.A.GateError):
                self.A.execute("not-used",self.A.BASE+"/fake-auth.json","0"*64,self.A.BASE+"/runs/test")

    def test_default_main_validate_has_no_execute_or_output_action(self):
        plan=self.fake_plan()
        with mock.patch.object(self.A,"validate_plan",return_value=(plan,{})) as validate, \
             mock.patch.object(self.A,"execute",side_effect=AssertionError("no execution")), \
             mock.patch.object(self.A,"compile_plan",side_effect=AssertionError("no plan writing")):
            self.assertEqual(self.A.main([]),0)
            validate.assert_called_once()
            self.assertEqual(self.A.main(["--output-directory",self.A.BASE+"/runs/ignored"]),2)

    def test_loo_summary_complete_sign_changes_without_tests(self):
        for a,b,counts,fullsign,reversals in [([0,0,0,10],[1]*4,(7,1,0),1,1),
                                             ([0,0,0,4],[1]*4,(3,1,4),0,0)]:
            records=[{"sample":f"s{i}","group":"SNP" if i<4 else "noSNP","Y":str(y)}
                     for i,y in enumerate(a+b)]
            rows,summary=self.A.leaveouts(records)
            expected=V.independent_loo(a,b)
            self.assertEqual(len(rows),8)
            self.assertEqual({row["omitted_sample"] for row in rows},{f"s{i}" for i in range(8)})
            for row,value in zip(rows,expected):
                self.assertAlmostEqual(float(row["delta"]),value,places=14)
            self.assertEqual(tuple(summary[k] for k in ("positive_count","negative_count","zero_count")),counts)
            self.assertEqual(summary["full_point_sign"],fullsign)
            self.assertEqual(summary["strict_sign_reversal_count"],reversals)
            self.assertTrue(all(set(row)=={"omitted_sample","delta"} for row in rows))

    def test_numeric_tail_underflow_never_forced_zero(self):
        from decimal import Decimal
        got=self.A.welch({"SNP":[Decimal(0),Decimal("1e-200"),Decimal("2e-200"),Decimal("3e-200")],
                          "noSNP":[Decimal(1)]*4})
        self.assertIsNone(got["P"])
        self.assertIsNotNone(got["CI95"])
        self.assertEqual(got["inference_status"],"CI_available_P_unavailable_numeric_tail")

    def test_nonfinite_quantile_leaves_point_only(self):
        from decimal import Decimal
        with mock.patch("scipy.stats.t.ppf",return_value=math.inf):
            got=self.A.welch({"SNP":[Decimal(0),Decimal(1),Decimal(2),Decimal(3)],"noSNP":[Decimal(0)]*4})
        self.assertEqual(float(got["delta"]),1.5)
        self.assertIsNone(got["P"])
        self.assertIsNone(got["CI95"])

    def test_exact_product_equal_point_zero_with_positive_variance(self):
        measurements,selection=self.make_donors()
        b_values=[1,2,3,5,0,3,5,5]
        for i,b in enumerate(b_values,1):
            measurements[f"sample{i:02d}"]={"B":b,"R":1000000,"B_positive_cells":int(b>0),"selected_cells":2}
        got=self.A.analyze_donors(measurements,selection)
        from decimal import Decimal
        self.assertEqual(Decimal(got["inference"]["delta"]),0)
        self.assertGreater(Decimal(got["inference"]["SE"]),0)
        self.assertEqual(Decimal(got["inference"]["P"]),1)
        self.assertEqual(got["LOO_summary"]["full_point_sign"],0)
        self.assertEqual(got["LOO_summary"]["strict_sign_reversal_count"],0)

    def test_genuine_fourth_order_unit_product_difference_not_clamped(self):
        from decimal import Decimal,localcontext
        measurements,selection=self.make_donors()
        offsets=[1,-1,8,-8,4,-4,7,-7]
        pairs=[(10**93+offset,10**99) for offset in offsets]
        for i,(b,r) in enumerate(pairs,1):
            measurements[f"sample{i:02d}"]={"B":b,"R":r,"B_positive_cells":1,"selected_cells":2}
        expected=V.independent_unit_points({"SNP":pairs[:4],"noSNP":pairs[4:]})
        got=self.A.analyze_donors(measurements,selection)
        point=Decimal(got["inference"]["delta"])
        self.assertLess(point,0)
        self.assertGreater(Decimal(got["inference"]["SE"]),0)
        self.assertEqual(got["LOO_summary"]["full_point_sign"],-1)
        with localcontext() as context:
            context.prec=1200
            self.assertLess(abs((point-expected["delta"])/expected["delta"]),Decimal("1e-300"))
            for row,value in zip(got["LOO_private"],expected["leaveouts"]):
                if value:
                    self.assertLess(abs((Decimal(row["delta"])-value)/value),Decimal("1e-300"))
                else:
                    self.assertEqual(Decimal(row["delta"]),0)

    def test_native_runtime_record_checks_installed_body(self):
        import base64,hashlib
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp)
            body=b"synthetic-installed-module\n"
            (directory/"module.py").write_bytes(body)
            digest=base64.urlsafe_b64encode(hashlib.sha256(body).digest()).decode().rstrip("=")
            record=f"module.py,sha256={digest},{len(body)}\nexample.dist-info/RECORD,,\n"
            class Distribution:
                version="synthetic"
                def read_text(self,name): return record
                def locate_file(self,name): return directory/name
            self.assertEqual(self.A.verify_distribution(Distribution())["verified_installed_files"],1)
            (directory/"module.py").write_bytes(body+b"tampered")
            with self.assertRaises(self.A.GateError):
                self.A.verify_distribution(Distribution())

    def test_entire_synthetic_analysis_ignores_ambient_decimal_context(self):
        from decimal import localcontext,ROUND_DOWN
        measurements,selection=self.make_donors()
        expected=self.A.analyze_donors(measurements,selection)
        with localcontext() as context:
            context.prec=6
            context.rounding=ROUND_DOWN
            context.Emin=-99
            got=self.A.analyze_donors(measurements,selection)
        self.assertEqual(expected,got)

    def test_wrong_authorization_file_role_rejected_before_any_open(self):
        wrong=["work/psc-bach2-cohort-qualification/sample01/members/matrix.txt",
               "data/raw/psc-bach2-qualification/cohort/sample01_filtered_feature_bc_matrix.tar.gz",
               "data/derived/psc-bach2-cohort-qualification.json",
               self.A.BASE+"/fake-auth.json","work/psc-bach2-naive-rna-freeze/../escape.json",
               "work/psc-bach2-naive-rna-freeze/nested/auth.json","work/psc-bach2-naive-rna-freeze/allowed.txt"]
        with mock.patch.object(Path,"open",side_effect=AssertionError("file-role rejection must precede open")):
            for path in wrong:
                with self.subTest(path=path),self.assertRaises(self.A.GateError):
                    self.A.execute(self.A.DEFAULT_PLAN,path,"0"*64,self.A.BASE+"/runs/synthetic-role")

    def test_auth_size_and_symlink_rejected_before_body_open(self):
        with tempfile.TemporaryDirectory() as tmp,mock.patch.object(self.A,"ROOT",Path(tmp)):
            directory=Path(tmp)/"work/psc-bach2-naive-rna-freeze"
            directory.mkdir(parents=True)
            (directory/"large.json").write_bytes(b"x"*8193)
            (directory/"target.json").write_bytes(b"{}")
            (directory/"link.json").symlink_to(directory/"target.json")
            with mock.patch.object(Path,"open",side_effect=AssertionError("size/link gate before body")):
                for name in ("large.json","link.json"):
                    with self.subTest(name=name),self.assertRaises(self.A.GateError):
                        self.A.read_root_authorization("work/psc-bach2-naive-rna-freeze/"+name)

    def test_wrong_plan_directory_role_rejected_before_any_open(self):
        wrong=["data/raw/psc-bach2-qualification", "work/psc-bach2-cohort-qualification/sample01/members",
               self.A.BASE+"/runs", self.A.BASE+"/runs/old-run", self.A.BASE+"/../old-result"]
        with mock.patch.object(Path,"open",side_effect=AssertionError("plan-role rejection must precede open")):
            for path in wrong:
                with self.subTest(path=path),self.assertRaises(self.A.GateError):
                    self.A.validate_plan(path)

    def test_completion_inventory_binds_exact_three_payloads_and_plan(self):
        import hashlib
        plan=self.fake_plan()
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp)
            bodies={"measurements.private.json":b'{"synthetic":true}\n',
                    "public-aggregate.json":b'{"synthetic":true}\n',"authorization.used.json":b'{}\n'}
            for name,body in bodies.items(): (directory/name).write_bytes(body)
            self.A.write_completion(directory,plan,"e"*64,"synthetic-output")
            record=json.loads((directory/"completion.json").read_text())
            self.assertEqual(set(p.name for p in directory.iterdir()),set(bodies)|{"completion.json"})
            self.assertEqual(set(record["allowed_inventory"]),set(bodies)|{"completion.json"})
            self.assertEqual(record["payloads"],{name:{"bytes":len(body),"sha256":hashlib.sha256(body).hexdigest()} for name,body in bodies.items()})
            self.assertEqual(record["provenance"]["plan_sha256"],hashlib.sha256((json.dumps(plan,sort_keys=True,indent=2,allow_nan=False)+"\n").encode()).hexdigest())
            self.A.validate_completed_run(directory,plan,"e"*64,"synthetic-output")
            for alternate in [(plan,"f"*64,"synthetic-output"),({**plan,"changed":True},"e"*64,"synthetic-output"),(plan,"e"*64,"other-output")]:
                with self.assertRaises(self.A.GateError): self.A.validate_completed_run(directory,*alternate)
            for change in ("missing","extra","modified","link"):
                if change=="missing": (directory/"public-aggregate.json").unlink()
                if change=="extra": (directory/"extra.txt").write_text("synthetic")
                if change=="modified": (directory/"public-aggregate.json").write_bytes(b"changed")
                if change=="link":
                    (directory/"public-aggregate.json").unlink()
                    (directory/"public-aggregate.json").symlink_to(directory/"measurements.private.json")
                with self.subTest(change=change),self.assertRaises(self.A.GateError):
                    self.A.validate_completed_run(directory,plan,"e"*64,"synthetic-output")
                if change=="extra": (directory/"extra.txt").unlink()
                if change=="link": (directory/"public-aggregate.json").unlink()
                (directory/"public-aggregate.json").write_bytes(bodies["public-aggregate.json"])

    def test_failed_completion_rename_removes_both_markers(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory=Path(tmp)
            for name in self.A.RUN_PAYLOADS: (directory/name).write_bytes(b"{}\n")
            with mock.patch.object(Path,"replace",side_effect=OSError("synthetic interrupted completion")):
                with self.assertRaises(OSError): self.A.write_completion(directory,self.fake_plan(),"e"*64,"synthetic-output")
            self.assertFalse((directory/"completion.json").exists())
            self.assertFalse((directory/"completion.pending").exists())

    def test_execute_failure_before_or_after_completion_cannot_be_complete(self):
        for stage in ("write_completion","validate_completed_run"):
            with self.subTest(stage=stage),tempfile.TemporaryDirectory() as tmp,mock.patch.object(self.A,"ROOT",Path(tmp)):
                plan={**self.fake_plan(),"source_manifest":{"matrix_expected_pins_NOT_read":[]}}
                measurements,selection=self.make_donors()
                private=self.A.analyze_donors(measurements,selection)
                output=self.A.BASE+"/runs/failure-fixture"
                auth={"schema_version":1,"status":self.A.AUTH_STATUS,"expression_authorized":True,"authorizer":"root",
                      "one_fixed_extraction_only":True,"output_directory":output,"plan_sha256":self.A.digest(plan),
                      "selection_sha256":plan["selection_sha256"],"source_manifest_sha256":plan["source_manifest_sha256"],
                      "implementation_sha256":plan["implementation"][0]["sha256"],"runtime_sha256":self.A.digest(plan["runtime"]),
                      "public_schema_sha256":plan["public_schema_sha256"],"candidate_sha256":self.A.CANDIDATE_SHA}
                relative="work/psc-bach2-naive-rna-freeze/synthetic.json"
                auth_path=Path(tmp)/relative
                auth_path.parent.mkdir(parents=True)
                auth_path.write_bytes(self.A.canonical(auth))
                with mock.patch.object(self.A,"validate_plan",return_value=(plan,selection)), \
                     mock.patch.object(self.A,"analyze_donors",return_value=private), \
                     mock.patch.object(self.A,stage,side_effect=self.A.GateError("synthetic output failure")):
                    with self.assertRaises(self.A.GateError): self.A.execute(self.A.DEFAULT_PLAN,relative,self.A.digest(auth),output)
                directory=Path(tmp)/output
                self.assertTrue((directory/"failure.json").is_file())
                self.assertFalse((directory/"completion.json").exists())
                self.assertFalse((directory/"completion.pending").exists())
                with self.assertRaises(self.A.GateError): self.A.validate_completed_run(directory,plan,self.A.digest(auth),output)


if __name__ == "__main__":
    unittest.main()
