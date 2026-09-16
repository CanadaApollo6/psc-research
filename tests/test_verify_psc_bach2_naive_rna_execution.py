"""Synthetic only. Production helpers are a fixture producer, NEVER the auditor.

Full native axes/55,460 metadata cells/16,659 selected/eight source labels are
fabricated in TemporaryDirectory. No real matrix or real endpoint is opened.
"""
import copy
import csv
from decimal import Decimal, Context, localcontext
from fractions import Fraction
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]

def load(relative, name):
    spec = importlib.util.spec_from_file_location(name, ROOT/relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

e = load("scripts/verify_psc_bach2_naive_rna_execution.py", "execution_audit_synthetic")
a = load("scripts/analyze_psc_bach2_naive_rna.py", "synthetic_fixture_producer_only")
v = load("scripts/verify_psc_bach2_naive_rna.py", "synthetic_annotation_primitives")


def simple_selection():
    return {"donors": [{"sample": s, "donor": "SYNTHETIC_DONOR_"+str(i), "group": "SNP" if i < 4 else "noSNP",
        "retained_count": 2, "state_counts": {"C0": 1, "C1": 1},
        "selected": [{"barcode": "AAAA-1", "column_1based": 1, "state": "C0"},
                     {"barcode": "CCCC-1", "column_1based": 2, "state": "C1"}]} for i, s in enumerate(e.SAMPLES)],
        "retained_keys": [[s, b] for s in e.SAMPLES for b in ("AAAA-1", "CCCC-1")]}


def units(values=None, selection=None):
    selection = selection or simple_selection()
    values = values or [(i, 100+i) for i in range(8)]
    return {s: {"B": b, "R": r, "B_positive_cells": int(b > 0),
                "selected_cells": len(selection["donors"][i]["selected"])}
            for i, (s, (b, r)) in enumerate(zip(e.SAMPLES, values))}


def fixture_plan():
    return {"selection_sha256": "1"*64, "source_manifest_sha256": "2"*64,
            "implementation": [{"sha256": "3"*64}], "runtime": {"synthetic": True},
            "public_schema_sha256": e.digest(e.PUBLIC_SCHEMA)}


def payloads(measurements, selection=None, plan=None):
    selection, plan = selection or simple_selection(), plan or fixture_plan()
    private = a.analyze_donors(measurements, selection)
    public = a.public_projection(private, plan, "4"*64)
    return private, public


def matrix(entries=None, cols=3, nnz=None):
    entries = entries if entries is not None else [(1,1,7), (12314,1,2), (36639,1,999), (2,2,11),
                                                   (12314,2,3), (1,3,9000)]
    lines = "".join(f"{r} {c} {x}\n" for r,c,x in entries)
    body = (e.HEADER+"\n%metadata_json: "+json.dumps(e.EXPORTER)+f"\n36639 {cols} {len(entries) if nnz is None else nnz}\n"+lines).encode()
    return body, {"rows":36639, "columns":cols, "bytes":len(body), "sha256":e.sha(body)}


class NumericAudit(unittest.TestCase):
    def audit(self, pairs):
        m = units(pairs)
        private, public = payloads(m)
        result = e.compare_payloads(private, public, m, simple_selection(), fixture_plan(), "4"*64)
        self.assertTrue(result["all_aggregate_numbers_status_QC_provenance_agree"])
        return private, public

    def test_ordinary_full_private_and_public_numbers(self):
        self.audit([(i,100+i) for i in range(8)])

    def test_zero_B_retained_one_constant_arm(self):
        private,_ = self.audit([(0,10)]*4 + [(i,100) for i in range(1,5)])
        self.assertEqual(Decimal(private["inference"]["df"]),3)
        self.assertEqual(private["zero_B_donors"],4)

    def test_both_constant_zero_or_distinct_point_only(self):
        for pairs in ([(0,10)]*8, [(1,10)]*4+[(2,10)]*4):
            with self.subTest(pairs=pairs):
                private,_ = self.audit(pairs)
                self.assertEqual(private["inference"]["inference_status"],"unavailable_zero_SE")
                self.assertIsNone(private["inference"]["P"])

    def test_exact_product_cancellation_positive_variance(self):
        q = [2,2,8,8,4,4,4,4]
        private,_ = self.audit([(x-1,1_000_000) for x in q])
        self.assertEqual(private["inference"]["delta"],"0")
        self.assertGreater(Decimal(private["inference"]["SE"]),0)
        self.assertEqual(private["inference"]["P"],"1.0")

    def test_tiny_variance_near_nonzero_large_denominators(self):
        n=10**99
        private,_ = self.audit([(n+i,2*n+i+1) for i in range(8)])
        self.assertGreater(Decimal(private["inference"]["SE"]),0)
        self.assertLess(Decimal(private["inference"]["SE"]),Decimal("1e-90"))

    def test_smallest_ratio_spacing_below_float_resolution(self):
        n=10**106
        private,_ = self.audit([(n+i,2*n+2*i+1) for i in range(8)])
        self.assertGreater(Decimal(private["inference"]["SE"]),0)
        self.assertLess(Decimal(private["inference"]["SE"]),Decimal("1e-205"))

    def test_tiny_fourth_order_cancelling_delta_preserved(self):
        n=10**94
        # Symmetric pairs share first three power sums: 1**2+8**2=4**2+7**2.
        # Their exact four-factor products differ only at fourth order.
        left=[-8,-1,1,8]; right=[-7,-4,4,7]
        pairs=[(n+i,2*n) for i in left+right]
        private,_=self.audit(pairs)
        self.assertNotEqual(Decimal(private["inference"]["delta"]),0)
        self.assertLess(abs(Decimal(private["inference"]["delta"])),Decimal("1e-370"))

    def test_positive_decimal_tail_when_native_probability_underflows(self):
        n=10**106
        pairs=[(n+i,2*n+2*i+1) for i in range(4)] + [(1,n)]*4
        private,_=self.audit(pairs)
        self.assertEqual(private["inference"]["inference_status"],"CI_available_P_unavailable_numeric_tail")
        ref=e.reference_endpoint(units(pairs),simple_selection())
        self.assertGreater(ref["student"]["P"],0)
        self.assertLess(ref["student"]["P"],Decimal("1e-308"))

    def test_nonfinite_native_quantile_corrobated_point_only(self):
        m=units()
        with patch("scipy.special.stdtrit", return_value=float("nan")):
            private,public=payloads(m)
            e.compare_payloads(private,public,m,simple_selection(),fixture_plan(),"4"*64)
            self.assertEqual(private["inference"]["inference_status"],"unavailable_nonfinite_inference")

    def test_probability_reference_matches_known_student_value(self):
        value=e.beta_tail(Decimal(1),Decimal(3))
        self.assertLess(abs(value-Decimal("0.3910022189557706")),Decimal("1e-14"))
        self.assertGreater(e.beta_tail(Decimal("1e200"),Decimal(3)),0)

    def test_strict_zero_sign_relative_no_absolute_epsilon(self):
        for actual,reference in (("0",Decimal("1e-1000")),("-1e-1000",Decimal("1e-1000")),("1e-1000",Decimal(0))):
            with self.subTest(actual=actual),self.assertRaises(e.AuditError):e.agree(actual,reference)
        e.agree("1e-1000",Decimal("1e-1000"))
        with self.assertRaises(e.AuditError):e.agree("1.0000001e-1000",Decimal("1e-1000"),"1e-10")

    def test_root_fixed_relative_boundary_has_no_absolute_floor(self):
        e.agree("1.00000000005e-1000",Decimal("1e-1000"))
        with self.assertRaises(e.AuditError):e.agree("1.0000000002e-1000",Decimal("1e-1000"))
        with self.assertRaises(e.AuditError):e.agree("1e-1000",Decimal(0))

    def test_4000_digit_log_covers_conservative_product_cancellation_bound(self):
        n=10**2663
        positive=e.direct_log2(Fraction(n+1,n))
        negative=e.direct_log2(Fraction(n,n+1))
        self.assertGreater(positive,0);self.assertLess(negative,0)
        with localcontext(Context(prec=4000)):
            self.assertLess(abs(positive/(-negative)-1),Decimal("1e-1000"))

    def test_standardized_CI_gate_rejects_error_even_at_tiny_scale(self):
        n=10**106;m=units([(n+i,2*n+2*i+1) for i in range(8)])
        private,public=payloads(m)
        ref=e.reference_endpoint(m,simple_selection())
        with localcontext(Context(prec=320)):
            private["inference"]["CI95"][0]=str(Decimal(private["inference"]["CI95"][0])-ref["SE"]*Decimal("1e-7"))
        with self.assertRaises(e.AuditError):e.inference_compare(private["inference"],ref)

    def test_all_eight_and_invalid_denominators_gate_entire_endpoint(self):
        for replacement in ({"B":1,"R":0},{"B":2,"R":1},{"B":-1,"R":10},{"B":1,"R":e.MAX_R}):
            m=units();m[e.SAMPLES[0]].update(replacement)
            with self.subTest(replacement=replacement),self.assertRaises(e.AuditError):e.reference_endpoint(m,simple_selection())
        m=units();del m[e.SAMPLES[0]]
        with self.assertRaises(e.AuditError):e.reference_endpoint(m,simple_selection())

    def test_private_quantity_Y_state_fraction_and_LOO_tampering(self):
        m=units();private,public=payloads(m)
        for path,value in ((["donors",0,"B_source_units"],"9"),(["donors",0,"Y"],"1"),
                           (["donors",0,"state_cell_fractions","C0"],"0.4"),
                           (["LOO_private",0,"delta"],"0"),(["inference","SE"],"0")):
            altered=copy.deepcopy(private);slot=altered
            for k in path[:-1]:slot=slot[k]
            slot[path[-1]]=value
            with self.subTest(path=path),self.assertRaises(e.AuditError):e.compare_payloads(altered,public,m,simple_selection(),fixture_plan(),"4"*64)

    def test_nested_public_private_leaks_refused(self):
        m=units();private,public=payloads(m)
        for path in ([],["groups","SNP"],["QC"],["provenance"],["leave_one_donor_out"]):
            altered=copy.deepcopy(public);slot=altered
            for k in path:slot=slot[k]
            slot["PRIVATE_DONOR"]="synthetic"
            with self.subTest(path=path),self.assertRaises(e.AuditError):e.compare_payloads(private,altered,m,simple_selection(),fixture_plan(),"4"*64)
        altered=copy.deepcopy(public);altered["groups"]["SNP"]["mean"]={"donor":"private"}
        with self.assertRaises(e.AuditError):e.public_shape(altered)

    def test_public_contrast_QC_provenance_and_complete_LOO_summary(self):
        m=units();private,public=payloads(m)
        for path,value in ((["contrast","P_two_sided"],"0"),(["QC","BACH2_positive_selected_cells"],0),
                           (["provenance","plan_sha256"],"0"*64),(["leave_one_donor_out","count"],7)):
            altered=copy.deepcopy(public);slot=altered
            for k in path[:-1]:slot=slot[k]
            slot[path[-1]]=value
            with self.subTest(path=path),self.assertRaises(e.AuditError):e.compare_payloads(private,altered,m,simple_selection(),fixture_plan(),"4"*64)

    def test_ambient_decimal_precision_does_not_change_reference(self):
        with localcontext(Context(prec=6)):
            self.audit([(i,100+i) for i in range(8)])


class StreamAudit(unittest.TestCase):
    def parse(self, entries=None, **kwargs):
        body,pin=matrix(entries,**kwargs)
        return e.stream_matrix(io.BytesIO(body),pin,[1,2],list(range(1,36602)))

    def test_full_RNA_axis_same_selected_cells_no_ADT_or_unselected_RNA(self):
        self.assertEqual(self.parse(),{"B":5,"R":23,"B_positive_cells":2,"selected_cells":2})

    def test_zero_B_retained_but_zero_R_unavailable(self):
        self.assertEqual(self.parse([(1,1,7)])["B"],0)
        with self.assertRaises(e.AuditError):self.parse([(36639,1,999)])

    def test_exact_integers_and_exponent_forms(self):
        self.assertEqual(e.integer_token("+0002.00e3"),2000)
        self.assertEqual(self.parse([(1,1,2**53+1)])["R"],2**53+1)
        self.assertEqual(e.integer_token("9e99"),9*10**99)

    def test_bad_tokens_include_invalid_ADT(self):
        for token in ("0","-1","1.1","NaN","Inf","1e100","1e99999999999","１","1"*101):
            with self.subTest(token=token),self.assertRaises(e.AuditError):self.parse([(1,1,1),(36639,1,token)])

    def test_bad_dimensions_nnz_bounds_order_duplicates(self):
        for entries in ([(1,1,1),(1,1,2)],[(2,1,1),(1,1,2)],[(36640,1,1)],[(1,0,1)]):
            with self.subTest(entries=entries),self.assertRaises(e.AuditError):self.parse(entries)
        with self.assertRaises(e.AuditError):self.parse(nnz=100)

    def test_complete_hash_bytes_header_comment_and_line_bounds(self):
        body,pin=matrix()
        for changed in (body+b"\n",body.replace(b"integer general",b"real general"),body.replace(b"cellranger-6.1.1",b"cellranger-3.0.2"),body+b"x"*513):
            p={**pin,"bytes":len(changed),"sha256":e.sha(changed)}
            with self.subTest(changed=changed[:50]),self.assertRaises(e.AuditError):e.stream_matrix(io.BytesIO(changed),p,[1,2],list(range(1,36602)))
        with self.assertRaises(e.AuditError):e.stream_matrix(io.BytesIO(body),{**pin,"sha256":"0"*64},[1,2],list(range(1,36602)))

    def test_duplicate_or_missing_selected_columns_and_RNA_axis(self):
        body,pin=matrix()
        for columns,rna in (([1,1],list(range(1,36602))),([],list(range(1,36602))),([1],[1,12314])):
            with self.subTest(columns=columns),self.assertRaises(e.AuditError):e.stream_matrix(io.BytesIO(body),pin,columns,rna)


class NativeFixture:
    """Generated sources, with native schema and the exact fixed axis sizes."""
    def __init__(self, root):
        self.root=Path(root)
        self.patches=[]
        def write(relative, body):
            path=self.root/relative;path.parent.mkdir(parents=True,exist_ok=True)
            path.write_bytes(body if isinstance(body,bytes) else body.encode())
            return {"path":relative,"bytes":path.stat().st_size,"sha256":e.sha(path.read_bytes())}
        self.write=write
        for relative in (".gitignore",a.CANDIDATE,a.REVIEW,a.DECISION,a.__file__.replace(str(ROOT)+"/",""),
                         "tests/test_psc_bach2_naive_rna.py",e.VERIFIER,"tests/test_verify_psc_bach2_naive_rna.py",e.ADAPTER,e.TESTS,e.DESIGN):
            write(relative,(ROOT/relative).read_bytes())
        axis=[(f"ENSG_SYNTHETIC_{i}",f"SYNTHETIC_{i}","Gene Expression" if i <=36601 else "Antibody Capture")
              for i in range(1,36640)]
        axis[12313]=e.TARGET
        feature_body=("\n".join("\t".join(row) for row in axis)+"\n").encode()
        refrows=[[e.TARGET[0],"60468",e.TARGET[1]]] + [[f"SYNTHETIC_REF_{i}",str(i),f"SYNTHETIC_SYMBOL_{i}"] for i in range(1,86411)]
        refrows += refrows[1:5338]
        reference="Gene stable ID\tNCBI gene (formerly Entrezgene) ID\tHGNC symbol\n" + "\n".join("\t".join(r) for r in refrows)+"\n[success]\n"
        refpin=write(a.REFERENCE,reference)
        def barcode(i):
            return "".join("ACGT"[(i >> shift)&3] for shift in range(18,-1,-2))+"-1"
        metaheader=["sample_name","sex","condition","seurat_clusters","paper_clusters","paper_clusters_short","UMAP_1","UMAP_2","barcode"]
        metabuf=io.StringIO();mw=csv.writer(metabuf,delimiter="\t",lineterminator="\n");mw.writerow(metaheader)
        sdrfheader=["Characteristics[individual]","Characteristics[sex]","Characteristics[disease]","Characteristics[genotype]","Factor Value[genotype]",
                    "Derived Array Data File","Derived Array Data File","Derived Array Data File"]
        sb=io.StringIO();sw=csv.writer(sb,delimiter="\t",lineterminator="\n");sw.writerow(sdrfheader)
        receipt_pins={};archive_pins={}
        self.expected={}
        for i,sample in enumerate(e.SAMPLES):
            group="SNP" if i<4 else "noSNP"
            cells=7414+int(i<3);retained=6932+int(i<4);c0=1165;c1=917+int(i<3);selected=c0+c1
            bars=[barcode(j) for j in range(cells)]
            sw.writerow([f"SYNTHETIC_DONOR_{i}","male","primary sclerosing cholangitis",group,group,
                         sample+"_filtered_feature_bc_matrix.tar.gz","SYNTHETIC_OTHER","SYNTHETIC_OTHER2"])
            for j in range(retained):
                short="C0" if j<c0 else "C1" if j<selected else "C2"
                label={"C0":"C0: CD4+ TN RTE","C1":"C1: CD4+ TN mature","C2":"C2: SYNTHETIC_UNSELECTED"}[short]
                mw.writerow([sample,"male",group,short[1:],label,short,"0","0",sample+"_"+bars[j]])
            base=f"work/psc-bach2-cohort-qualification/{sample}/members"
            entries=[]
            for col in range(1,selected+1):
                entries.append((1,col,1))
                if col==1:
                    entries.extend([(12314,col,i+1),(36639,col,999)])
            entries.extend([(1,cells,999999),(36639,cells,999999)])
            body,_=matrix(entries,cols=cells)
            for kind,data in (("features",feature_body),("barcodes",("\n".join(bars)+"\n").encode()),("matrix",body)):
                p=write(base+"/"+kind+".txt",data)
                receipt={"decoded_path":p["path"],"decoded_bytes":p["bytes"],"decoded_sha256":p["sha256"],
                         "nested_gzip_CRC_ISIZE_and_EOF_passed":True,"source_payload_matches_tar":True}
                rp=write(base+"/"+kind+"-integrity.json",e.canonical(receipt));receipt_pins[rp["path"]]=rp["sha256"]
            archive_path=f"data/raw/psc-bach2-qualification/{'counts' if i==4 else 'cohort'}/{sample}_filtered_feature_bc_matrix.tar.gz"
            # Intentionally not an archive: metadata review must NEVER open it.
            ap=write(archive_path,b"SYNTHETIC_OPAQUE_NOT_A_REAL_ARCHIVE\n"+sample.encode())
            archive_pins[archive_path]={k:ap[k] for k in ("bytes","sha256")}
            receipt={"complete":True,"body_bytes":ap["bytes"],"expected_bytes":ap["bytes"],"authorized_archive_bytes":ap["bytes"],
                     "archive_path":archive_path,"sha256":ap["sha256"],"url":"https://invalid.example/SYNTHETIC", "completed_at_utc":"SYNTHETIC_ONLY"}
            rp=write("data/raw/psc-bach2-qualification/counts/acquisition-receipt.json" if i==4 else archive_path+".receipt.json",e.canonical(receipt))
            receipt_pins[rp["path"]]=rp["sha256"]
            self.expected[sample]={"B":i+1,"R":selected+i+1,"selected_cells":selected,"B_positive_cells":1}
        metapin=write(a.METADATA,metabuf.getvalue());sdrfpin=write(a.SDRF,sb.getvalue())
        patches=[patch.object(a,"ROOT",self.root),patch.object(a,"META_SHA",metapin["sha256"]),
                 patch.object(a,"SDRF_SHA",sdrfpin["sha256"]),patch.object(a,"FEATURE_SHA",e.sha(feature_body)),
                 patch.object(a,"REFERENCE_SHA",refpin["sha256"]),
                 patch.object(v,"FIXED",{**v.FIXED,a.REFERENCE:refpin["sha256"]}),patch.object(v,"RECEIPT_PINS",receipt_pins),
                 patch.object(v,"ARCHIVE_RELEASE_PINS",archive_pins),
                 patch.object(v,"METADATA_PIN",(a.METADATA,metapin["sha256"],metapin["bytes"])),
                 patch.object(v,"SDRF_PIN",(a.SDRF,sdrfpin["sha256"],sdrfpin["bytes"])),patch.object(e,"load_verifier",return_value=v)]
        for p in patches:p.start();self.patches.append(p)
        self.plan,self.selection=a.build_source_plan()
        self.plan_dir=a.DEFAULT_PLAN
        for name,value in (("plan.json",self.plan),("selection.private.json",self.selection),("public-schema.json",e.PUBLIC_SCHEMA)):
            write(self.plan_dir+"/"+name,e.canonical(value))

    def close(self):
        for p in reversed(self.patches):p.stop()


class NativeMetadataAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory()
        cls.fixture=NativeFixture(cls.temp.name)

    @classmethod
    def tearDownClass(cls):
        cls.fixture.close();cls.temp.cleanup()

    def test_full_native_annotation_selection_without_any_matrix_or_archive_open(self):
        fixture=self.fixture
        original=Path.open
        def guarded(path,*args,**kwargs):
            if path.name=="matrix.txt" or path.name.endswith(".tar.gz"):
                raise AssertionError("matrix/archive read in metadata-only audit")
            return original(path,*args,**kwargs)
        with patch.object(Path,"open",guarded):
            plan,selection,receipt=e.metadata_only(fixture.plan_dir,e.digest(fixture.plan),fixture.root)
        self.assertEqual(len(selection["retained_keys"]),55460)
        self.assertEqual(sum(len(d["selected"]) for d in selection["donors"]),16659)
        self.assertEqual(receipt["real_matrices_opened"],0)
        self.assertFalse(receipt["expression_authorized"])

    def test_full_native_streams_independently_derive_all_eight_B_R(self):
        fixture=self.fixture
        actual={}
        for pin in fixture.plan["source_manifest"]["matrix_expected_pins_NOT_read"]:
            donor=next(d for d in fixture.selection["donors"] if d["sample"]==pin["sample"])
            with (fixture.root/pin["path"]).open("rb") as handle:
                actual[pin["sample"]]=e.stream_matrix(handle,pin,[c["column_1based"] for c in donor["selected"]],fixture.plan["RNA_rows_1based"])
        self.assertEqual(actual,fixture.expected)
        private,public=payloads(actual,fixture.selection,fixture.plan)
        e.compare_payloads(private,public,actual,fixture.selection,fixture.plan,"4"*64)

    def test_metadata_allowlist_refuses_disguised_numeric_source_before_read(self):
        plan=copy.deepcopy(self.fixture.plan)
        plan["source_manifest"]["annotation_and_receipt_pins"][0]["path"]="work/psc-bach2-cohort-qualification/sample01/members/matrix.txt"
        with self.assertRaises(e.AuditError):e.annotation_review(plan,self.fixture.selection,v,self.fixture.root)

    def test_changed_complete_axis_membership_and_selection_hash(self):
        plan=copy.deepcopy(self.fixture.plan);plan["RNA_rows_1based"][0]=36639
        with self.assertRaises(e.AuditError):e.annotation_review(plan,self.fixture.selection,v,self.fixture.root)
        selection=copy.deepcopy(self.fixture.selection);selection["donors"][0]["selected"][0]["column_1based"]=True
        with self.assertRaises(e.AuditError):e.annotation_review(self.fixture.plan,selection,v,self.fixture.root)


class GuardAudit(unittest.TestCase):
    def test_default_and_missing_acknowledgment_never_open_numeric_files(self):
        with patch.object(e,"metadata_only",return_value=({}, {}, {"status":"synthetic_metadata"})) as metadata, \
             patch.object(e,"audit_real",side_effect=AssertionError("real access")):
            with patch("sys.stdout",new_callable=io.StringIO):self.assertEqual(e.main(["--plan-sha256","1"*64]),0)
            metadata.assert_called_once()
        with patch.object(Path,"open",side_effect=AssertionError("any read")):
            with self.assertRaises(e.AuditError):e.audit_real()
        with patch.object(e,"metadata_only",side_effect=AssertionError("metadata read")),patch("sys.stderr",new_callable=io.StringIO):
            self.assertEqual(e.main(["--authorization","matrix.txt"]),2)

    def test_root_roles_refuse_matrix_traversal_and_oversized_before_parsing(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for path in ("work/psc-bach2-cohort-qualification/sample01/members/matrix.txt","work/psc-bach2-naive-rna-freeze/../matrix.json",
                         "work//psc-bach2-naive-rna-freeze/x.json","/tmp/auth.json"):
                with patch.object(Path,"open",side_effect=AssertionError("wrong role open")),self.assertRaises(e.AuditError):
                    e.root_json(root,path,"1"*64,8192)
            p=root/"work/psc-bach2-naive-rna-freeze/test.json";p.parent.mkdir(parents=True);p.write_bytes(b"x"*9000)
            with patch.object(e,"decode_json",side_effect=AssertionError("oversized parse")),self.assertRaises(e.AuditError):
                e.root_json(root,str(p.relative_to(root)),e.sha(p.read_bytes()),8192)

    def test_all_external_pins_required_before_read(self):
        with patch.object(Path,"open",side_effect=AssertionError("read")):
            with self.assertRaises(e.AuditError):e.audit_real(acknowledge=True)

    def test_symlinks_hardlinks_aliases_traversal_reused_or_source_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/".gitignore").write_text("work/\n")
            (root/"body").write_text("synthetic");os.link(root/"body",root/"hardlink")
            with self.assertRaises(e.AuditError):e.path_at(root,"hardlink",regular=True)
            (root/"link").symlink_to(root/"body")
            with self.assertRaises(e.AuditError):e.path_at(root,"link",regular=True)
            for name in ("a/../body","./body","a//body","/tmp/body","a\\body"):
                with self.subTest(name=name),self.assertRaises(e.AuditError):e.path_at(root,name)
            for out in ("data/raw/output",e.OWN+"/../source",e.OWN+"/a/b"):
                with self.subTest(out=out),self.assertRaises(e.AuditError):e.fresh_output(root,out)
            out=e.OWN+"/once";e.fresh_output(root,out,create=True)
            with self.assertRaises(e.AuditError):e.fresh_output(root,out)

    def test_duplicate_JSON_keys_noncanonical_and_nested_numeric_objects(self):
        for body in (b'{"x":1,"x":2}',b'{"x":NaN}',b'{}'):
            with self.subTest(body=body),self.assertRaises(e.AuditError):e.decode_json(body)
        for value in ({"private":"SYNTHETIC"}, True, "NaN", "1e10001"):
            with self.subTest(value=value),self.assertRaises(e.AuditError):e.number(value)


class CompletedRunAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory();cls.fixture=NativeFixture(cls.temp.name)
        f=cls.fixture
        cls.run_dir=a.BASE+"/runs/SYNTHETIC_PRIMARY"
        cls.auth_path="work/psc-bach2-naive-rna-freeze/SYNTHETIC_authorization.json"
        cls.auth=e.author_authorization(f.plan,cls.run_dir)
        cls.auth_sha=e.digest(cls.auth)
        f.write(cls.auth_path,e.canonical(cls.auth))
        # Genuine native synthetic producer execute; no numerical auditor helper.
        a.execute(f.plan_dir,cls.auth_path,cls.auth_sha,cls.run_dir)
        cls.completion_sha=e.sha((f.root/cls.run_dir/"completion.json").read_bytes())

    @classmethod
    def tearDownClass(cls):
        cls.fixture.close();cls.temp.cleanup()

    def bound_args(self,name):
        f=self.fixture
        output=e.OWN+"/SYNTHETIC_"+name
        seal=e.seal_bindings(f.plan,f.plan_dir,self.run_dir,output,self.auth_path,self.auth_sha,self.completion_sha,f.root)
        seal_path="work/psc-bach2-naive-rna-freeze/SYNTHETIC_seal_"+name+".json"
        f.write(seal_path,e.canonical(seal))
        return dict(acknowledge=True,plan_directory=f.plan_dir,plan_sha=e.digest(f.plan),run_directory=self.run_dir,
                    authorization_path=self.auth_path,authorization_sha=self.auth_sha,seal_path=seal_path,
                    seal_sha=e.digest(seal),completion_sha=self.completion_sha,output=output,root=f.root)

    def test_native_four_artifact_completed_workflow_all_eight_and_aggregate_only_receipt(self):
        args=self.bound_args("complete")
        result=e.audit_real(**args)
        self.assertEqual(result["comparisons"]["donor_records_checked"],8)
        text=e.canonical(result).decode()
        for forbidden in ("SYNTHETIC_DONOR", "sample01", "B_source_units", "RNA_source_units", "omitted_sample", "CI95", "SNP_minus_noSNP"):
            self.assertNotIn(forbidden,text)
        with self.assertRaises(e.AuditError):e.audit_real(**args)

    def test_completion_authentication_does_not_parse_numeric_JSON(self):
        f=self.fixture
        original=e.decode_json
        calls=[]
        def guard(body):
            calls.append(len(body));value=original(body)
            self.assertNotIn("donors",value);self.assertNotIn("contrast",value)
            return value
        with patch.object(e,"decode_json",side_effect=guard):
            payloads=e.authenticate_completed(f.plan,self.run_dir,self.completion_sha,self.auth_sha,f.root)
        self.assertEqual(len(calls),1);self.assertEqual(set(payloads),e.PAYLOADS)

    def test_missing_extra_partial_pending_inventory_refused_before_numeric_decode(self):
        f=self.fixture
        for name in ("failure.json","completion.pending","extra.private.json"):
            p=f.root/self.run_dir/name;p.write_text("SYNTHETIC")
            try:
                with patch.object(e,"decode_json",side_effect=AssertionError("partial decode")),self.assertRaises(e.AuditError):
                    e.authenticate_completed(f.plan,self.run_dir,self.completion_sha,self.auth_sha,f.root)
            finally:p.unlink()
        p=f.root/self.run_dir/"public-aggregate.json";saved=p.read_bytes();st=p.stat();p.unlink()
        try:
            with self.assertRaises(e.AuditError):e.authenticate_completed(f.plan,self.run_dir,self.completion_sha,self.auth_sha,f.root)
        finally:p.write_bytes(saved);os.utime(p,ns=(st.st_atime_ns,st.st_mtime_ns))

    def test_all_payload_hashes_pass_before_any_numerical_parse(self):
        f=self.fixture;p=f.root/self.run_dir/"measurements.private.json";saved=p.read_bytes();st=p.stat()
        p.write_bytes(saved+b" ")
        try:
            with self.assertRaises(e.AuditError):e.authenticate_completed(f.plan,self.run_dir,self.completion_sha,self.auth_sha,f.root)
        finally:p.write_bytes(saved);os.utime(p,ns=(st.st_atime_ns,st.st_mtime_ns))

    def test_stale_completion_and_changed_source_implementation_auth_seal_protocol(self):
        f=self.fixture
        with self.assertRaises(e.AuditError):e.authenticate_completed(f.plan,self.run_dir,"0"*64,self.auth_sha,f.root)
        for field in ("protocol_sha256","runtime_sha256","selection_sha256","adapter_sha256","implementation_sha256","authorization_sha256","completion_sha256"):
            args=self.bound_args("bad_"+field)
            seal=e.json_pin(f.root,args["seal_path"]);seal[field]="0"*64
            f.write(args["seal_path"],e.canonical(seal));args["seal_sha"]=e.digest(seal)
            with patch.object(e,"stream_matrix",side_effect=AssertionError("unbound matrix read")),self.assertRaises(e.AuditError):e.audit_real(**args)

    def test_completion_not_last_written_and_inventory_symlink(self):
        f=self.fixture;p=f.root/self.run_dir/"public-aggregate.json";st=p.stat()
        os.utime(p,ns=(st.st_atime_ns,(f.root/self.run_dir/"completion.json").stat().st_mtime_ns+1_000_000))
        try:
            with self.assertRaises(e.AuditError):e.authenticate_completed(f.plan,self.run_dir,self.completion_sha,self.auth_sha,f.root)
        finally:os.utime(p,ns=(st.st_atime_ns,st.st_mtime_ns))
        saved=p.read_bytes();p.unlink();p.symlink_to(f.root/self.run_dir/"measurements.private.json")
        try:
            with self.assertRaises(e.AuditError):e.authenticate_completed(f.plan,self.run_dir,self.completion_sha,self.auth_sha,f.root)
        finally:p.unlink();p.write_bytes(saved);os.utime(p,ns=(st.st_atime_ns,st.st_mtime_ns))


if __name__ == "__main__":
    unittest.main()
