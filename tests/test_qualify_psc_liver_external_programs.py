"""Synthetic source/parser guards; no tissue data or network required."""
import copy
import csv
import hashlib
import io
from pathlib import Path
import tempfile
import unittest

from scripts import qualify_psc_liver_external_programs as q


def fixture():
    program = {
        "program_id": "reactome_ecm_organization", "official_program_name": "TEST_ECM", "msigdb_systematic_id": "M610",
        "exact_source_id": "R-HSA-1474244", "organism": "Homo sapiens", "source_identifier_namespace": "Human_Ensembl_Gene_ID",
        "msigdb_version": "2026.1.Hs", "reactome_release": "95", "collection": "C2:CP:REACTOME",
        "source_disposition_header": 'attachment; filename="TEST_ECM.v2026.1.Hs.json";',
        "direction": "UNSIGNED", "membership_weight": "1", "published_response_directions_available": False,
        "scoring_authorized": False, "expected_source_identifier_count": 3, "expected_member_count": 2,
        "member_namespace": "MSigDB_geneSymbols", "source_keys": {"html": "html", "json": "json"},
    }
    metadata = {
        "Standard name": "TEST_ECM", "Systematic name": "M610", "Exact source": "R-HSA-1474244",
        "Source species": "Homo sapiens", "Source platform or<br/>identifier namespace": "Human_Ensembl_Gene_ID",
        "Version history": "2026.1.Hs: Updated to Reactome 95",
    }
    html = "".join(f"<tr><th>{k}</th><td>{v}</td></tr>" for k, v in metadata.items())
    html += '(show 3 source identifiers mapped to 2 genes)<div id="geneListing"><table><tr>'
    html += "".join("<th>" + value + "</th>" for value in q.HEADERS) + "</tr>"
    for index, ncbi, symbol in [(1, "1", "AAA"), (2, "1", "AAA"), (3, "2", "BBB")]:
        html += f"<tr><td>ENSG{index:011}</td><td>{ncbi}</td><td>{symbol}</td><td>synthetic</td></tr>"
    html += "</table></div>"
    doc = {"TEST_ECM": {"systematicName": "M610", "exactSource": "R-HSA-1474244", "collection": "C2:CP:REACTOME", "geneSymbols": ["BBB", "AAA"]}}
    return html, doc, program


class TestExternalProgramParser(unittest.TestCase):
    def setUp(self):
        self.html, self.doc, self.program = fixture()

    def parse(self):
        return q.parse_sources(self.html, self.doc, self.program)

    def test_complete_projection_retains_all_original_ids_and_json_order(self):
        rows, counts = self.parse()
        self.assertEqual([r["symbol"] for r in rows], ["BBB", "AAA"])
        self.assertEqual(rows[1]["source_ids"], ["ENSG00000000001", "ENSG00000000002"])
        self.assertEqual(rows[1]["source_indices"], [1, 2])
        self.assertEqual(counts["source_ids_preserved_once"], 3)
        self.assertEqual(counts["many_to_one_symbol_groups"], 1)

    def test_no_double_weight_for_many_source_ids(self):
        rows, _ = self.parse()
        payload = q.render_membership(rows, self.program, {"html": {"sha256": "a"*64}, "json": {"sha256": "b"*64}})
        result = list(csv.DictReader(io.StringIO(payload.decode())))
        self.assertEqual(len(result), 2)
        self.assertEqual([r["membership_weight"] for r in result], ["1", "1"])
        self.assertEqual([r["published_direction"] for r in result], ["UNSIGNED", "UNSIGNED"])
        self.assertEqual(result[1]["original_source_ids"], "ENSG00000000001;ENSG00000000002")
        self.assertTrue(payload.endswith(b"\n"))
        self.assertNotIn(b"\r", payload)

    def test_html_entities_and_breaks_only(self):
        self.assertEqual(q.text("NCBI&nbsp;(Entrez)<br/>Gene&nbsp;Id"), "NCBI (Entrez) Gene Id")
        self.assertEqual(q.text("<a>abc-1</a>"), "abc-1")

    def test_missing_original_row_fails(self):
        self.html = self.html.replace("<tr><td>ENSG00000000003</td><td>2</td><td>BBB</td><td>synthetic</td></tr>", "")
        with self.assertRaisesRegex(ValueError, "Source row count"):
            self.parse()

    def test_duplicate_original_id_fails(self):
        self.html = self.html.replace("ENSG00000000002", "ENSG00000000001")
        with self.assertRaisesRegex(ValueError, "Repeated original"):
            self.parse()

    def test_invalid_source_id_fails(self):
        self.html = self.html.replace("ENSG00000000002", "ENSG00000000002.1")
        with self.assertRaisesRegex(ValueError, "Invalid/missing source ID"):
            self.parse()

    def test_missing_symbol_fails(self):
        self.html = self.html.replace("<td>BBB</td>", "<td></td>")
        with self.assertRaisesRegex(ValueError, "Invalid/missing symbol"):
            self.parse()

    def test_missing_ncbi_id_fails(self):
        self.html = self.html.replace("<td>2</td>", "<td></td>")
        with self.assertRaisesRegex(ValueError, "Invalid/missing NCBI"):
            self.parse()

    def test_conflicting_symbol_ncbi_projection_fails(self):
        self.html = self.html.replace("ENSG00000000002</td><td>1</td>", "ENSG00000000002</td><td>3</td>")
        with self.assertRaisesRegex(ValueError, "Conflicting symbol-to-NCBI"):
            self.parse()

    def test_conflicting_ncbi_symbol_projection_fails(self):
        self.html = self.html.replace("ENSG00000000003</td><td>2</td>", "ENSG00000000003</td><td>1</td>")
        with self.assertRaisesRegex(ValueError, "Conflicting NCBI-to-symbol"):
            self.parse()

    def test_missing_official_symbol_fails(self):
        self.doc["TEST_ECM"]["geneSymbols"].pop()
        with self.assertRaisesRegex(ValueError, "Repeated/missing official"):
            self.parse()

    def test_duplicate_official_symbol_fails(self):
        self.doc["TEST_ECM"]["geneSymbols"] = ["AAA", "AAA"]
        with self.assertRaisesRegex(ValueError, "Repeated/missing official"):
            self.parse()

    def test_replaced_official_symbol_fails(self):
        self.doc["TEST_ECM"]["geneSymbols"] = ["AAA", "CCC"]
        with self.assertRaisesRegex(ValueError, "Official symbol set differs"):
            self.parse()

    def test_wrong_declared_count_fails(self):
        self.html = self.html.replace("3 source identifiers", "4 source identifiers")
        with self.assertRaisesRegex(ValueError, "Source count pin"):
            self.parse()

    def test_truncated_listing_fails(self):
        self.html = self.html.replace("</table>", "")
        with self.assertRaisesRegex(ValueError, "Incomplete gene listing"):
            self.parse()

    def test_wrong_release_fails(self):
        self.html = self.html.replace("Updated to Reactome 95", "Updated to Reactome 94")
        with self.assertRaisesRegex(ValueError, "Metadata mismatch"):
            self.parse()

    def test_wrong_download_release_fails(self):
        self.program["source_disposition_header"] = "2025.1.Hs.json"
        with self.assertRaisesRegex(ValueError, "disposition"):
            self.parse()

    def test_no_invented_response_direction(self):
        self.program["direction"] = "UP"
        with self.assertRaisesRegex(ValueError, "no signed response"):
            self.parse()

    def test_no_invented_published_direction_availability(self):
        self.program["published_response_directions_available"] = True
        with self.assertRaisesRegex(ValueError, "invent published"):
            self.parse()

    def test_no_scoring_authorization(self):
        self.program["scoring_authorized"] = True
        with self.assertRaisesRegex(ValueError, "does not authorize scoring"):
            self.parse()

    def test_repeated_metadata_fails(self):
        self.html += "<th>Standard name</th><td>TEST_ECM</td>"
        with self.assertRaisesRegex(ValueError, "Missing/repeated metadata"):
            self.parse()

    def test_wrong_header_fails(self):
        self.html = self.html.replace("<th>Gene Symbol</th>", "<th>Alias Symbol</th>")
        with self.assertRaisesRegex(ValueError, "Unexpected gene-list schema"):
            self.parse()

    def test_extra_program_fails(self):
        self.doc["OTHER"] = copy.deepcopy(self.doc["TEST_ECM"])
        with self.assertRaisesRegex(ValueError, "Unexpected program"):
            self.parse()


class TestExternalProgramPinsAndHold(unittest.TestCase):
    def test_byte_hash_and_relative_path_checks(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            data = b"public synthetic metadata\n"
            (root/"source.txt").write_bytes(data)
            pin = {"path": "source.txt", "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
            self.assertEqual(q.pinned_bytes(root, pin), data)
            for changed in [dict(pin, bytes=len(data)+1), dict(pin, sha256="0"*64), dict(pin, path="../source.txt"), dict(pin, path="/etc/passwd")]:
                with self.assertRaises(ValueError):
                    q.pinned_bytes(root, changed)

    def test_unavailable_il17_is_not_zero_or_score(self):
        config = {"source_pass_closed": True, "programs": [{"program_id": "reactome_ecm_organization"}],
                  "unavailable_requested_programs": [{"program_id": "epithelial_il17a_alone_response",
                    "status": "HOLD_no_complete_qualified_response_source", "member_count": 0, "scoring_authorized": False}]}
        with self.assertRaisesRegex(ValueError, "cannot become an empty list"):
            q.qualify(config, Path("."))
        config["unavailable_requested_programs"][0].update(member_count=None, scoring_authorized=True)
        with self.assertRaisesRegex(ValueError, "cannot become an empty list"):
            q.qualify(config, Path("."))

    def test_no_endpoint_substitution(self):
        config = {"source_pass_closed": True, "programs": [{"program_id": "IL17_pathway_fallback"}]}
        with self.assertRaisesRegex(ValueError, "no alternative endpoint"):
            q.qualify(config, Path("."))

    def test_source_pass_must_be_closed(self):
        with self.assertRaisesRegex(ValueError, "must be closed"):
            q.qualify({"source_pass_closed": False}, Path("."))


if __name__ == "__main__":
    unittest.main()
