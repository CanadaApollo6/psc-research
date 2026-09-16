"""Offline, source-only qualification of one fixed unsigned ECM gene set.

No network, tissue values, scores, models or IL17 fallback are implemented.
The unavailable epithelial IL17A perturbation endpoint remains unavailable.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from html.parser import HTMLParser
import io
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config/psc-liver-external-program-sources.json"
COLUMNS = [
    "program_id", "official_program_name", "msigdb_systematic_id", "exact_source_id",
    "msigdb_version", "reactome_release", "source_species", "source_member_index",
    "source_value", "source_identifier_namespace", "published_ncbi_gene_id",
    "original_source_ids", "original_source_id_namespace", "original_source_row_indices",
    "original_source_id_count", "published_direction", "membership_weight",
    "mapping_authority", "membership_json_sha256", "source_id_table_html_sha256",
]
HEADERS = ["Source Id", "NCBI (Entrez) Gene Id", "Gene Symbol", "Gene Description"]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


class VisibleText(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.parts: list[str] = []

    def handle_starttag(self, tag, attrs):
        if tag == "br":
            self.parts.append(" ")

    def handle_data(self, data):
        self.parts.append(data)


def text(fragment: str) -> str:
    parser = VisibleText()
    parser.feed(fragment)
    parser.close()
    return " ".join("".join(parser.parts).split())


def metadata(html: str, label: str) -> str:
    found = re.findall(r"<th>" + re.escape(label) + r"</th>\s*<td>(.*?)</td>", html, re.S)
    require(len(found) == 1, "Missing/repeated metadata: " + label)
    return text(found[0])


def parse_sources(html: str, document: dict, program: dict) -> tuple[list[dict], dict]:
    """Return the complete official gene-symbol projection, never a filtered subset."""
    name = program["official_program_name"]
    require(list(document) == [name], "Unexpected program(s) in official JSON")
    official = document[name]
    expected_meta = {
        "Standard name": name,
        "Systematic name": program["msigdb_systematic_id"],
        "Exact source": program["exact_source_id"],
        "Source species": program["organism"],
        "Source platform or<br/>identifier namespace": program["source_identifier_namespace"],
        "Version history": f"{program['msigdb_version']}: Updated to Reactome {program['reactome_release']}",
    }
    for label, expected in expected_meta.items():
        require(metadata(html, label) == expected, "Metadata mismatch: " + label)
    for field, key in [("systematicName", "msigdb_systematic_id"), ("exactSource", "exact_source_id"),
                       ("collection", "collection")]:
        require(official[field] == program[key], "Official JSON metadata mismatch: " + field)
    require(f".v{program['msigdb_version']}.json" in program["source_disposition_header"],
            "JSON response disposition does not pin release")
    require(program["direction"] == "UNSIGNED" and program["membership_weight"] == "1",
            "This source has no signed response direction or response weights")
    require(program["published_response_directions_available"] is False,
            "Cannot invent published response directions")
    require(program["scoring_authorized"] is False, "Qualification does not authorize scoring")

    claims = re.findall(r"(\d+) source identifiers mapped to (\d+) genes", html)
    require(len(claims) == 1, "Missing/repeated complete-membership declaration")
    declared_sources, declared_symbols = map(int, claims[0])
    require(declared_sources == program["expected_source_identifier_count"], "Source count pin mismatch")
    require(declared_symbols == program["expected_member_count"], "Member count pin mismatch")
    require(html.count('<div id="geneListing"') == 1, "Missing/repeated gene listing")
    start = html.index('<div id="geneListing"')
    require("</table>" in html[start:], "Incomplete gene listing")
    table = html[start:html.index("</table>", start) + len("</table>")]
    rows = re.findall(r"<tr\b[^>]*>(.*?)</tr>", table, re.S)
    require(bool(rows), "Missing table rows")
    headers = [text(value) for value in re.findall(r"<th\b[^>]*>(.*?)</th>", rows[0], re.S)]
    require(headers == HEADERS, "Unexpected gene-list schema")
    source_rows = []
    for index, row in enumerate(rows[1:], 1):
        cells = [text(value) for value in re.findall(r"<td\b[^>]*>(.*?)</td>", row, re.S)]
        require(len(cells) == 4, f"Incomplete source row {index}")
        source_id, ncbi, symbol, description = cells
        require(bool(re.fullmatch(r"ENSG\d{11}", source_id)), f"Invalid/missing source ID at {index}")
        require(bool(re.fullmatch(r"[1-9]\d*", ncbi)), f"Invalid/missing NCBI ID at {index}")
        require(bool(symbol) and not re.search(r"\s|;", symbol), f"Invalid/missing symbol at {index}")
        source_rows.append(dict(index=index, source_id=source_id, ncbi=ncbi, symbol=symbol))
    require(len(source_rows) == declared_sources, "Source row count mismatch")
    require(len({row["source_id"] for row in source_rows}) == declared_sources,
            "Repeated original source IDs must not be silently collapsed")

    symbols = official["geneSymbols"]
    require(isinstance(symbols, list) and all(isinstance(v, str) and v for v in symbols),
            "Missing/invalid official geneSymbols list")
    require(len(symbols) == len(set(symbols)) == declared_symbols, "Repeated/missing official symbols")
    require(set(symbols) == {row["symbol"] for row in source_rows}, "Official symbol set differs from HTML")
    groups = {symbol: [row for row in source_rows if row["symbol"] == symbol] for symbol in symbols}
    ncbi_to_symbol = {}
    projection = []
    for index, (symbol, group) in enumerate(groups.items(), 1):
        ids = {row["ncbi"] for row in group}
        require(len(ids) == 1, "Conflicting symbol-to-NCBI source mapping: " + symbol)
        ncbi = next(iter(ids))
        require(ncbi not in ncbi_to_symbol, "Conflicting NCBI-to-symbol source mapping: " + ncbi)
        ncbi_to_symbol[ncbi] = symbol
        projection.append(dict(index=index, symbol=symbol, ncbi=ncbi,
                               source_ids=[row["source_id"] for row in group],
                               source_indices=[row["index"] for row in group]))
    counts = {"source_identifier_rows": len(source_rows), "official_member_symbols": len(projection),
              "source_ids_preserved_once": sum(len(row["source_ids"]) for row in projection),
              "many_to_one_symbol_groups": sum(len(group) > 1 for group in groups.values()),
              "extra_source_ids_beyond_official_symbols": len(source_rows) - len(projection),
              "missing_source_ids": 0, "missing_symbols": 0, "missing_ncbi_ids": 0,
              "source_mapping_conflicts": 0, "direction": "UNSIGNED", "membership_weight": "1"}
    return projection, counts


def render_membership(projection: list[dict], program: dict, sources: dict) -> bytes:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=COLUMNS, lineterminator="\n")
    writer.writeheader()
    for row in projection:
        writer.writerow({
            "program_id": program["program_id"], "official_program_name": program["official_program_name"],
            "msigdb_systematic_id": program["msigdb_systematic_id"], "exact_source_id": program["exact_source_id"],
            "msigdb_version": program["msigdb_version"], "reactome_release": program["reactome_release"],
            "source_species": program["organism"], "source_member_index": row["index"],
            "source_value": row["symbol"], "source_identifier_namespace": program["member_namespace"],
            "published_ncbi_gene_id": row["ncbi"], "original_source_ids": ";".join(row["source_ids"]),
            "original_source_id_namespace": program["source_identifier_namespace"],
            "original_source_row_indices": ";".join(map(str, row["source_indices"])),
            "original_source_id_count": len(row["source_ids"]), "published_direction": "UNSIGNED",
            "membership_weight": "1", "mapping_authority": "MSigDB_published_projection_not_independently_reconciled",
            "membership_json_sha256": sources[program["source_keys"]["json"]]["sha256"],
            "source_id_table_html_sha256": sources[program["source_keys"]["html"]]["sha256"],
        })
    return output.getvalue().encode("utf-8")


def pinned_bytes(root: Path, pin: dict) -> bytes:
    relative = Path(pin["path"])
    require(not relative.is_absolute() and ".." not in relative.parts, "Source path leaves repository")
    resolved = (root / relative).resolve()
    require(resolved.is_relative_to(root.resolve()), "Source symlink leaves repository")
    payload = resolved.read_bytes()
    require(len(payload) == pin["bytes"], "Byte-count mismatch: " + pin["path"])
    require(hashlib.sha256(payload).hexdigest() == pin["sha256"], "SHA-256 mismatch: " + pin["path"])
    return payload


def qualify(config: dict, root: Path) -> tuple[bytes, dict]:
    require(config["source_pass_closed"] is True, "Source pass must be closed")
    require([p["program_id"] for p in config["programs"]] == ["reactome_ecm_organization"],
            "This bounded qualifier has no alternative endpoint")
    absent = config["unavailable_requested_programs"]
    require(len(absent) == 1 and absent[0]["program_id"] == "epithelial_il17a_alone_response",
            "Missing IL17 HOLD record")
    require(absent[0]["status"] == "HOLD_no_complete_qualified_response_source", "IL17 endpoint is not qualified")
    require(absent[0]["member_count"] is None and absent[0]["scoring_authorized"] is False,
            "Unavailable IL17 response cannot become an empty list or score")
    sources = config["sources"]
    require(set(sources) == {"ecm_html", "ecm_json", "nograles_bibliographic_search"}, "Unexpected source inputs")
    for pin in sources.values():
        require(pin["path"].startswith("data/raw/psc-liver-external-programs/"), "Input outside owned source area")
    for pin in config.get("ignored_audit_pins", []):
        require(pin["path"].startswith(("data/raw/psc-liver-external-programs/", "work/psc-liver-external-programs/")),
                "Audit input outside owned source area")
    payloads = {key: pinned_bytes(root, pin) for key, pin in sources.items()}
    for pin in config.get("ignored_audit_pins", []):
        pinned_bytes(root, pin)
    program = config["programs"][0]
    hkey, jkey = program["source_keys"]["html"], program["source_keys"]["json"]
    for key in [hkey, jkey]:
        require(sources[key]["status"] == 200 and sources[key]["complete"] is True,
                "Incomplete membership capture: " + key)
    projection, counts = parse_sources(payloads[hkey].decode("iso-8859-1"),
                                        json.loads(payloads[jkey].decode("utf-8")), program)
    for key in ["many_to_one_symbol_groups", "extra_source_ids_beyond_official_symbols"]:
        require(counts[key] == program[key], "Projection count pin mismatch: " + key)
    require(counts["official_member_symbols"] == config["membership_output"]["rows"], "Output row pin mismatch")
    output = render_membership(projection, program, sources)
    output_pin = config["membership_output"]
    if output_pin["sha256"] is not None:
        require(hashlib.sha256(output).hexdigest() == output_pin["sha256"], "Rebuilt membership hash mismatch")
        require(len(output) == output_pin["bytes"], "Rebuilt membership byte-count mismatch")
    result = {"status": "source_only_ECM_qualified_IL17_response_HOLD", "counts": counts,
              "membership_sha256": hashlib.sha256(output).hexdigest(), "membership_bytes": len(output),
              "source_pins_checked": len(sources), "ignored_audit_pins_checked": len(config.get("ignored_audit_pins", [])),
              "network_requests": 0, "liver_or_organoid_expression_reads": 0,
              "no_normalization_scoring_or_hypothesis_tests": True}
    return output, result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--build", action="store_true", help="Rebuild only the owned membership CSV offline")
    group.add_argument("--verify", action="store_true", help="Read-only pin and exact CSV replay checks")
    args = parser.parse_args()
    config = json.loads(CONFIG.read_text())
    output, result = qualify(config, ROOT)
    destination = ROOT / config["membership_output"]["path"]
    require(destination == ROOT / "data/derived/psc-liver-external-program-membership.csv", "Unexpected output path")
    if args.build:
        destination.write_bytes(output)
        result["action"] = "built_only_owned_membership_CSV"
    else:
        require(config["membership_output"]["sha256"] is not None, "Output must be pinned for verification")
        require(destination.read_bytes() == output, "Existing CSV differs from exact offline reconstruction")
        result["action"] = "verified_read_only_exact_membership_replay"
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
