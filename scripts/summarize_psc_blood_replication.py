#!/usr/bin/env python3
"""Direction-aware intersection-union summaries for PSC blood replication.

This module does not infer disease causality or equate cross-platform effect scales.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def adjust_pvalues(values, method="BH"):
    """Correct the declared family; unavailable tests get 1 only for correction."""
    p = np.asarray(values, dtype=float)
    if p.ndim != 1:
        raise ValueError("P values must be one-dimensional")
    if np.any(np.isinf(p)) or np.any((p[np.isfinite(p)] < 0) | (p[np.isfinite(p)] > 1)):
        raise ValueError("P value outside [0,1]")
    if method not in {"BH", "BY"}:
        raise ValueError("Unknown correction")
    n = len(p)
    if not n:
        return p.copy()
    p = np.where(np.isnan(p), 1.0, p)
    order = np.argsort(p, kind="stable")
    factor = np.sum(1.0 / np.arange(1, n + 1)) if method == "BY" else 1.0
    q = np.minimum.accumulate((p[order] * n * factor / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


def conjunction(pvalues, effects):
    """Max two-sided P with every nonzero effect in the same direction.

    Missing source evidence is separately labeled, not declared a measured null.
    The conservative conjunction P is 1 for unavailable or discordant evidence.
    """
    p = np.asarray(pvalues, dtype=float)
    b = np.asarray(effects, dtype=float)
    if p.ndim != 2 or p.shape != b.shape or p.shape[1] < 2:
        raise ValueError("Need matching gene-by-study arrays for >=2 tests")
    finite_p = p[np.isfinite(p)]
    if np.any(np.isinf(p)) or np.any((finite_p < 0) | (finite_p > 1)):
        raise ValueError("Invalid P value")
    available = np.isfinite(p).all(axis=1) & np.isfinite(b).all(axis=1)
    concordant = available & ((b > 0).all(axis=1) | (b < 0).all(axis=1))
    result = np.ones(p.shape[0])
    result[concordant] = p[concordant].max(axis=1)
    status = np.where(~available, "unavailable_constituent", np.where(concordant, "concordant", "discordant_or_zero_direction"))
    return result, status


def reciprocal_crosswalk(frame):
    """Recompute identity from full relationships, never from a test-selected set."""
    required = {"ensembl_gene_id", "entrez_gene_id", "hgnc_symbol"}
    if not required.issubset(frame.columns):
        raise ValueError("Missing crosswalk fields")
    e_to_t = defaultdict(set)
    t_to_e = defaultdict(set)
    symbols = defaultdict(set)
    ens_universe = set()
    for row in frame.fillna("").to_dict("records"):
        e, t, s = (str(row[k]).strip() for k in ["ensembl_gene_id", "entrez_gene_id", "hgnc_symbol"])
        if e:
            ens_universe.add(e)
        if e and s:
            symbols[e].add(s)
        if e and t:
            if not (e.startswith("ENSG") and e[4:].isdigit() and t.isdigit()):
                raise ValueError("Malformed identity in prior full crosswalk")
            e_to_t[e].add(t)
            t_to_e[t].add(e)
    rows = []
    for e in sorted(ens_universe):
        ts = e_to_t[e]
        t = next(iter(ts)) if len(ts) == 1 else ""
        if not ts:
            status = "no_entrez_relation"
        elif len(ts) != 1:
            status = "ambiguous_ensembl_to_entrez"
        elif len(t_to_e[t]) != 1:
            status = "ambiguous_entrez_to_ensembl"
        else:
            status = "reciprocal_one_to_one"
        rows.append({"ensembl_gene_id": e, "entrez_gene_id": t if status == "reciprocal_one_to_one" else "", "hgnc_symbol": next(iter(symbols[e])) if len(symbols[e]) == 1 else "", "source_symbols": ";".join(sorted(symbols[e])), "identity_status": status})
    return pd.DataFrame(rows)


def build_conjunctions(frame, constituents):
    out = frame.copy()
    for family, columns in constituents.items():
        p = out[[x[0] for x in columns]].to_numpy(dtype=float)
        effects = out[[x[1] for x in columns]].to_numpy(dtype=float)
        value, status = conjunction(p, effects)
        out[f"{family}_pvalue"] = value
        out[f"{family}_status"] = status
        out[f"{family}_qvalue_bh"] = adjust_pvalues(value, "BH")
        out[f"{family}_qvalue_by"] = adjust_pvalues(value, "BY")
    return out



def assemble_results(rna, array, crosswalk, schema):
    """Join all eligible shared genes; never define family by significant/nonmissing P."""
    rna = rna.copy()
    array = array.copy()
    gene_col = schema.get("rna_gene_id", "gene_id")
    filter_col = schema.get("rna_filter", "count_filter_pass")
    rna[gene_col] = rna[gene_col].astype(str)
    if rna[gene_col].duplicated().any():
        raise ValueError("Duplicate RNA stable gene IDs")
    if filter_col not in rna:
        raise ValueError("No explicit RNA count-filter accounting")
    flags = rna[filter_col].astype(str).str.lower()
    if not flags.isin(["true", "false"]).all():
        raise ValueError("Invalid count-filter flags")
    mapping = reciprocal_crosswalk(crosswalk)
    audit = rna[[gene_col, filter_col]].merge(mapping, left_on=gene_col, right_on="ensembl_gene_id", how="left", validate="one_to_one")
    audit["identity_status"] = audit["identity_status"].fillna("not_in_mapping_snapshot")
    array_id = schema.get("array_gene_id", "entrez_id")
    array_contrast = schema.get("array_contrast", "contrast")
    array[array_id] = array[array_id].astype(str)
    if array.duplicated([array_id, array_contrast]).any():
        raise ValueError("Duplicate array gene/contrast row")
    expected_contrasts = schema["array_contrasts"]
    if set(array[array_contrast]) != set(expected_contrasts.values()):
        raise ValueError("Unexpected/incomplete array contrasts")
    sets = {key: set(array.loc[array[array_contrast].eq(value), array_id]) for key, value in expected_contrasts.items()}
    if any(value != next(iter(sets.values())) for value in sets.values()):
        raise ValueError("Array family changed across contrasts")
    array_genes = next(iter(sets.values()))
    audit["array_gene_measured"] = audit["entrez_gene_id"].isin(array_genes)
    audit["eligible_shared_family"] = flags.eq("true").to_numpy() & audit["identity_status"].eq("reciprocal_one_to_one") & audit["array_gene_measured"]
    audit["family_status"] = np.select([flags.eq("false").to_numpy(), ~audit["identity_status"].eq("reciprocal_one_to_one"), ~audit["array_gene_measured"]], ["rna_low_count", "no_unambiguous_cross_platform_identity", "not_measured_on_array"], default="eligible_shared_family")
    selected = audit.loc[audit["eligible_shared_family"], ["ensembl_gene_id", "entrez_gene_id", "hgnc_symbol", "source_symbols"]]
    # Sorting fixes public output order independently of source row or dataframe join order.
    selected = selected.sort_values("ensembl_gene_id", kind="stable").reset_index(drop=True)
    rcols = {schema["rna_pvalue"]: "rna_pvalue", schema["rna_effect"]: "rna_effect"}
    rcols.update(schema.get("rna_extra_columns", {}))
    rr = rna[[gene_col, *rcols]].rename(columns={gene_col: "ensembl_gene_id", **rcols})
    selected = selected.merge(rr, on="ensembl_gene_id", how="left", validate="one_to_one")
    constituents = {}
    for short, contrast in expected_contrasts.items():
        acolumns = {schema["array_pvalue"]: f"array_{short}_pvalue", schema["array_effect"]: f"array_{short}_effect"}
        acolumns.update({source: f"array_{short}_{target}" for source,target in schema.get("array_extra_columns", {}).items()})
        aa = array.loc[array[array_contrast].eq(contrast), [array_id,*acolumns]].rename(columns={array_id: "entrez_gene_id", **acolumns})
        selected = selected.merge(aa, on="entrez_gene_id", how="left", validate="one_to_one")
    core = [("rna_pvalue", "rna_effect"), ("array_control_pvalue", "array_control_effect")]
    constituents["replication"] = core
    constituents["uc_support"] = core + [("array_uc_pvalue", "array_uc_effect")]
    constituents["all_disease_support"] = constituents["uc_support"] + [("array_pbc_pvalue", "array_pbc_effect"), ("array_cd_pvalue", "array_cd_effect")]
    return build_conjunctions(selected, constituents), audit


def run(plan_path, output_dir):
    plan_path = Path(plan_path)
    plan = json.loads(plan_path.read_text())
    if plan.get("status") != "frozen" or not plan.get("allow_integration"):
        raise ValueError("A frozen integration plan is required")
    root = Path(__file__).resolve().parents[1]
    for rel, expected in plan["input_sha256"].items():
        if file_sha256(root / rel) != expected:
            raise ValueError(f"Input hash mismatch: {rel}")
    rna = pd.read_csv(root / plan["rna_results"], dtype={plan["schema"]["rna_gene_id"]: str})
    array = pd.read_csv(root / plan["array_results"], dtype={plan["schema"]["array_gene_id"]: str})
    crosswalk = pd.read_csv(root / plan["crosswalk"], dtype=str, keep_default_na=False)
    result, audit = assemble_results(rna, array, crosswalk, plan["schema"])
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for filename, table in [("cross-cohort-gene-results.csv.gz", result), ("cross-platform-coverage.csv.gz", audit)]:
        table.to_csv(output_dir / filename, index=False, float_format="%.17g", compression={"method": "gzip", "mtime": 0})
        outputs[filename] = file_sha256(output_dir / filename)
    summary = {"shared_gene_family": len(result), "rna_coverage_status": audit["family_status"].value_counts().to_dict(), "families": {}, "output_sha256": outputs, "plan_sha256": file_sha256(plan_path)}
    for family in ["replication", "uc_support", "all_disease_support"]:
        passing = result[f"{family}_qvalue_bh"] <= 0.05
        summary["families"][family] = {"q_bh_le_0_05": int(passing.sum()), "q_by_le_0_05": int((result[f"{family}_qvalue_by"] <= 0.05).sum()), "passing_up_in_psc": int((passing & result.rna_effect.gt(0)).sum()), "passing_down_in_psc": int((passing & result.rna_effect.lt(0)).sum()), "status_counts": result[f"{family}_status"].value_counts().to_dict()}
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.plan, args.output_dir), indent=2))


if __name__ == "__main__":
    main()
