#!/usr/bin/env python3
"""Aggregate diagnostics and benchmark reporting; no new hypothesis tests."""
from pathlib import Path
import gzip
import hashlib
import json
import re
from collections import Counter
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = {"UBASH3A":"ENSG00000160185", "ETS2":"ENSG00000157557", "PRKD2":"ENSG00000105287", "PFKFB3":"ENSG00000170525", "BCL2L11":"ENSG00000153094", "IFIT1":"ENSG00000185745", "G0S2":"ENSG00000123689"}


def main():
    out = ROOT / "data/derived/psc-blood-replication"
    result = pd.read_csv(out / "cross-cohort-gene-results.csv.gz")
    rna = pd.read_csv(ROOT / "data/derived/psc-blood-rnaseq/rnaseq-primary.csv")
    leave = pd.read_csv(ROOT / "data/derived/psc-blood-rnaseq/rnaseq-leave-one-plate-out.csv")
    qc = json.loads((ROOT / "data/derived/psc-blood-rnaseq/rnaseq-qc.json").read_text())
    original_map = pd.read_csv(ROOT / "data/derived/gse84161-gene-map.csv.gz", dtype=str, keep_default_na=False)
    for symbol, gene in BENCHMARKS.items():
        ids = set(original_map.loc[original_map.hgnc_symbol.eq(symbol), "ensembl_gene_id"])
        if ids != {gene}:
            raise ValueError(f"Benchmark source identity mismatch: {symbol}")
    joined = result.merge(leave[["gene_id", "n_estimable_omissions", "n_same_nonzero_sign", "n_opposite_sign"]], left_on="ensembl_gene_id", right_on="gene_id", how="left", validate="one_to_one")
    if joined.n_estimable_omissions.isna().any():
        raise ValueError("Missing prespecified leave-plate diagnostics")
    flagcols = ["rna_LFC_converged", "rna_MAP_converged", "rna_genewise_converged"]
    joined["any_reported_nb_nonconvergence"] = (~joined[flagcols].eq(True)).any(axis=1)
    joined["rna_ols_matches_nb_direction"] = joined.rna_effect.mul(joined.rna_ols_effect).gt(0)
    joined["psc_alone_ols_matches_nb_direction"] = joined.rna_effect.mul(joined.rna_psc_alone_ols_effect).gt(0)
    joined["pscuc_ols_matches_nb_direction"] = joined.rna_effect.mul(joined.rna_pscuc_ols_effect).gt(0)
    hits = joined.loc[joined.replication_qvalue_bh.le(0.05)].sort_values(["replication_pvalue", "ensembl_gene_id"])
    bench = joined.loc[joined.ensembl_gene_id.isin(BENCHMARKS.values())].copy()
    if set(bench.ensembl_gene_id) != set(BENCHMARKS.values()):
        raise ValueError("Benchmark is outside shared family; report coverage before selection")
    bench["prespecified_benchmark"] = bench.ensembl_gene_id.map({v:k for k,v in BENCHMARKS.items()})
    for name, frame in [("bh-candidates-with-diagnostics.csv.gz", hits), ("prespecified-benchmarks.csv.gz", bench)]:
        frame.to_csv(out / name, index=False, float_format="%.17g", compression={"method":"gzip", "mtime":0})
    warnings = Counter()
    log = (ROOT / "work/research-path/rnaseq-run.log").read_text()
    for line in log.splitlines():
        m = re.search(r"(\w+Warning): (.+)$", line)
        if m:
            warnings[m.group(0)] += 1
    outlier_cutoff = qc["negative_binomial"]["cooks_distance_cutoff"]
    all_cooks = rna.loc[rna.gene_id.isin(hits.ensembl_gene_id), "max_cooks_distance_all_samples"]
    flag_hits = hits.loc[hits.any_reported_nb_nonconvergence]
    diag = {
        "scope":"Descriptive accounting of frozen results/sensitivities; no new significance family or amended primary analysis.",
        "shared_family":len(joined), "bh_candidates":len(hits),
        "rna_primary_bh_q_le_0_05":int((rna.eligible & rna.padj.le(0.05)).sum()),
        "rna_ols_primary_bh_q_le_0_05":int((rna.eligible & rna.ols_padj_bh.le(0.05)).sum()),
        "all_shared_fit_flags":{c:int(joined[c].eq(False).sum()) for c in flagcols},
        "bh_candidate_fit_flags":{c:int(hits[c].eq(False).sum()) for c in flagcols},
        "bh_candidate_any_nonconvergence":int(hits.any_reported_nb_nonconvergence.sum()),
        "flagged_bh_candidates":flag_hits[["ensembl_gene_id","hgnc_symbol",*flagcols]].to_dict("records"),
        "bh_candidate_ols_direction_match":int(hits.rna_ols_matches_nb_direction.sum()),
        "bh_candidate_ols_bh_q_le_0_05_and_direction_match":int((hits.rna_ols_qvalue_bh.le(0.05)&hits.rna_ols_matches_nb_direction).sum()),
        "bh_candidate_all_five_plate_omissions_match_full_ols_direction":int((hits.n_estimable_omissions.eq(5)&hits.n_same_nonzero_sign.eq(5)).sum()),
        "bh_candidate_both_subgroup_ols_directions_match_nb":int((hits.psc_alone_ols_matches_nb_direction&hits.pscuc_ols_matches_nb_direction).sum()),
        "bh_candidate_no_fit_flags_and_both_ols_and_five_omission_sign_checks":int((~hits.any_reported_nb_nonconvergence & hits.rna_ols_matches_nb_direction & hits.n_same_nonzero_sign.eq(5) & hits.n_estimable_omissions.eq(5)).sum()),
        "bh_candidate_max_cooks_any_sample_above_package_cutoff":int(all_cooks.gt(outlier_cutoff).sum()),
        "cooks_caveat":"Package filters only 34 samples with >=3 identical full design rows including continuous age; all-sample Cook values are diagnostics, not retrospectively imposed exclusions.",
        "lowest_primary_bh_q":float(joined.replication_qvalue_bh.min()),
        "lowest_primary_by_q":float(joined.replication_qvalue_by.min()),
        "lowest_uc_support_bh_q":float(joined.uc_support_qvalue_bh.min()),
        "lowest_all_disease_support_bh_q":float(joined.all_disease_support_qvalue_bh.min()),
        "benchmark_source_identity":BENCHMARKS,
        "benchmark_role":"Already named project/source benchmarks, not discoveries. All seven have unambiguous reciprocal IDs and pass measured-family eligibility.",
        "worker_runtime_warning_counts":dict(sorted(warnings.items())),
        "warnings_caveat":"Joblib worker warnings appear in the outer stderr log; the RNA QC main-process captured-warning list is empty and does not imply a warning-free fit.",
    }
    (out / "diagnostics.json").write_text(json.dumps(diag,indent=2,allow_nan=False)+"\n")
    print(json.dumps(diag,indent=2,allow_nan=False))


if __name__ == "__main__":
    main()
