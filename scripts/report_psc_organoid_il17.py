#!/usr/bin/env python3
"""Summarize and plot frozen aggregate organoid outputs. No models are refitted."""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ID = "source_feature_id"
TABLES = {
    "primary": "organoid-primary.csv.gz",
    "within": "organoid-within-group.csv.gz",
    "welch": "organoid-welch.csv.gz",
    "leaveout": "organoid-leave-one-donor-out.csv.gz",
    "allocations": "organoid-allocations.csv.gz",
    "universe": "organoid-feature-universe.csv.gz",
}
ANNOTATION = [ID, "source_gene_label", "source_feature_type", "eligible",
              "libraries_count_ge_threshold", "libraries_nonzero",
              "all_zero_released", "mean_raw_count"]


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def number(values: pd.Series) -> pd.Series:
    return pd.to_numeric(values.replace("", np.nan), errors="raise")


def truth(values: pd.Series) -> pd.Series:
    require(set(values).issubset({"True", "False"}), "Invalid boolean literal")
    return values.eq("True")


def validate_tables(tables: dict[str, pd.DataFrame], qc: dict) -> None:
    p = tables["primary"]
    require(not p[ID].duplicated().any(), "Duplicate primary feature ID")
    require(len(p) == qc["released_features"], "Source coverage mismatch")
    require(int(truth(p["eligible"]).sum()) == qc["eligible_features"], "Family mismatch")
    for name, frame in tables.items():
        if name == "within":
            require(set(frame["contrast"]) == {"PSC_treatment", "nonPSC_treatment"}, "Within contrasts mismatch")
            groups = [part for _, part in frame.groupby("contrast", sort=True)]
        else:
            groups = [frame]
        for part in groups:
            require(not part[ID].duplicated().any(), f"Duplicate IDs in {name}")
            require(part[ANNOTATION].reset_index(drop=True).equals(p[ANNOTATION].reset_index(drop=True)),
                    f"Annotation, order or universe mismatch: {name}")
    for name, qcols in (("primary", ["q_bh", "q_by"]), ("within", ["q_bh", "q_by"]),
                        ("welch", ["q_bh_sensitivity", "q_by_sensitivity"])):
        frame = tables[name]
        valid = frame.loc[truth(frame["eligible"])]
        for column in qcols:
            q = number(valid[column])
            require(bool(q.notna().all() and q.between(0, 1).all()), f"Invalid q: {name}/{column}")
        absent = number(valid["pvalue"]).isna()
        for column in qcols:
            require(bool((number(valid.loc[absent, column]) == 1).all()), "Missing P is not correction-only q1")


def contrast_counts(frame: pd.DataFrame, q_bh: str = "q_bh", q_by: str = "q_by",
                    effect: str = "log2FoldChange") -> dict:
    use = frame.loc[truth(frame["eligible"])]
    result = {"family_size": len(use), "finite_pvalues": int(np.isfinite(number(use["pvalue"])).sum())}
    for label, col in (("BH", q_bh), ("BY", q_by)):
        selected = use.loc[number(use[col]) <= 0.05]
        effect_values = number(selected[effect])
        result[label] = {"count": len(selected), "positive": int((effect_values > 0).sum()),
                         "negative": int((effect_values < 0).sum())}
    result["minimum_BH"] = float(number(use[q_bh]).min())
    return result


def summarize(tables: dict[str, pd.DataFrame], qc: dict, verification: dict) -> tuple[dict, dict[str, pd.DataFrame]]:
    validate_tables(tables, qc)
    p = tables["primary"]
    eligible = truth(p["eligible"])
    candidates = p.loc[eligible & (number(p["q_bh"]) <= 0.05)].copy()
    candidates = candidates.assign(_q=number(candidates["q_bh"])).sort_values(["_q", ID]).drop(columns="_q").set_index(ID, drop=False)
    additions = {
        "welch": ["mean_change_PSC", "mean_change_nonPSC", "effect_log2_cpm", "se", "df", "pvalue",
                  "q_bh_sensitivity", "q_by_sensitivity", "ci95_lower", "ci95_upper"],
        "leaveout": ["full_effect_log2_cpm", "minimum_leaveout_effect", "maximum_leaveout_effect",
                     "n_estimable_omissions", "n_same_nonzero_sign", "n_opposite_sign", "n_zero_effect"],
        "allocations": ["n_allocations", "n_as_or_more_extreme", "descriptive_tail_fraction"],
    }
    for name, columns in additions.items():
        other = tables[name].set_index(ID)[columns].add_prefix(name + "_")
        candidates = candidates.join(other, validate="one_to_one")
    for label in ("PSC", "nonPSC"):
        other = tables["within"].loc[tables["within"]["contrast"] == label + "_treatment"].set_index(ID)
        columns = ["log2FoldChange", "lfcSE", "pvalue", "q_bh", "q_by", "ci95_lower_log2fc", "ci95_upper_log2fc"]
        candidates = candidates.join(other[columns].add_prefix("within_" + label + "_"), validate="one_to_one")
    candidates = candidates.reset_index(drop=True)
    flagged_mask = eligible & (p["genewise_converged"].eq("False") | p["MAP_converged"].eq("False") | p["LFC_converged"].eq("False"))
    masked_frames = []
    for name in ("primary", "within"):
        frame = tables[name].copy()
        if name == "primary":
            frame["contrast"] = "interaction"
        masked_frames.append(frame.loc[truth(frame["eligible"]) & frame["count_supported_contrast"].eq("False")])
    mask_columns = ANNOTATION + ["contrast", "status", "zero_required_count_arms", "log2FoldChange", "pvalue", "q_bh", "q_by",
        "pydeseq2_unmasked_log2FoldChange", "pydeseq2_unmasked_lfcSE", "pydeseq2_unmasked_stat",
        "pydeseq2_unmasked_pvalue", "pydeseq2_unmasked_finite_family_padj"]
    masked = pd.concat([frame[mask_columns] for frame in masked_frames], ignore_index=True)
    panel = verification["conditional_fixed_dispersion_NB"]["results"]
    def panel_max(key: str):
        values = [float(row[key]) for row in panel if key in row]
        return max(values) if values else None
    summary = {
        "schema_version": 1,
        "status": "aggregate_reporting_only_no_new_models",
        "plan_sha256": qc["plan_sha256"],
        "feature_coverage": {"released": len(p), "eligible": int(eligible.sum()),
                             "excluded_low_count": int((~eligible).sum()),
                             "all_zero_released": int(truth(p["all_zero_released"]).sum()),
                             "normalization_basis": qc["ratio_basis_positive_in_all_16_libraries"]},
        "units": {key: qc[key] for key in ("donors", "libraries", "released_cells")},
        "interaction": contrast_counts(p),
        "within_group": {label: contrast_counts(tables["within"].loc[tables["within"]["contrast"] == label])
                         for label in ("PSC_treatment", "nonPSC_treatment")},
        "welch_sensitivity": contrast_counts(tables["welch"], "q_bh_sensitivity", "q_by_sensitivity", "effect_log2_cpm"),
        "diagnostics": {
            "genewise_nonconverged": int((eligible & p["genewise_converged"].eq("False")).sum()),
            "MAP_nonconverged": int((eligible & p["MAP_converged"].eq("False")).sum()),
            "LFC_nonconverged": int((eligible & p["LFC_converged"].eq("False")).sum()),
            "any_native_nonconvergence": int(flagged_mask.sum()),
            "eligible_with_cook_above_diagnostic_cutoff": int((eligible & (number(p["n_libraries_cooks_above_diagnostic_cutoff"]) > 0)).sum()),
            "eligible_with_nonfinite_cook": int((eligible & (number(p["n_nonfinite_cooks"]) > 0)).sum()),
            "automatic_cook_exclusions": int((eligible & p["cooks_automatic_exclusion"].eq("True")).sum()),
            "eligible_with_zero_total_donor_pair": int((eligible & ((number(p["n_zero_total_donor_pairs_PSC"]) > 0) | (number(p["n_zero_total_donor_pairs_nonPSC"]) > 0))).sum()),
            "boundary_masks_by_contrast": {key: int(value) for key, value in Counter(masked["contrast"]).items()},
            "masked_raw_package_q_le_0_05": int((number(masked["pydeseq2_unmasked_finite_family_padj"]) <= 0.05).sum()),
            "fit_warning_count": qc["fit_warning_count"],
            "dispersion_fit_type": qc["negative_binomial"]["dispersion_fit_type"],
        },
        "BH_candidate_sensitivities": {
            "same_NB_and_Welch_nonzero_direction": int(((np.sign(number(candidates["log2FoldChange"])) == np.sign(number(candidates["welch_effect_log2_cpm"]))) & (number(candidates["welch_effect_log2_cpm"]) != 0)).sum()),
            "all_eight_leaveout_signs_match_full_paired_change": int((number(candidates["leaveout_n_same_nonzero_sign"]) == 8).sum()),
            "Welch_BH": int((number(candidates["welch_q_bh_sensitivity"]) <= 0.05).sum()),
            "any_native_nonconvergence": int((candidates["genewise_converged"].eq("False") | candidates["MAP_converged"].eq("False") | candidates["LFC_converged"].eq("False")).sum()),
            "any_Cook_above_diagnostic_cutoff": int((number(candidates["n_libraries_cooks_above_diagnostic_cutoff"]) > 0).sum()),
            "allocation_fraction_min": float(number(candidates["allocations_descriptive_tail_fraction"]).min()) if len(candidates) else None,
            "allocation_fraction_max": float(number(candidates["allocations_descriptive_tail_fraction"]).max()) if len(candidates) else None,
        },
        "verification": {
            "status": verification["status"], "checks": len(verification["checks"]),
            "conditional_panel_selected": len(panel),
            "conditional_panel_status": dict(Counter(row["status"] for row in panel)),
            "conditional_optimizer_success_false": sum(row.get("optimizer_success") is False for row in panel),
            "maximum_contrast_difference_in_native_SE": panel_max("maximum_contrast_difference_in_native_SE"),
            "maximum_twice_objective_improvement": panel_max("twice_objective_improvement"),
            "maximum_newton_decrement": panel_max("newton_decrement"),
            "full_dispersion_model_refitted": verification["full_dispersion_model_refitted"],
        },
        "interpretation_limits": [
            "Four PSC and four procedure-control donors; not cell-level replication or healthy-volunteer controls.",
            "NB interaction is a log2 ratio of treatment ratios; positive does not alone mean induction in PSC.",
            "Welch compares equal-donor log2(CPM+1) changes; it is a different estimand, not independent replication.",
            "Leaveout sign stability is descriptive and is not an NB refit, confidence interval or new discovery rule.",
            "Allocation fractions are descriptive; disease labels were not randomized.",
            "BH/BY depend on valid constituent model P values; small-n calibration and clinical response confounding remain.",
            "Native optimizer flags are retained, not used for a new post-effect family or exclusion.",
            "Conditional NB coefficient agreement does not refit dispersions or establish biological replication.",
            "No inference about PSC causation, disease specificity, treatment benefit or equivalence.",
        ],
    }
    return summary, {"bh-candidates-with-diagnostics.csv.gz": candidates,
                     "masked-contrasts.csv.gz": masked,
                     "native-flagged-features.csv.gz": p.loc[flagged_mask].copy()}


def read_verified_inputs(directory: Path, verification_path: Path) -> tuple[dict, dict, dict]:
    manifest_path = directory / "organoid-output-manifest.json"
    manifest = json.loads(manifest_path.read_text())
    verification = json.loads(verification_path.read_text())
    require(verification["status"] == "independent_numerical_verification_passed", "Verification not passed")
    require(digest(manifest_path) == verification["public_manifest_sha256"], "Verified manifest changed")
    require(set(manifest["outputs"]) == set(TABLES.values()) | {"organoid-qc.json"}, "Unexpected public inventory")
    for filename, pin in manifest["outputs"].items():
        path = directory / filename
        require(path.stat().st_size == pin["bytes"] and digest(path) == pin["sha256"], f"Changed input: {filename}")
    qc = json.loads((directory / "organoid-qc.json").read_text())
    require(qc["effects_computed"] is True and manifest["source_values_mutated"] is False, "Wrong execution state")
    require(qc["plan_sha256"] == manifest["plan_sha256"] == verification["plan_sha256"], "Plan mismatch")
    tables = {name: pd.read_csv(directory / filename, dtype=str, keep_default_na=False) for name, filename in TABLES.items()}
    return tables, qc, verification


def destinations(output: Path, figure_prefix: Path, inputs: Path) -> tuple[Path, list[Path]]:
    figures = [figure_prefix.with_suffix(suffix) for suffix in (".png", ".svg")]
    for path in [output] + figures:
        require(not path.exists() and not path.is_symlink(), f"Destination already exists: {path}")
        resolved = path.resolve()
        require(not any(resolved.is_relative_to(p.resolve()) for p in (inputs, ROOT / "data/raw", ROOT / "work", ROOT / "scripts", ROOT / "tests", ROOT / "config")), "Protected output destination")
    require(not output.resolve().is_relative_to(figure_prefix.resolve()), "Overlapping destinations")
    return output.resolve(), figures


def csv_gz(frame: pd.DataFrame, path: Path) -> None:
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            compressed.write(frame.to_csv(index=False, lineterminator="\n", na_rep="", float_format="%.17g").encode())


def make_figure(tables: dict, candidates: pd.DataFrame, paths: list[Path], plan_hash: str) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams.update({"svg.hashsalt": plan_hash, "font.family": "DejaVu Sans", "font.size": 9})
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), layout="constrained")
    p = tables["primary"]
    p = p.loc[truth(p["eligible"]) & np.isfinite(number(p["pvalue"]))]
    classes = [(number(p["q_bh"]) > .05, "Not BH-discovered", "#bbbbbb", 7),
               ((number(p["q_bh"]) <= .05) & (number(p["q_by"]) > .05), "BH only", "#c77c14", 24),
               (number(p["q_by"]) <= .05, "BH and BY", "#087f8c", 28)]
    ax = axes[0, 0]
    for mask, label, color, size in classes:
        part = p.loc[mask]
        ax.scatter(number(part["log2FoldChange"]), -np.log10(np.maximum(number(part["pvalue"]), np.finfo(float).tiny)), s=size, c=color, alpha=.65, linewidths=0, label=f"{label} (n={len(part):,})")
    ax.axvline(0, color="black", lw=.6)
    ax.set(xlabel="NB interaction: log2 ratio of treatment ratios (PSC / nonPSC)", ylabel="−log10 raw model P", title=f"A. Full supported interaction family (n={len(p):,})")
    ax.legend(frameon=False, fontsize=8)
    labels = candidates["source_gene_label"].tolist()
    positions = np.arange(len(candidates))
    colors = np.where(number(candidates["q_by"]) <= .05, "#087f8c", "#c77c14")
    for ax in (axes[0, 1], axes[1, 0], axes[1, 1]):
        ax.set_yticks(positions, labels)
        ax.invert_yaxis()
        ax.axvline(0, color="black", lw=.6)
        ax.grid(axis="x", alpha=.2)
    ax = axes[0, 1]
    for i, row in candidates.iterrows():
        ax.plot([float(row["ci95_lower_log2fc"]), float(row["ci95_upper_log2fc"])], [i, i], color=colors[i], lw=1.4)
        ax.scatter(float(row["log2FoldChange"]), i, color=colors[i], s=26)
    ax.set(title=f"B. {len(candidates)} primary BH candidates: NB estimates and 95% Wald CIs", xlabel="log2 ratio of treatment ratios")
    ax = axes[1, 0]
    for i, row in candidates.iterrows():
        ax.plot([float(row["welch_ci95_lower"]), float(row["welch_ci95_upper"])], [i, i], color="#555555", lw=1.4)
        ax.scatter(float(row["welch_effect_log2_cpm"]), i, color="#555555", s=26)
    welch_bh = contrast_counts(tables["welch"], "q_bh_sensitivity", "q_by_sensitivity", "effect_log2_cpm")["BH"]["count"]
    ax.set(title=f"C. Paired-change sensitivity: {welch_bh} BH discoveries", xlabel="PSC − nonPSC mean donor change in log2(CPM+1); 95% t CIs")
    ax = axes[1, 1]
    for label, shift, color, marker in (("nonPSC", -.13, "#56568c", "o"), ("PSC", .13, "#b45a3a", "^")):
        for i, row in candidates.iterrows():
            prefix = "within_" + label + "_"
            ax.plot([float(row[prefix + "ci95_lower_log2fc"]), float(row[prefix + "ci95_upper_log2fc"])], [i + shift, i + shift], color=color, lw=1)
        ax.scatter(number(candidates["within_" + label + "_log2FoldChange"]), positions + shift, color=color, marker=marker, s=26, label=label)
    ax.set(title="D. Within-group model responses (not a separate interaction test)", xlabel="IL17 / CNT log2 fold change; 95% Wald CIs")
    ax.legend(frameon=False)
    fig.suptitle("GSE239283: 4 PSC + 4 procedure-control donors; model-dependent ex-vivo response associations", fontsize=12)
    fig.supxlabel("Panels B–D show the same post-result BH-selected genes. CIs are unadjusted. Leaveout stability is not replication. No causal or clinical inference.", fontsize=9)
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)
        metadata = {"Date": None, "Creator": "PSC organoid aggregate report"} if path.suffix == ".svg" else {"Software": "PSC organoid aggregate report"}
        fig.savefig(path, dpi=180, metadata=metadata)
    plt.close(fig)
    return matplotlib.__version__


def run(inputs: Path, verification_path: Path, output: Path, figure_prefix: Path) -> dict:
    output, figure_paths = destinations(output, figure_prefix, inputs)
    tables, qc, verification = read_verified_inputs(inputs, verification_path)
    summary, exports = summarize(tables, qc, verification)
    output.mkdir(parents=True)
    for filename, frame in exports.items():
        csv_gz(frame, output / filename)
    mpl_version = make_figure(tables, exports["bh-candidates-with-diagnostics.csv.gz"], figure_paths, qc["plan_sha256"])
    summary["reporting_software"] = {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "matplotlib": mpl_version}
    with (output / "summary.json").open("x") as stream:
        json.dump(summary, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    manifest = {"report_script_sha256": digest(Path(__file__)), "plan_sha256": qc["plan_sha256"],
                "source_public_manifest_sha256": digest(inputs / "organoid-output-manifest.json"),
                "independent_verification_sha256": digest(verification_path),
                "outputs": {p.name: {"sha256": digest(p), "bytes": p.stat().st_size} for p in list(output.iterdir()) + figure_paths},
                "patient_level_inputs_read": False, "new_tests_or_models": False}
    with (output / "report-output-manifest.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--verification-json", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figure-prefix", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.input_dir, args.verification_json, args.output_dir, args.figure_prefix)
    print(json.dumps({key: result[key] for key in ("status", "interaction", "welch_sensitivity")}, sort_keys=True))


if __name__ == "__main__":
    main()
