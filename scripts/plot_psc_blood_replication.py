#!/usr/bin/env python3
"""Plot aggregate PSC blood replication; the two assay effect axes differ."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results", default="data/derived/psc-blood-replication/cross-cohort-gene-results.csv.gz")
    p.add_argument("--summary", default="data/derived/psc-blood-replication/summary.json")
    p.add_argument("--output-prefix", default="reports/figures/psc-blood-replication")
    args = p.parse_args()
    data = pd.read_csv(args.results)
    summary = json.loads(Path(args.summary).read_text())
    available = np.isfinite(data.rna_effect) & np.isfinite(data.array_control_effect)
    hit = data.replication_qvalue_bh.le(0.05)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
    ax = axes[0]
    ordinary = available & ~hit
    ax.scatter(data.loc[ordinary, "rna_effect"], data.loc[ordinary, "array_control_effect"], s=5, color="#8c959f", alpha=0.25, rasterized=True, label="Other shared genes")
    selected = available & hit
    ax.scatter(data.loc[selected, "rna_effect"], data.loc[selected, "array_control_effect"], s=13, color="#bb4430", alpha=0.75, rasterized=True, label="Conjunction BH q ≤ 0.05")
    ax.axhline(0, color="black", linewidth=0.6)
    ax.axvline(0, color="black", linewidth=0.6)
    ax.set(xlabel="Norway: adjusted RNA log2 fold change", ylabel="Poland: log2 normalized-intensity difference", title="PSC minus healthy controls")
    ax.legend(fontsize=8, frameon=False, loc="best")
    ax.text(0.02, 0.02, "Different assay scales; no 1:1 effect comparison", transform=ax.transAxes, fontsize=8)
    family_keys = ["replication", "uc_support", "all_disease_support"]
    labels = ["Replicated\nPSC–control", "+ PSC–UC\ncontrast support", "+ PSC–PBC/CD\ncontrast support"]
    bh = [summary["families"][key]["q_bh_le_0_05"] for key in family_keys]
    by = [summary["families"][key]["q_by_le_0_05"] for key in family_keys]
    x = np.arange(3)
    axes[1].bar(x - 0.18, bh, width=0.36, color="#bb4430", label="BH q ≤ 0.05")
    axes[1].bar(x + 0.18, by, width=0.36, color="#2c5777", label="BY q ≤ 0.05")
    for xx, yy in zip(x - 0.18, bh):
        axes[1].annotate(str(yy), (xx, yy), xytext=(0, 3), textcoords="offset points", ha="center", fontsize=8)
    for xx, yy in zip(x + 0.18, by):
        axes[1].annotate(str(yy), (xx, yy), xytext=(0, 3), textcoords="offset points", ha="center", fontsize=8)
    axes[1].set_xticks(x, labels)
    axes[1].set_ylabel("Genes passing direction-consistent conjunction")
    axes[1].set_title(f"Complete shared family: {len(data):,} genes")
    axes[1].legend(fontsize=8, frameon=False)
    axes[1].margins(y=0.18)
    fig.suptitle("PSC blood cross-cohort reanalysis — association, not causality", fontsize=13)
    out = Path(args.output_prefix)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out) + ".png", dpi=180)
    fig.savefig(str(out) + ".svg")
    plt.close(fig)


if __name__ == "__main__":
    main()
