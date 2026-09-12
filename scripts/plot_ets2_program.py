"""Plot recorded aggregate ETS2 program contrasts; never read donor records."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "work/matplotlib"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

PROGRAMS = (
    ("ets2_g1_dn", "A  ETS2-dependent program"),
    ("ets2_g1_dn_without_comparators", "B  After removing overlap"),
    ("inflammation", "C  Inflammatory response"),
)
CONTEXTS = tuple((dataset, population) for dataset in ("sc", "sn")
                 for population in ("ActMac", "LAM-like", "Monocyte"))


def main() -> None:
    source = ROOT / "data/derived/ets2-program-contrasts.csv"
    rows = pd.read_csv(source)
    required = {"dataset_key", "left_population", "right_population", "program_id",
                "threshold", "status", "paired_donors", "median_delta",
                "minimum_delta", "maximum_delta"}
    if not required.issubset(rows.columns):
        raise ValueError("Recorded contrast table lacks required columns")
    selected = rows[rows.program_id.isin([p for p, _ in PROGRAMS])]
    keys = ["dataset_key", "left_population", "program_id", "threshold"]
    if selected.duplicated(keys).any() or len(selected) != 36:
        raise ValueError("Figure requires exactly one assay per declared contrast and 36 rows")
    if not selected.right_population.eq("Kupffer").all():
        raise ValueError("Unexpected contrast baseline")
    finite = selected[selected.status.eq("evaluable_descriptive")]
    limit = float(np.nanmax(np.abs(finite[["minimum_delta", "maximum_delta"]].to_numpy())))
    limit = max(0.1, limit) * 1.15
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.spines.left": False, "axes.titleweight": "bold"})
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 7.1), sharex=True, sharey=True)
    colors = {20: "#176C82", 50: "#C47A25"}
    for ax, (program, title) in zip(axes, PROGRAMS):
        ax.set_title(title, loc="left", fontsize=11, pad=15)
        ax.axvline(0, color="#687681", lw=1, zorder=0)
        ax.axhspan(2.5, 5.5, color="#F0F3F5", zorder=-2)
        for y, (dataset, population) in enumerate(CONTEXTS):
            for threshold, offset in ((20, -0.12), (50, 0.12)):
                hit = selected[(selected.program_id == program)
                               & (selected.dataset_key == dataset)
                               & (selected.left_population == population)
                               & (selected.threshold == threshold)]
                if len(hit) != 1:
                    raise ValueError("Missing planned figure row")
                row = hit.iloc[0]
                n = int(row.paired_donors)
                if row.status == "evaluable_descriptive":
                    lo, mid, hi = float(row.minimum_delta), float(row.median_delta), float(row.maximum_delta)
                    if not np.isfinite([lo, mid, hi]).all() or not lo <= mid <= hi:
                        raise ValueError("Invalid recorded range")
                    ax.hlines(y + offset, lo, hi, color=colors[threshold], lw=1.5)
                    ax.scatter([mid], [y + offset], s=30, color=colors[threshold],
                               facecolor=colors[threshold] if threshold == 20 else "white", zorder=3)
                else:
                    ax.text(0.04, y + offset, "unavailable", transform=ax.get_yaxis_transform(),
                            ha="left", va="center", fontsize=8, color="#697782")
                ax.text(1.02, y + offset, f"n={n}", transform=ax.get_yaxis_transform(),
                        ha="left", va="center", fontsize=8, color=colors[threshold])
        ax.set_xlim(-limit, limit)
        ax.set_ylim(5.6, -0.6)
        ax.grid(axis="x", alpha=0.15)
        ax.tick_params(axis="y", length=0)
        ax.set_xlabel("Within-donor score difference", labelpad=10)
    axes[0].set_yticks(range(len(CONTEXTS)), [
        "Cells · activated macrophages", "Cells · LAM-like macrophages", "Cells · monocytes",
        "Nuclei · activated macrophages", "Nuclei · LAM-like macrophages", "Nuclei · monocytes",
    ])
    fig.suptitle("PSC liver: paired program contrasts against Kupffer cells", x=0.03,
                 y=0.97, ha="left", fontsize=16, weight="bold")
    fig.text(0.03, 0.915, "Points: donor-pair median  •  Lines: full donor range  •  n: eligible donor pairs",
             fontsize=10, color="#43525E")
    fig.legend(handles=[Line2D([0], [0], color=colors[t], marker="o",
                               markerfacecolor=colors[t] if t == 20 else "white",
                               label=f"At least {t} cells per population") for t in (20, 50)],
               loc="lower left", bbox_to_anchor=(0.027, 0.086), frameon=False, ncol=2)
    fig.text(0.03, 0.055,
             "All contrasts are against Kupffer cells in the same donor and assay. Ranges are not confidence intervals.",
             fontsize=9, color="#43525E")
    fig.text(0.03, 0.025,
             "Scores use mean log2(CPM + 1) relative to expression-matched genes. These are descriptive differences, not fold changes or treatment effects.",
             fontsize=9, color="#43525E")
    fig.subplots_adjust(left=0.225, right=0.945, top=0.83, bottom=0.24, wspace=0.3)
    out = ROOT / "reports/figures"
    out.mkdir(parents=True, exist_ok=True)
    files = [out / "ets2-program-contrasts.png", out / "ets2-program-contrasts.pdf"]
    fig.savefig(files[0], dpi=190, facecolor="white", metadata={"Software": "psc-research"})
    fig.savefig(files[1], facecolor="white", metadata={"CreationDate": None, "ModDate": None,
                                                     "Creator": "psc-research"})
    plt.close(fig)
    hashes = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    audit = {"source_file": str(source.relative_to(ROOT)), "source_sha256": hashes(source),
             "script_sha256": hashes(Path(__file__)), "plotted_rows": len(selected),
             "selection": "Primary, nonoverlap sensitivity and inflammation programs; all 3 contrasts, both modalities and thresholds",
             "outputs_sha256": {str(p.relative_to(ROOT)): hashes(p) for p in files}}
    (ROOT / "reports/ets2-program-figure-audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({"figure_files": list(audit["outputs_sha256"]), "plotted_rows": len(selected)}))


if __name__ == "__main__":
    main()
