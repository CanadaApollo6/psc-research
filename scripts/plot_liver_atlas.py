"""Plot selected PSC cell contexts from the complete five-gene atlas summary.

The displayed cell labels were selected for immune/biliary/stromal relevance;
all author labels, sparse groups and excluded groups remain in the full tables.
"""

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "work/matplotlib-liver-atlas"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

GENES = ["ETS2", "PRKD2", "UBASH3A", "PFKFB3", "BCL2L11"]
CONTEXTS = ["CD4T", "Monocyte", "Kupffer", "LAM-like", "ActMac", "Chol",
            "CholMucus", "Stellate", "aStellate", "Fibroblast"]
LABELS = {"CD4T": "CD4 T cells", "Monocyte": "Monocytes", "Kupffer": "Kupffer cells",
          "LAM-like": "LAM-like macrophages", "ActMac": "Activated macrophages",
          "Chol": "Cholangiocytes", "CholMucus": "Mucus cholangiocytes",
          "Stellate": "Stellate cells", "aStellate": "Activated stellate cells",
          "Fibroblast": "Fibroblasts"}


def main():
    table = pd.read_csv(ROOT / "data/derived/liver-atlas-expression-summary.csv")
    table = table[(table.cohort == "PSC") & (table.cell_threshold == 20)]
    selection = table[table.author_cell_type.isin(CONTEXTS)]
    if selection.duplicated(["dataset_key", "author_cell_type", "gene"]).any():
        raise ValueError("Multiple assay strata require separate plotting; do not pool them")
    maximum = max(1., float(selection.donor_weighted_median_log2_cpm1.max()))
    norm = Normalize(0, np.ceil(maximum))
    cmap = plt.get_cmap("viridis")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 8.4), gridspec_kw={"wspace": .70})
    fig.subplots_adjust(left=.20, right=.96, top=.80, bottom=.24)
    for ax, mode, title in zip(axes, ["sc", "sn"], ["Whole cells · 10x 5′ v2", "Nuclei · 10x 3′ v3"]):
        subset = selection[selection.dataset_key == mode]
        labels = []
        for y, context in enumerate(CONTEXTS):
            rows = subset[subset.author_cell_type == context]
            n = int(rows.eligible_donors.max()) if len(rows) else 0
            labels.append(f"{LABELS[context]}  (n={n})")
            for x, gene in enumerate(GENES):
                row = rows[rows.gene == gene]
                if len(row) != 1 or n < 3:
                    ax.text(x, y, "–", color="#9a9a9a", ha="center", va="center")
                    continue
                r = row.iloc[0]
                detection = float(r.donor_weighted_median_detection_fraction)
                abundance = float(r.donor_weighted_median_log2_cpm1)
                if not np.isfinite(detection) or not np.isfinite(abundance):
                    ax.text(x, y, "–", color="#9a9a9a", ha="center", va="center")
                elif detection == 0:
                    ax.scatter(x, y, s=10, facecolors="none", edgecolors="#7c7c7c")
                else:
                    ax.scatter(x, y, s=12 + 420*detection, c=[abundance], cmap=cmap, norm=norm,
                               edgecolors="#ffffff", linewidths=.5)
        ax.set_title(title, loc="left", fontsize=12, fontweight="bold", pad=13)
        ax.set_yticks(range(len(CONTEXTS)), labels=labels, fontsize=9)
        ax.set_xticks(range(len(GENES)), labels=GENES, rotation=35, ha="right")
        ax.set_xlim(-.6, len(GENES)-.4)
        ax.set_ylim(len(CONTEXTS)-.4, -.6)
        ax.set_axisbelow(True)
        ax.grid(axis="y", color="#ededed", linewidth=.7)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(length=0)
    fig.text(.04, .95, "Where our candidate genes are measured in PSC liver", fontsize=17, weight="bold")
    fig.text(.04, .906, "RNA detection and abundance across donors in selected cell populations", fontsize=11, color="#555555")
    fig.text(.04, .864, "Public atlas reanalysis · late-stage PSC · descriptive cell context", fontsize=10, color="#666666")
    color_ax = fig.add_axes([.24, .137, .24, .018])
    fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=color_ax, orientation="horizontal")
    color_ax.set_xlabel("Median donor log₂(CPM + 1)", fontsize=9)
    legend_ax = fig.add_axes([.59, .10, .36, .085])
    legend_ax.set_xlim(0, 4); legend_ax.set_ylim(0, 1); legend_ax.axis("off")
    for x, fraction in zip([.5, 1.5, 2.5], [.05, .25, .75]):
        legend_ax.scatter(x, .65, s=12+420*fraction, color="#777777")
        legend_ax.text(x, .12, f"{int(fraction*100)}%", ha="center", fontsize=9)
    legend_ax.text(1.5, -.25, "Median donor fraction of cells with detected RNA", ha="center", fontsize=9)
    fig.text(.04, .045, "n = donors with ≥20 cells in that population. Dashes: absent population or fewer than 3 eligible donors. Hollow dots: median detection is zero.\n"
             "Panels share some donors and use different assays; they are not independent replications. Full tables retain all author labels.",
             fontsize=8, color="#555555", linespacing=1.6)
    fig.text(.04, .012, "Source: Andrews et al., J Hepatol (2024); CELLxGENE versions ff713455 / 34702b40. Reanalysis: 2026-09-12.",
             fontsize=8, color="#555555")
    for extension in ("png", "pdf"):
        fig.savefig(ROOT / f"reports/liver-atlas-candidate-contexts.{extension}", dpi=180,
                    metadata={"Creator": "PSC public-data research", "CreationDate": None,
                              "ModDate": None} if extension == "pdf" else None)
    plt.close(fig)
    print("Wrote reports/liver-atlas-candidate-contexts.png and .pdf")


if __name__ == "__main__":
    main()
