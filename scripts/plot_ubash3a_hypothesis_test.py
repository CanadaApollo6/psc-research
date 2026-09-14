"""Render the RNA-structure comparison and deposited CD4 support values."""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports/figures"


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "svg.hashsalt": "ubash3a-hypothesis-test-v1"})
    fig = plt.figure(figsize=(14, 9), facecolor="#fbfcfe")
    grid = fig.add_gridspec(2, 1, height_ratios=[1.35, 1], left=.19, right=.96, top=.82, bottom=.16, hspace=.55)
    ax = fig.add_subplot(grid[0])
    ax.set_xlim(-.3, 10.1)
    ax.set_ylim(-.5, 2.7)
    ax.axis("off")
    navy, gold, teal, gray = "#305b89", "#e7ac3e", "#087f8c", "#9aa7b4"
    def exon(x, y, width=.55, color=navy):
        ax.add_patch(Rectangle((x-width/2, y-.12), width, .24, facecolor=color, edgecolor="none", zorder=3))
    def line(a, b, y, color=gray):
        ax.plot([a, b], [y, y], color=color, lw=1.6, zorder=1)
    positions = [.4, 2, 4.1, 5.9, 7.2, 9]
    for y, label in [(2, "Canonical-based +29\nconditional reference"), (1, "Boundary +29\nconditional reference")]:
        ax.text(-.7, y, label, ha="right", va="center", fontsize=11, color="#182b3d")
        line(positions[0], positions[-1], y)
        xs = positions if y == 2 else [positions[0], positions[1], positions[2], positions[-1]]
        for x in xs:
            exon(x, y)
        ax.add_patch(Rectangle((positions[0]+.10, y-.12), .175, .24, facecolor=gold, edgecolor="none", zorder=4))
        ax.plot(4.05, y, marker="X", color="#a33a43", markersize=9, zorder=5)
    ax.text(.4, 2.42, "+29 splice choice", ha="center", fontsize=10)
    ax.text(2, 2.42, "Shared exon", ha="center", fontsize=10)
    ax.text(4.1, 2.42, "Predicted stop", ha="center", fontsize=10)
    ax.text(6.55, 2.42, "Extra downstream exons", ha="center", fontsize=10)
    ax.text(9, 2.42, "Final exon", ha="center", fontsize=10)
    ax.text(6.7, .72, "Skip to final exon: 49 nt from stop end to last junction", ha="center", fontsize=9, color=navy)
    ax.text(-.7, 0, "TRAILS +29\ncatalog candidate", ha="right", va="center", fontsize=11, color="#182b3d")
    line(.4, 4.85, 0, teal)
    exon(.4, 0, color=teal)
    exon(2, 0, color=teal)
    ax.add_patch(Rectangle((.5, -.12), .175, .24, facecolor=gold, edgecolor="none", zorder=4))
    exon(4.85, 0, width=2.45, color=teal)
    ax.text(4.85, -.33, "Different terminal exon: chr21:42437855–42440816", ha="center", fontsize=9, color=teal)
    ax.text(7.6, 0, "Ends before the reference stop region", ha="left", va="center", fontsize=10, color=teal)
    ax.text(0, 1.1, "A  •  RNA endings differ", transform=ax.transAxes, fontsize=13, fontweight="bold", ha="left")
    ax.text(0, 1.005, "Schematic: upstream exons omitted; exon and intron lengths are not to scale.", transform=ax.transAxes, fontsize=9, color="#627487")
    counts = [r for r in csv.DictReader((ROOT / "data/derived/ubash3a-hypothesis-test/trails-locus-source-counts.csv").open()) if r["event"] == "plus29" and r["primary_CD4_context"] == "True"]
    bx = fig.add_subplot(grid[1])
    labels = [r["cell_subset"] for r in counts]
    values = [float(r["source_count_value"]) for r in counts]
    bx.bar(range(len(values)), values, color=teal, width=.62)
    bx.set_xticks(range(len(labels)), labels, rotation=28, ha="right", fontsize=9)
    bx.set_ylim(0, max(values)*1.25)
    bx.set_ylabel("Source count units*", fontsize=10)
    bx.spines[["top", "right"]].set_visible(False)
    bx.spines[["left", "bottom"]].set_color("#bbc7d1")
    bx.grid(axis="y", color="#e4eaf0", zorder=0)
    bx.set_axisbelow(True)
    for i, value in enumerate(values):
        bx.text(i, value+1.2, str(int(value)), ha="center", fontsize=9)
    bx.text(0, 1.2, "B  •  Catalog candidate has assigned support in 10 of 11 CD4 subsets", transform=bx.transAxes, fontsize=12, fontweight="bold")
    bx.text(0, 1.08, "All subsets are from one donor. This is not an NMD-perturbation experiment.", transform=bx.transAxes, fontsize=10, color="#627487")
    fig.text(.05, .95, "The +29 catalog candidate uses a different ending", fontsize=22, weight="bold", color="#182b3d")
    fig.text(.05, .908, "No inspected annotation contains the predicted boundary chain, with or without +29.", fontsize=12, color="#4c6275")
    fig.text(.05, .073, "*Original TRAILS matrix values. The precise export/normalization step is not fully documented; these are not unique molecules or donors.", fontsize=9, color="#536779")
    fig.text(.05, .044, "Direct UPF1 test: eight libraries, one unspliced alignment overlapping UBASH3A, no target-junction reads. Decay contrasts are not estimable.", fontsize=9, color="#536779")
    OUT.mkdir(parents=True, exist_ok=True)
    png = OUT / "ubash3a-hypothesis-test.png"
    svg = OUT / "ubash3a-hypothesis-test.svg"
    fig.savefig(png, dpi=150, facecolor=fig.get_facecolor(), metadata={"Software": "PSC public-data research"})
    fig.savefig(svg, facecolor=fig.get_facecolor(), metadata={"Date": None, "Creator": "PSC public-data research"})
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
    plt.close(fig)
    print(str(png))


if __name__ == "__main__":
    main()
