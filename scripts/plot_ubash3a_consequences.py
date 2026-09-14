"""Plot the verified UBASH3A sequence consequences and conditional NMD features."""

import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/derived/ubash3a-consequences"


def main():
    audit = json.loads((DATA / "audit.json").read_text())
    for filename in ("consequences.csv", "protein-features.csv"):
        assert hashlib.sha256((DATA / filename).read_bytes()).hexdigest() == audit["files_sha256"][filename]
    rows = list(csv.DictReader((DATA / "consequences.csv").open()))
    features = list(csv.DictReader((DATA / "protein-features.csv").open()))
    canonical = [r for r in rows if r["is_canonical"] == "1"]
    reference = next(r for r in canonical if r["scenario"] == "normal_A")
    altered = next(r for r in canonical if r["scenario"] == "plus29_A")
    length, prefix, product = (int(reference["protein_length_aa"]), int(altered["unchanged_prefix_aa"]), int(altered["protein_length_aa"]))
    palette = {"UBA": "#167D7F", "SH3": "#356AAC", "Phosphatase-like": "#7960A8"}
    novel, foreground, gray = "#C65C28", "#172D40", "#D9E1E7"
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "text.color": foreground, "axes.labelcolor": foreground, "xtick.color": foreground, "ytick.color": foreground, "svg.hashsalt": "ubash3a-consequences-v1"})
    fig = plt.figure(figsize=(12.8, 9.1), facecolor="white")
    top = fig.add_axes([0.24, 0.55, 0.68, 0.29])
    bottom = fig.add_axes([0.24, 0.17, 0.68, 0.24])
    fig.text(0.06, 0.951, "UBASH3A: what follows if 29 RNA bases are retained?", fontsize=20, weight="bold")
    fig.text(0.06, 0.91, "Sequence reconstruction of five transcript models  •  GRCh38 / Ensembl 116  •  14 September 2026", fontsize=10.5)
    fig.text(0.06, 0.866, "A", fontsize=15, weight="bold")
    fig.text(0.09, 0.866, "Canonical protein sequence: 661 → 536 amino acids", fontsize=14, weight="bold")
    for y, stop in [(2, length), (1, product), (0, product)]:
        top.add_patch(Rectangle((0, y - .17), stop, .34, facecolor=gray, edgecolor="none"))
        if y != 2:
            top.add_patch(Rectangle((prefix, y - .17), stop - prefix, .34, facecolor=novel, edgecolor="none", hatch="///"))
        for feature in [f for f in features if f["scenario_id"] == reference["scenario_id"]]:
            lo, hi = int(feature["uniprot_start_1based"]) - 1, int(feature["uniprot_end_1based"])
            hi = min(hi, prefix) if y != 2 else hi
            if hi > lo:
                top.add_patch(Rectangle((lo, y - .17), hi - lo, .34, facecolor=palette[feature["feature_name"]], edgecolor="none"))
            if y == 2:
                top.text((lo + hi) / 2, 2.34, f"{feature['feature_name']}\n{lo + 1}–{hi}", ha="center", fontsize=10, color=palette[feature["feature_name"]])
        top.plot([stop, stop], [y - .25, y + .25], color=foreground, lw=1.4)
        top.text(stop + 9, y, f"{stop} aa", va="center", fontsize=10)
    top.text(165, .51, f"First {prefix} amino acids unchanged", fontsize=10, ha="center")
    top.text(prefix + (product - prefix) / 2, -.48, f"{product - prefix}-aa new tail", ha="center", fontsize=10, color=novel)
    top.annotate("Tail begins M with A, L with C", xy=(prefix + 15, 1), xytext=(565, .64), fontsize=9, ha="left", arrowprops={"arrowstyle": "-", "color": foreground, "lw": .7})
    top.set_yticks([2, 1, 0], ["Normal A or C", "+29 bases, allele A", "+29 bases, allele C"])
    top.set_xlim(0, 730)
    top.set_ylim(-.77, 2.92)
    top.set_xticks([1, 100, 200, 300, 400, 465, 536, 661])
    top.set_xlabel("Amino-acid position in the constructed product", labelpad=8)
    top.tick_params(axis="y", length=0, pad=12)
    for side in ("top", "right", "left"):
        top.spines[side].set_visible(False)
    top.spines["bottom"].set_color("#BECBD5")
    fig.text(0.06, .459, "B", fontsize=15, weight="bold")
    fig.text(.09, .459, "RNA-decay context differs across transcript models", fontsize=14, weight="bold")
    selected = [next(r for r in rows if r["transcript_id"] == tid and r["scenario"] == "plus29_A") for tid in ("ENST00000319294", "ENST00001131803", "ENST00000291535", "ENST00000398367", "ENST00000635325")]
    bottom.axvspan(0, 55, color="#FCF1E7", zorder=0)
    bottom.axvline(55, color=novel, linestyle="--", lw=1.2)
    bottom.text(61, -.71, "Pre-set flag: >55 nt", color=novel, fontsize=10)
    for i, row in enumerate(selected):
        distance = int(row["stop_end_to_final_junction_nt"])
        color = novel if distance <= 55 else "#356AAC"
        bottom.hlines(i, 0, distance, color=color, lw=3, alpha=.7)
        bottom.scatter([distance], [i], color=color, s=45, zorder=3)
        suffix = "  ·  boundary case" if distance <= 55 else "  ·  flag met"
        bottom.text(distance + 8, i, f"{distance} nt{suffix}", va="center", fontsize=10, color=color)
    bottom.set_yticks(range(5), [r["transcript_id"] + "." + r["transcript_version"] + (" *" if r["biotype"] == "nonsense_mediated_decay" else "") for r in selected])
    bottom.set_xlim(0, 610)
    bottom.set_ylim(4.6, -.9)
    bottom.set_xticks([0, 50, 100, 200, 300, 400, 500])
    bottom.set_xlabel("Distance from stop codon’s final base to the final exon junction (nt)", labelpad=8)
    bottom.tick_params(axis="y", length=0, pad=12)
    bottom.xaxis.grid(True, color="#EAF0F4", lw=.7)
    bottom.set_axisbelow(True)
    for side in ("top", "right", "left"):
        bottom.spines[side].set_visible(False)
    bottom.spines["bottom"].set_color("#BECBD5")
    fig.text(.06, .075, "A: Colors show UniProt reference regions and their unchanged sequence; the hatched tail uses a new reading frame.", fontsize=9)
    fig.text(.06, .055, "B: A and C have the same stop positions. * Already annotated for nonsense-mediated decay before +29 retention.", fontsize=9)
    fig.text(.06, .028, "Conditional sequence products only. Neither protein accumulation nor RNA degradation has been measured here.", fontsize=10, weight="bold")
    destination = ROOT / "reports/figures"
    destination.mkdir(exist_ok=True)
    for suffix in ("png", "svg"):
        output = destination / ("ubash3a-transcript-consequences." + suffix)
        fig.savefig(output, dpi=190, metadata={"Date": None} if suffix == "svg" else {"Software": "Matplotlib; UBASH3A consequence audit"})
        if suffix == "svg":
            output.write_text("\n".join(line.rstrip() for line in output.read_text().splitlines()) + "\n")
    plt.close(fig)
    print("Saved reports/figures/ubash3a-transcript-consequences.{png,svg}")


if __name__ == "__main__":
    main()
