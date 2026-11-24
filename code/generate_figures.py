#!/usr/bin/env python3
"""
code/generate_figures.py

Data-driven figure generator for the ReScience C SADDLE replication manuscript.
Every plotted value is loaded directly from canonical_results.json and
saddle_384_dimer_pairs.tsv.
"""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def setup_style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 8.5,
        "axes.titlesize": 9.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 8.5,
        "axes.edgecolor": "#222222",
        "axes.linewidth": 0.8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7.5,
        "legend.frameon": True,
        "legend.edgecolor": "#CCCCCC",
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def generate_fig1_multipanel(canon: dict, out_dir: Path) -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.4), gridspec_kw={"width_ratios": [1.0, 1.35]})

    ps_ids = ["PS1", "PS2", "PS3"]
    pub_dimers = [canon["panels"][p]["published_reference"]["dimer_pct"] for p in ps_ids]
    bt2_dimers = [canon["panels"][p]["pipeline_a_bowtie2"]["pct_of_classified_ipynb"]["dimer_pct"] for p in ps_ids]
    py_dimers = [canon["panels"][p]["pipeline_b_edlib"]["method2_edlib_75"]["dimer_pct_classified"] for p in ps_ids]

    x = np.arange(len(ps_ids))
    w = 0.25
    b1 = ax1.bar(x - w, pub_dimers, width=w, label="Published (Xie et al., Fig. 3a)", color="#222222", edgecolor="black", linewidth=0.5)
    b2 = ax1.bar(x, bt2_dimers, width=w, label="Pipeline A: Python + Bowtie2", color="#1F77B4", edgecolor="black", linewidth=0.5)
    b3 = ax1.bar(x + w, py_dimers, width=w, label="Pipeline B: Python/edlib", color="#2CA02C", edgecolor="black", linewidth=0.5)

    for bars in (b1, b2, b3):
        for bar in bars:
            h = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width() / 2.0, h + 1.4, f"{h:.1f}%", ha="center", va="bottom", fontsize=7.0)

    ax1.set_xticks(x)
    ax1.set_xticklabels(["PS1\n(Unoptimized)", "PS2\n(Intermediate)", "PS3\n(Optimized)"])
    ax1.set_ylabel("Primer-Dimer Fraction (% of Classified Reads)")
    ax1.set_ylim(0, 108)
    ax1.set_title("A. 96-Plex Dimer Suppression (Fig. 3a Replication)", loc="left", pad=8)
    ax1.grid(axis="y", linestyle="--", alpha=0.35)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.legend(loc="upper right")

    e19 = canon["panels"]["E19"]
    pub = e19["published_reference"]
    pa = e19["pipeline_a_bowtie2"]["pct_of_classified_ipynb"]
    pb = e19["pipeline_b_edlib"]

    methods = [
        ("Published (Fig. 3f)", [pub["ontarget_pct"], pub["nonspecific_pct"], pub["dimer_pct"]], "#222222"),
        ("Pipeline A (Bowtie2)", [pa["ontarget_pct"], pa["nonspecific_pct"], pa["dimer_pct"]], "#1F77B4"),
        ("Pipeline B M1 (Primer Pair)", [pb["method1_primer_pairing"]["ontarget_pct_classified"], pb["method1_primer_pairing"]["nonspecific_pct_classified"], pb["method1_primer_pairing"]["dimer_pct_classified"]], "#6BAED6"),
        ("Pipeline B M2 (edlib ≥75%)", [pb["method2_edlib_75"]["ontarget_pct_classified"], pb["method2_edlib_75"]["nonspecific_pct_classified"], pb["method2_edlib_75"]["dimer_pct_classified"]], "#2CA02C"),
        ("Pipeline B M2 (edlib ≥80%)", [pb["method2_edlib_80"]["ontarget_pct_classified"], pb["method2_edlib_80"]["nonspecific_pct_classified"], pb["method2_edlib_80"]["dimer_pct_classified"]], "#98DF8A"),
        ("Pipeline B M2 (edlib ≥85%)", [pb["method2_edlib_85"]["ontarget_pct_classified"], pb["method2_edlib_85"]["nonspecific_pct_classified"], pb["method2_edlib_85"]["dimer_pct_classified"]], "#FF7F0E"),
        ("Ungapped (≥75% Ref)", [pb["ungapped_positional_75"]["ontarget_pct_classified"], pb["ungapped_positional_75"]["nonspecific_pct_classified"], pb["ungapped_positional_75"]["dimer_pct_classified"]], "#D62728"),
    ]

    cats = ["On-Target\nAmplicons", "Nonspecific\nProducts", "Primer-\nDimers"]
    xc = np.arange(len(cats))
    n_m = len(methods)
    bw = 0.80 / n_m

    for i, (label, vals, color) in enumerate(methods):
        offset = (i - (n_m - 1) / 2.0) * bw
        bars = ax2.bar(xc + offset, vals, width=bw, label=label, color=color, edgecolor="black", linewidth=0.35)
        for j, bar in enumerate(bars):
            h = bar.get_height()
            if j < 2:
                ax2.text(bar.get_x() + bar.get_width() / 2.0, h + 1.0, f"{h:.1f}", ha="center", va="bottom", fontsize=5.8, rotation=90)
            elif i in (0, 1, 3):
                ax2.text(bar.get_x() + bar.get_width() / 2.0, h + 1.0, f"{h:.2f}", ha="center", va="bottom", fontsize=5.8, rotation=90)

    ax2.set_xticks(xc)
    ax2.set_xticklabels(cats)
    ax2.set_ylabel("Share of Classified Reads (%)")
    ax2.set_ylim(0, 84)
    ax2.set_title("B. 384-Plex E19 Read Partitions & Threshold Sensitivity (Fig. 3f)", loc="left", pad=8)
    ax2.grid(axis="y", linestyle="--", alpha=0.35)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.legend(loc="upper right", ncol=2, fontsize=6.6)

    plt.tight_layout(w_pad=2.0)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "fig1_multipanel_replication.pdf", bbox_inches="tight")
    fig.savefig(out_dir / "fig1_multipanel_replication.png", bbox_inches="tight")
    plt.close(fig)


def generate_fig2_dimer_and_matrix(canon: dict, dimer_tsv: Path, out_dir: Path) -> None:
    counts = []
    with open(dimer_tsv, "r", encoding="utf-8") as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            c = int(parts[3])
            if c > 0:
                counts.append(c)
    counts.sort(reverse=True)
    ranks = np.arange(1, len(counts) + 1)

    e19_dist = canon["panels"]["E19"]["pipeline_b_edlib"]["dimer_distribution"]
    top1 = e19_dist["top5_pairs"][0]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.2), gridspec_kw={"width_ratios": [1.2, 1.0]})

    ax1.plot(ranks, counts, color="#1F77B4", linewidth=1.6, marker="o", markersize=2.2, alpha=0.85)
    ax1.scatter([1], [counts[0]], color="#D62728", s=45, zorder=5)
    ax1.annotate(
        f"Rank 1: {top1['primer_1']} × {top1['primer_2']}\n"
        f"{top1['dimer_reads']:,} reads ({top1['pct_of_dimer_reads']:.2f}% of all dimers)",
        xy=(1, counts[0]),
        xytext=(5.5, counts[0] * 0.45),
        fontsize=7.5,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#FFF5F5", edgecolor="#D62728", linewidth=0.8),
        arrowprops=dict(arrowstyle="->", color="#D62728", lw=1.0),
    )
    ax1.text(
        0.04, 0.08,
        f"Candidate pairs: {e19_dist['total_candidate_pairs']:,}\n"
        f"Zero-dimer pairs: {e19_dist['zero_dimer_pairs']:,} ({e19_dist['zero_dimer_pct']:.2f}%)\n"
        f"Active pairs (≥1 read): {e19_dist['active_dimer_pairs_ge1']:,} ({e19_dist['active_dimer_pct_ge1']:.2f}%)\n"
        f"Top-10 pairs share: {e19_dist['top10_pct_of_dimer_reads']:.2f}%",
        transform=ax1.transAxes,
        fontsize=7.3,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#F8F9FA", edgecolor="#BBBBBB", linewidth=0.7),
    )
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("Active Primer-Pair Rank (log scale)")
    ax1.set_ylabel("Observed Primer-Dimer Reads (log scale)")
    ax1.set_title("A. 384-Plex E19 Heavy-Tailed Dimer Concentration", loc="left", pad=8)
    ax1.grid(True, which="both", linestyle="--", alpha=0.3)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    ps_ids = ["PS1", "PS2", "PS3"]
    pa_r_conc = [canon["panels"][p]["pipeline_a_bowtie2"]["moesm7_matrix_validation"]["conc_optimized"]["pearson_r_unordered"] for p in ps_ids]
    pb_r_conc = [canon["panels"][p]["pipeline_b_edlib"]["moesm7_matrix_validation"]["conc_optimized"]["pearson_r_unordered"] for p in ps_ids]
    pa_r_notopt = [canon["panels"][p]["pipeline_a_bowtie2"]["moesm7_matrix_validation"]["not_optimized"]["pearson_r_unordered"] for p in ps_ids]

    x = np.arange(len(ps_ids))
    w = 0.25
    b1 = ax2.bar(x - w, pa_r_conc, width=w, label="Pipeline A vs MOESM7 (conc. opt.)", color="#1F77B4", edgecolor="black", linewidth=0.5)
    b2 = ax2.bar(x, pb_r_conc, width=w, label="Pipeline B vs MOESM7 (conc. opt.)", color="#2CA02C", edgecolor="black", linewidth=0.5)
    b3 = ax2.bar(x + w, pa_r_notopt, width=w, label="Pipeline A vs MOESM7 (not opt.)", color="#BDBDBD", edgecolor="black", linewidth=0.5)

    for bars in (b1, b2, b3):
        for bar in bars:
            h = bar.get_height()
            lbl = f"{h:.5f}" if h > 0.999 else (f"{h:.4f}" if h > 0.99 else f"{h:.3f}")
            ax2.text(bar.get_x() + bar.get_width() / 2.0, h + 0.018, lbl, ha="center", va="bottom", fontsize=6.5, rotation=90 if h > 0.99 else 0)

    ax2.set_xticks(x)
    ax2.set_xticklabels(ps_ids)
    ax2.set_ylabel("Pairwise Matrix Pearson r (18,528 unordered pairs)")
    ax2.set_ylim(0, 1.48)
    ax2.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax2.set_title("B. 96-Plex Dimer Matrix Validation vs Published MOESM7", loc="left", pad=8)
    ax2.grid(axis="y", linestyle="--", alpha=0.35)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.legend(loc="upper center", ncol=1, fontsize=7.0)

    plt.tight_layout(w_pad=2.0)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "fig2_dimer_concentration_and_matrix.pdf", bbox_inches="tight")
    fig.savefig(out_dir / "fig2_dimer_concentration_and_matrix.png", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate data-driven replication figures.")
    parser.add_argument("--data-dir", type=str, default=None)
    parser.add_argument("--out-dir", type=str, default=None)
    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parent.parent
    data_dir = Path(args.data_dir).resolve() if args.data_dir else (base_dir / "data")
    fig_dir = Path(args.out_dir).resolve() if args.out_dir else (base_dir / "article" / "figures")

    canon = json.loads((data_dir / "canonical_results.json").read_text(encoding="utf-8"))
    dimer_tsv = data_dir / "saddle_384_dimer_pairs.tsv"

    setup_style()
    generate_fig1_multipanel(canon, fig_dir)
    generate_fig2_dimer_and_matrix(canon, dimer_tsv, fig_dir)
    print(f"Generated data-driven figures in {fig_dir}")


if __name__ == "__main__":
    main()
