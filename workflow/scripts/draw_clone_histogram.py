import argparse
import sys
import os
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd


def load_palette(path):
    with open(path) as fh:
        colors = [ln.strip() for ln in fh if ln.strip()]
    if not colors:
        sys.exit(f"Missing palette.")
    return colors[1:]


def load_data(path):
    phy_df = pd.read_csv(path, sep="\t")

    required = {"mutation_id", "clone_id", "sample_id", "clonal_prev"}
    missing = required - set(phy_df.columns)
    if missing:
        sys.exit(f"Missing mandatory columns: {', '.join(sorted(missing))}")

    phy_df = phy_df[phy_df["clone_id"] != -1].copy()

    counts = (
        phy_df.groupby(["sample_id", "clone_id"])["mutation_id"]
        .nunique()
        .reset_index(name="n_mutations")
    )

    totals = (
        counts.groupby("sample_id")["n_mutations"]
        .sum()
        .reset_index(name="sample_total_mutations")
    )

    counts = counts.merge(totals, on="sample_id", how="left")
    counts["pct_mutations"] = (
        counts["n_mutations"] / counts["sample_total_mutations"] * 100
    )

    cp = (
        phy_df.groupby(["sample_id", "clone_id"])["clonal_prev"]
        .max()
        .reset_index()
    )

    df = counts.merge(cp, on=["sample_id", "clone_id"], how="left")

    if df["clonal_prev"].dropna().max() <= 1.0:
        df["clonal_prev_pct"] = df["clonal_prev"] * 100
    else:
        df["clonal_prev_pct"] = df["clonal_prev"]

    df = df.sort_values(["sample_id", "clone_id"]).reset_index(drop=True)
    return df


def sanitize_filename(text):
    text = str(text)
    text = re.sub(r"[^\w.-]+", "_", text)
    return text.strip("_")


def plot_one_sample(sample_df, sample_id, clone_colors, output_file):
    """
    Plot clonal histogram,  one per sample.
    Args:
        sample_df (str): Pandas dataframe with sample_id, clones, colors and clonal prevalence.
        sample_id (str): sample_id for the plot title.
        clone colors (str): Dictionary with the colors asigned to clones (Dict).
        out_file (str): Path to output file (PNG).
        
    """
    BG_COLOR = "white"
    GRID_COLOR = "#e8e8e8"
    FONT = "DejaVu Sans"

    sample_df = sample_df.sort_values("clone_id").reset_index(drop=True).copy()
    sample_df["clone_color"] = [
        clone_colors[int(cid) % len(clone_colors)] for cid in sample_df["clone_id"]
    ]

    n = len(sample_df)
    fig_h = max(3.0, 0.65 * n + 1.6)
    scale = 0.6
    fig, ax = plt.subplots(figsize=(8.5 * scale, fig_h * scale), dpi=300, facecolor=BG_COLOR)
    ax.set_facecolor(BG_COLOR)

    y = np.arange(n)
    bar_h = 0.34

    bars_prev = ax.barh(
        y - bar_h / 2,
        sample_df["clonal_prev_pct"],
        height=bar_h,
        color=sample_df["clone_color"],
        edgecolor="none",
        zorder=3,
        label="Clonal prevalence"
    )

    bars_mut = ax.barh(
        y + bar_h / 2,
        sample_df["pct_mutations"],
        height=bar_h,
        color="white",
        edgecolor=sample_df["clone_color"],
        linewidth=1.6,
        zorder=3,
        label="% mutations"
    )

    xmax = 75
    
    ax.set_xlim(0, xmax)
    ax.xaxis.set_major_locator(mticker.MultipleLocator(10))
    ax.xaxis.set_major_formatter(mticker.PercentFormatter(xmax=100, decimals=0))
    ax.set_xlabel("Percentage", fontsize=8, color="#555", fontfamily=FONT)

    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"Clone {int(c)}" for c in sample_df["clone_id"]],
        fontsize=7,
        fontfamily=FONT
    )
    ax.invert_yaxis()

    for tick, c in zip(ax.get_yticklabels(), sample_df["clone_color"]):
        tick.set_color(c)
        tick.set_fontweight("bold")

    ax.tick_params(axis="x", labelsize=7, colors="#777")
    ax.tick_params(axis="y", length=0, pad=6)

    ax.xaxis.grid(True, color=GRID_COLOR, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)

    for bar, val, c in zip(bars_prev, sample_df["clonal_prev_pct"], sample_df["clone_color"]):
        ax.text(
            bar.get_width() + xmax * 0.01,
            bar.get_y() + bar.get_height() / 2,
            f"{val:.1f}%",
            ha="left",
            va="center",
            fontsize=6,
            color=c,
            fontfamily=FONT
        )

    for bar, val, nmut, c in zip(bars_mut, sample_df["pct_mutations"], sample_df["n_mutations"], sample_df["clone_color"]):
        ax.text(
            bar.get_width() + xmax * 0.01,
            bar.get_y() + bar.get_height() / 2,
            f"{val:.1f}% #mut={int(nmut)}",
            ha="left",
            va="center",
            fontsize=6,
            color=c,
            fontfamily=FONT
        )

    ax.set_title(
        f"Sample {sample_id}",
        fontsize=9,
        color="#333",
        fontfamily=FONT,
        pad=16
    )

    for sp in ["top", "right", "left"]:
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(GRID_COLOR)

    ax.legend(
        loc="upper left",
        bbox_to_anchor=(0.0, 1.05),
        ncol=2,
        frameon=False,
        fontsize=6,
        borderaxespad=0.1
    )

    plt.tight_layout(pad=1.2)
    fig.savefig(output_file, bbox_inches="tight", facecolor=BG_COLOR)
    plt.close(fig)


def plot(df, clone_colors, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for sample_id, sample_df in df.groupby("sample_id", sort=False):
        safe_sample = sanitize_filename(sample_id)
        output_file = os.path.join(output_dir, f"{safe_sample}_mutation_contribution.png")
        plot_one_sample(sample_df, sample_id, clone_colors, output_file)


if __name__ == "__main__":
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--tree_df", action="store", required=True)
    input_parser.add_argument("--palette", action="store", required=True)
    input_parser.add_argument("--out_dir", action="store", required=True)
    args = input_parser.parse_args()

    clone_colors = load_palette(args.palette)
    counts       = load_data(args.tree_df)
    plot(counts, clone_colors, args.out_dir)