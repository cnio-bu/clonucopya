import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd


def load_palette(path):
    # Skip first color (parent cell)
    with open(path) as fh:
        colors = [ln.strip() for ln in fh if ln.strip()]
    if not colors:
        sys.exit(f"[ERROR] La paleta '{path}' está vacía.")
    return colors[1:]


def load_data(path):
    df = pd.read_csv(path, sep="\t")

    counts = (
        df.groupby("clone_id")["mutation_id"]
          .count().reset_index()
          .rename(columns={"mutation_id": "n_alterations"})
    )
    meta = df.groupby("clone_id")[["ccf", "clonal_prev"]].mean().reset_index()
    counts = counts.merge(meta, on="clone_id")
    counts["pct"] = counts["n_alterations"] / counts["n_alterations"].sum() * 100
    counts = counts.sort_values("clone_id").reset_index(drop=True)
    return counts


def plot(counts, clone_colors, output):
    BG_COLOR   = "#fafafa"
    GRID_COLOR = "#e8e8e8"
    FONT       = "DejaVu Sans"

    # ~half A4 page width at 300 dpi — fits cleanly in a PDF report
    FIGSIZE = (8.0, 5.0)
    DPI     = 300

    n      = len(counts)
    colors = [clone_colors[i % len(clone_colors)] for i in range(n)]
    x      = np.arange(n)
    max_v  = counts["n_alterations"].max()

    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI, facecolor=BG_COLOR)
    ax.set_facecolor(BG_COLOR)

    bars = ax.bar(x, counts["n_alterations"], width=0.52,
                  color=colors, zorder=3, linewidth=0)

    # Labels above bars
    for bar, row, color in zip(bars, counts.itertuples(), colors):
        yv = bar.get_height()
        bx = bar.get_x() + bar.get_width() / 2

        ax.text(bx, yv + max_v * 0.014,
                f"{int(yv):,}".replace(",", "\u202f"),
                ha="center", va="bottom",
                fontsize=11, fontweight="bold",
                color=color, fontfamily=FONT)

        pct_str = f"{row.pct:.1f}%" if row.pct >= 0.1 else "<0.1%"
        ax.text(bx, yv + max_v * 0.014 + max_v * 0.052,
                pct_str,
                ha="center", va="bottom",
                fontsize=8.5, color="#999", fontfamily=FONT)

    # X axis
    xlabels = [
        f"Clone {int(r.clone_id)}\nCCF {r.ccf:.2f} · prev {r.clonal_prev:.0%}"
        for r in counts.itertuples()
    ]
    ax.set_xticks(x)
    ax.set_xticklabels(xlabels, fontsize=9, fontfamily=FONT,
                       color="#444", multialignment="center")
    ax.tick_params(axis="x", length=0, pad=6)

    # Y axis
    ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins=7))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(
        lambda v, _: f"{int(v/1000)}k" if v >= 1000 else str(int(v))
    ))
    ax.set_ylabel("Number of alterations", fontsize=9.5, color="#555", fontfamily=FONT)
    ax.tick_params(axis="y", labelsize=9, colors="#777")
    ax.set_ylim(0, max_v * 1.30)

    # Grid and spines
    ax.yaxis.grid(True, color=GRID_COLOR, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for sp in ["top", "right", "left"]:
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(GRID_COLOR)

    plt.tight_layout(pad=1.5)
    fig.savefig(output, dpi=DPI, bbox_inches="tight", facecolor=BG_COLOR)
    plt.close()


if __name__ == "__main__":
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--input",   "-i", action="store", required=True)
    input_parser.add_argument("--palette", "-p", action="store", required=True)
    input_parser.add_argument("--output",  "-o", action="store", default="clone_alterations.png")
    args = input_parser.parse_args()

    clone_colors = load_palette(args.palette)
    counts       = load_data(args.input)
    plot(counts, clone_colors, args.output)
