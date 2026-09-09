import argparse
import sys
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_palette(path):
    with open(path) as fh:
        colors = [ln.strip() for ln in fh if ln.strip()]
    if not colors:
        sys.exit("Missing palette.")
    return colors[1:]


def load_data(path):
    """
    Aggregate clonal data across the whole study (all samples combined),
    instead of per-sample.
    """
    phy_df = pd.read_csv(path, sep="\t")

    required = {"mutation_id", "clone_id", "sample_id", "clonal_prev"}
    missing = required - set(phy_df.columns)
    if missing:
        sys.exit("Missing mandatory columns: " + ", ".join(sorted(missing)))

    phy_df = phy_df[phy_df["clone_id"] != -1].copy()

    counts = (
        phy_df.groupby("clone_id")["mutation_id"]
        .nunique()
        .reset_index(name="n_mutations")
    )

    total_mut = counts["n_mutations"].sum()
    counts["pct_mutations"] = counts["n_mutations"] / total_mut * 100

    cp = (
        phy_df.groupby("clone_id")["clonal_prev"]
        .mean()
        .reset_index()
    )

    df = counts.merge(cp, on="clone_id", how="left")

    if df["clonal_prev"].dropna().max() <= 1.0:
        df["clonal_prev_pct"] = df["clonal_prev"] * 100
    else:
        df["clonal_prev_pct"] = df["clonal_prev"]

    # Clonal prevalence normalization to adjust width of bars 360 grades
    total_prev = df["clonal_prev_pct"].sum()
    df["prev_frac_norm"] = df["clonal_prev_pct"] / total_prev

    df = df.sort_values("clone_id").reset_index(drop=True)
    return df


def plot_study(df, clone_colors, output_file):
    BG_COLOR = "white"
    FONT = "DejaVu Sans"

    df = df.sort_values("clone_id").reset_index(drop=True).copy()
    df["clone_color"] = [
        clone_colors[int(cid) % len(clone_colors)] for cid in df["clone_id"]
    ]

    widths = df["prev_frac_norm"].to_numpy(dtype=float) * 2 * np.pi

    # Set internal radius of the plot and width of bars
    r_inner = 0.35
    r_max_extra = 1.0
    # max_mut = df["n_mutations"].max()
    # heights = r_inner + (df["n_mutations"].to_numpy(dtype=float) / max_mut) * r_max_extra
    min_bar_height = 0.025

    nmut = df["n_mutations"].to_numpy(dtype=float)
    max_mut = nmut.max()
    
    heights = (nmut / max_mut) * r_max_extra
    heights = np.maximum(heights, min_bar_height)

    
    fig, ax = plt.subplots(
        figsize=(6.2, 5.4), dpi=300, facecolor=BG_COLOR,
        subplot_kw=dict(projection="polar"),
    )
    ax.set_facecolor(BG_COLOR)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    # Set angle of the stacked bars
    starts = np.concatenate(([0.0], np.cumsum(widths)[:-1]))
    mids = starts + widths / 2.0

    ax.bar(
        x=mids,
        height=heights,
        width=widths,
        bottom=r_inner,
        color=df["clone_color"],
        edgecolor="white",
        linewidth=1.6,
        align="center",
    )

    # Create labelling
    for i, (cid, nmut, mid, h, color) in enumerate(
        zip(df["clone_id"], df["n_mutations"], mids, heights, df["clone_color"])
    ):
        prev_value = df.loc[i, "clonal_prev_pct"] / 100.0

        mid_deg = 90 - np.degrees(mid)
        rot = mid_deg % 360
        if 90 < rot < 270:
            rot += 180
        rot = rot % 360

        # height of the bar
        r_outer = r_inner + h

        # Clonal prevalence label
        r_label_in = r_inner + h * 0.45

        ax.text(
            mid, r_label_in,
            f"c{int(cid)}\n{prev_value:.2f}",
            ha="center", va="center", linespacing=1.1,
            fontsize=5, color="black", fontfamily=FONT, fontweight="bold",
            bbox=dict(
                boxstyle="round,pad=0.20",
                facecolor="white",
                edgecolor=color,
                linewidth=1.1,
            ),
        )

        # Mutation number label
        r_label_out = r_outer + 0.26

        ax.text(
            mid, r_label_out,
            f"{int(nmut)} mut",
            ha="center", va="center",
            fontsize=8, color="#333333", fontfamily=FONT, fontweight="bold",
            rotation=rot, rotation_mode="anchor",
        )

    # ax.set_ylim(0, r_inner + r_max_extra + 0.45)
    label_gap = 0.26
    extra_margin = 0.08
    ax.set_ylim(0, r_inner + r_max_extra + label_gap + extra_margin)
    
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["polar"].set_visible(False)
    ax.grid(False)


    # plt.tight_layout(pad=1.2)
    # fig.savefig(output_file, bbox_inches="tight", facecolor=BG_COLOR)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

    fig.savefig(
        output_file,
        bbox_inches="tight",
        pad_inches=0.01,
        facecolor=BG_COLOR,
    )
    plt.close(fig)


def plot(study, df, clone_colors, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{study}_clonal_composition.png")
    plot_study(df, clone_colors, output_file)

if __name__ == "__main__":
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--tree_df", action="store", required=True)
    input_parser.add_argument("--study", action="store", required=True)
    input_parser.add_argument("--palette", action="store", required=True)
    input_parser.add_argument("--out_dir", action="store", required=True)
    args = input_parser.parse_args()

    clone_colors = load_palette(args.palette)
    counts       = load_data(args.tree_df)
    plot(args.study, counts, clone_colors, args.out_dir)
