#!/usr/bin/env python3
"""Generate the publication-style BC-HI inter-subunit interface summary figure.

The numerical values are precomputed measurements used for the manuscript
figure; this script does not calculate interface properties from coordinates.
"""
from __future__ import annotations
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.0,
    "figure.dpi": 600,
})

# Native structure: 4FTB / ugf2
native = [
    (1717 + 1669) + (1722 + 1666),  # interface area (Å²)
    4 + 4,                           # salt bridges
    10 + 10,                         # hydrogen bonds
    184 + 184,                       # non-bonded contacts
]

# Early-disassembly intermediate: 9LZL / uge6
intermediate = [
    (1652 + 1649) + (1649 + 1638),  # interface area (Å²)
    3 + 3,                           # salt bridges
    5 + 5,                           # hydrogen bonds
    87 + 80,                         # non-bonded contacts
]

relative_reduction = [
    (n - i) / n * 100.0 for n, i in zip(native, intermediate)
]


def create_figure(output: str = "Nature_Style_BCHI_Summed_BW.png") -> None:
    """Create and save the BC-HI interface summary figure."""
    fig = plt.figure(figsize=(180 / 25.4, 120 / 25.4))
    gs = fig.add_gridspec(2, 3, wspace=0.35, hspace=0.45)

    panel_labels = ["a", "b", "c", "d", "e"]
    titles = [
        r"Total interface area ($\AA^2$)",
        "Total salt bridges",
        "Total hydrogen bonds",
        "Total non-bonded contacts",
        "Relative reduction (%)",
    ]

    for i in range(5):
        row, col = divmod(i, 3)
        ax = fig.add_subplot(gs[row, col])

        if i < 4:
            bars = ax.bar(
                ["Native", "Intermediate"],
                [native[i], intermediate[i]],
                width=0.6,
                edgecolor="black",
                linewidth=0.8,
            )
            bars[1].set_hatch("////")
        else:
            bars = ax.bar(
                ["Area", "SB", "HB", "NBC"],
                relative_reduction,
                width=0.6,
                edgecolor="black",
                linewidth=0.8,
                hatch="////",
            )
            ax.set_ylim(0, 75)
            for bar in bars:
                h = bar.get_height()
                ax.text(bar.get_x() + bar.get_width() / 2, h + 2,
                        f"{h:.1f}%", ha="center", fontsize=6, fontweight="bold")

        ax.set_title(titles[i], fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="x", length=0)
        for tick in ax.get_xticklabels() + ax.get_yticklabels():
            tick.set_fontweight("bold")
        ax.text(-0.25, 1.1, panel_labels[i], transform=ax.transAxes,
                fontsize=10, fontweight="bold", va="top", ha="right")

    axes = fig.axes
    axes[5].axis("off")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved to {output}")


if __name__ == "__main__":
    create_figure()
