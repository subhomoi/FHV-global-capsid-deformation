#!/usr/bin/env python3
"""
Plot geometrical deformation of the FHV 2-fold subassembly during early disassembly.

Input:
    FHV_Global_Deformation_per_residue.csv

Output folder:
    FHV_Global_Geometrical_Deformation_Figure/

Outputs:
    FHV_Global_Geometrical_Deformation_Nature.png
    FHV_Global_Geometrical_Deformation_Nature.pdf
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# USER SETTINGS
# ============================================================

INPUT_DIR = Path(".")

INPUT_CSV = INPUT_DIR / "FHV_Global_Deformation_per_residue.csv"

OUTPUT_DIR = INPUT_DIR / "FHV_Global_Geometrical_Deformation_Figure"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_PNG = OUTPUT_DIR / "FHV_Global_Geometrical_Deformation_Nature.png"
OUTPUT_PDF = OUTPUT_DIR / "FHV_Global_Geometrical_Deformation_Nature.pdf"


# ============================================================
# NATURE-STYLE PLOT SETTINGS
# ============================================================

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9,
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 8,
        "axes.linewidth": 1.0,
        "xtick.major.width": 1.0,
        "ytick.major.width": 1.0,
        "xtick.major.size": 4,
        "ytick.major.size": 4,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.dpi": 600,
    }
)


# ============================================================
# LOAD DATA
# ============================================================

if not INPUT_CSV.exists():
    raise FileNotFoundError(
        f"Input CSV was not found:\n{INPUT_CSV.resolve()}"
    )

df = pd.read_csv(INPUT_CSV)

required_columns = [
    "native_radius_A",
    "target_radius_A",
    "angular_change_deg",
]

missing_columns = [
    column for column in required_columns if column not in df.columns
]

if missing_columns:
    raise ValueError(
        "The following required columns are missing from the CSV:\n"
        + "\n".join(f"  - {column}" for column in missing_columns)
        + "\n\nAvailable columns are:\n"
        + "\n".join(f"  - {column}" for column in df.columns)
    )


plot_df = df[required_columns].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True)

if plot_df.empty:
    raise ValueError("No complete numerical observations remain after filtering the required columns.")

native_radius = plot_df["native_radius_A"].to_numpy(dtype=float)
target_radius = plot_df["target_radius_A"].to_numpy(dtype=float)
angular_change = plot_df["angular_change_deg"].to_numpy(dtype=float)


# ============================================================
# CALCULATE SUMMARY STATISTICS
# ============================================================

native_mean_radius = np.mean(native_radius)
target_mean_radius = np.mean(target_radius)
mean_radius_change = target_mean_radius - native_mean_radius

angular_mean = np.mean(angular_change)
angular_median = np.median(angular_change)
angular_max = np.max(angular_change)

radial_change = target_radius - native_radius
fraction_outward = np.mean(radial_change > 0) * 100


# ============================================================
# CREATE FIGURE
# ============================================================

fig, axes = plt.subplots(
    nrows=1,
    ncols=2,
    figsize=(8.5, 4.5), # Adjusted width for 2 panels
)

ax_a, ax_b = axes


# ============================================================
# PANEL A — RADIAL EXPANSION
# ============================================================

violin_data = [
    native_radius,
    target_radius,
]

violin = ax_a.violinplot(
    violin_data,
    positions=[1, 2],
    widths=0.72,
    showmeans=False,
    showmedians=False,
    showextrema=True,
)

for body in violin["bodies"]:
    body.set_facecolor("#9ecae1")
    body.set_edgecolor("#6b6b6b")
    body.set_alpha(0.78)
    body.set_linewidth(0.9)

for key in ["cbars", "cmins", "cmaxes"]:
    violin[key].set_color("#6b6b6b")
    violin[key].set_linewidth(0.9)

# Mean values
ax_a.scatter(
    [1, 2],
    [native_mean_radius, target_mean_radius],
    s=85,
    color="#1f77b4",
    edgecolor="white",
    linewidth=0.8,
    zorder=5,
)

# Mean connector
ax_a.plot(
    [1, 2],
    [native_mean_radius, target_mean_radius],
    color="#ff7f0e",
    linestyle="--",
    linewidth=1.4,
    zorder=4,
)

# Mean-radius labels
ax_a.text(
    1,
    native_mean_radius + 2.0,
    f"{native_mean_radius:.2f} Å",
    ha="center",
    va="bottom",
    fontsize=9,
)

ax_a.text(
    2,
    target_mean_radius + 2.0,
    f"{target_mean_radius:.2f} Å",
    ha="center",
    va="bottom",
    fontsize=9,
)

# Delta annotation placed above the mean labels
ax_a.text(
    1.5,
    max(native_mean_radius, target_mean_radius) + 5.0,
    f"Δ = {mean_radius_change:+.2f} Å",
    ha="center",
    va="center",
    fontsize=10,
    fontweight="bold",
)

ax_a.set_title(
    "Radial expansion",
    loc="left",
    fontweight="bold",
    pad=10,
)

# Updated Y-axis label
ax_a.set_ylabel("Radius from 2-fold center of mass (Å)")

ax_a.set_xticks([1, 2])
ax_a.set_xticklabels(
    [
        "Native",
        "Early intermediate",
    ]
)

ax_a.set_xlim(0.55, 2.45)
ax_a.set_ylim(0, max(target_radius.max(), native_radius.max()) + 8)

ax_a.spines["top"].set_visible(False)
ax_a.spines["right"].set_visible(False)


# ============================================================
# PANEL B — ANGULAR DEFORMATION
# ============================================================

ax_b.hist(
    angular_change,
    bins=30,
    color="#6baed6",
    edgecolor="#3f3f3f",
    linewidth=0.6,
)

ax_b.axvline(
    angular_mean,
    color="#1f77b4",
    linestyle="--",
    linewidth=1.5,
)

ax_b.text(
    0.98,
    0.96,
    (
        f"Mean     {angular_mean:.2f}°\n"
        f"Median   {angular_median:.2f}°\n"
        f"Maximum  {angular_max:.2f}°"
    ),
    transform=ax_b.transAxes,
    ha="right",
    va="top",
    fontsize=9,
    linespacing=1.5,
)

ax_b.set_title(
    "Angular deformation",
    loc="left",
    fontweight="bold",
    pad=10,
)

ax_b.set_xlabel("Angular change (°)")
ax_b.set_ylabel("Frequency")

if np.all(angular_change >= 0):
    ax_b.set_xlim(left=0)

ax_b.spines["top"].set_visible(False)
ax_b.spines["right"].set_visible(False)


# ============================================================
# PANEL LABELS
# ============================================================

for axis, label in zip(axes, ["A.", "B."]):
    axis.text(
        -0.16,
        1.08,
        label,
        transform=axis.transAxes,
        fontsize=14,
        fontweight="bold",
        va="top",
        ha="left",
    )


# ============================================================
# MAIN FIGURE TITLE
# ============================================================

# Updated Main Title
fig.suptitle(
    "Geometrical deformation of the FHV 2-fold subassembly during early disassembly",
    fontsize=16,
    fontweight="bold",
    y=0.995,
)

# Updated Subtitle
fig.text(
    0.5,
    0.952,
    "Native FHV 2-fold axis (PDB 4FTB) versus the early-disassembly intermediate",
    ha="center",
    va="center",
    fontsize=10,
)


# ============================================================
# FINAL LAYOUT AND SAVE
# ============================================================

plt.tight_layout(
    rect=[0.0, 0.0, 1.0, 0.90]
)

fig.savefig(
    OUTPUT_PNG,
    dpi=600,
    bbox_inches="tight",
    facecolor="white",
)

fig.savefig(
    OUTPUT_PDF,
    bbox_inches="tight",
    facecolor="white",
)

plt.show()

print("\nFigure generation completed successfully.")
print(f"Input file:       {INPUT_CSV.resolve()}")
print(f"Output directory: {OUTPUT_DIR.resolve()}")
print(f"PNG saved to:     {OUTPUT_PNG.resolve()}")
print(f"PDF saved to:     {OUTPUT_PDF.resolve()}")