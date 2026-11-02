"""Performance evaluation figure for the GWI Pinpoint experiments.

Precision of the base model (Claude Opus 4.8) versus Pinpoint on two tasks
(network analysis and fine-mapping), each repeated three times with
different backbone models. Scores are pooled per condition and summarised
as boxplots in two side-by-side panels; means are annotated above each
box. Scores are reported as percentages.

Outputs a raster figure at the requested path plus a vector PDF twin
(recommended for LaTeX inclusion: small file, crisp at any scale, fast
compilation).
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np

# Publishable TrueType fonts embedded in the PDF (editable text, small files)
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 0.8,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

DPI = 300  # raster output resolution (use the PDF twin for LaTeX)

C_BASE = "#3B6BA5"   # base model
C_PIN = "#D55C00"    # Pinpoint (matches the network figure)
C_MED = "#303030"    # median line
BOX_ALPHA = 0.30     # boxplot fill transparency

# --------------------------------------------------------------------------- #
# Precision scores from three repeats (Claude Opus 4.8, GLM 5.3 and         #
# DeepSeek 4.1). Scores are pooled per condition: they are deliberately      #
# not mapped to the models on the figure.                                    #
# --------------------------------------------------------------------------- #
NET_BASE = np.array([0.5469, 0.797, 0.7656])
NET_PIN = np.array([0.7431, 0.8947, 0.7857])
FIN_BASE = np.array([0.7157, 0.45, 0.4210])
FIN_PIN = np.array([0.8515, 0.7777, 0.80])

# scores are displayed as percentages
NET = (NET_BASE * 100, NET_PIN * 100)
FIN = (FIN_BASE * 100, FIN_PIN * 100)

TASKS = [
    ("Network analysis", NET),
    ("Fine-mapping", FIN),
]


def save_figure(fig, path, dpi=DPI):
    """Save to `path`; if not already a PDF, also write a vector PDF twin."""
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    root, ext = os.path.splitext(path)
    if ext.lower() != ".pdf":
        fig.savefig(root + ".pdf", bbox_inches="tight", facecolor="white")


def style_boxplot(bp, colors):
    """Colour boxplot elements per condition: first box uses colors[0]
    (base model), the second colors[1] (Pinpoint)."""
    lw = 1.1
    for box, col in zip(bp["boxes"], colors):
        box.set(facecolor=col, alpha=BOX_ALPHA, edgecolor=col, linewidth=lw)
    for med in bp["medians"]:
        med.set(color=C_MED, linewidth=1.6)
    for whisk, col in zip(bp["whiskers"], colors * 2):
        whisk.set(color=col, linewidth=lw)
    for cap, col in zip(bp["caps"], colors * 2):
        cap.set(color=col, linewidth=lw)


def plot_performance(tasks=TASKS):
    """Render the two-panel boxplot figure."""
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), sharey=True)
    for ax, (title, (base, pin)) in zip(axes, tasks):
        data = [base, pin]
        bp = ax.boxplot(
            data, positions=[0, 1], widths=0.42, patch_artist=True,
            showfliers=False, whiskerprops=dict(linewidth=1.1),
            capprops=dict(linewidth=1.1), medianprops=dict(linewidth=1.6),
        )
        style_boxplot(bp, [C_BASE, C_PIN])

        # mean value above each box
        for x, vals, col in zip([0, 1], data, [C_BASE, C_PIN]):
            ax.text(x, max(vals) + 3.5, f"mean = {vals.mean():.0f}%",
                    ha="center", fontsize=7.5, color=col)

        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Base model", "Pinpoint"], fontsize=8)
        ax.set_xlim(-0.55, 1.55)
        ax.set_ylim(30, 105)
        ax.set_yticks([30, 50, 70, 90])
        ax.tick_params(axis="both", labelsize=7.5)
        ax.set_ylabel("Precision (%)", fontsize=8.5)
        ax.set_title(title, fontsize=9, fontweight="bold", pad=8)

    fig.subplots_adjust(wspace=0.18, top=0.88, bottom=0.12)
    return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-path", help="Path to save the performance figure at",
        default="performance.png",
    )
    args = parser.parse_args()

    fig = plot_performance()
    save_figure(fig, args.output_path)
    plt.show()
