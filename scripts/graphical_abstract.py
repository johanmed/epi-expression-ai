"""Graphical abstract for the GWI Pinpoint paper (Cell-style, methodology focus).

Layout: title -> 4 columns (1 GWI model & datasets, 2 Pinpoint agentic LLM,
3 genomic tasks incl. PheWAS, 4 disease insight cascade).

Outputs a raster figure at the requested path plus a vector PDF twin
(recommended for LaTeX inclusion: small file, crisp at any scale, fast
compilation).
"""

import argparse
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Ellipse

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

INK = "#1a1a1a"
C_ORANGE = "#D55C00"
C_BLUE = "#3B6BA5"
C_TEAL = "#008B8B"
C_GREEN = "#4A9C3D"
C_PURPLE = "#7B4FA3"
C_BROWN = "#8B6914"
C_GREY = "#6B7280"
C_GREY_L = "#9AA0A6"
C_LIGHT = "#F2F4F7"

W, H = 140, 92

DPI = 300  # raster output resolution (use the PDF twin for LaTeX)


def save_figure(fig, path, dpi=DPI):
    """Save to `path`; if not already a PDF, also write a vector PDF twin."""
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    root, ext = os.path.splitext(path)
    if ext.lower() != ".pdf":
        fig.savefig(root + ".pdf", bbox_inches="tight", facecolor="white")


def rbox(ax, x0, y0, x1, y1, fc="white", ec="#5B6470", lw=1.0, r=1.6, z=2,
         alpha=1.0, ls="-"):
    p = FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                       boxstyle=f"round,pad=0,rounding_size={r}",
                       fc=fc, ec=ec, lw=lw, zorder=z, alpha=alpha, linestyle=ls)
    ax.add_patch(p)
    return p


def arrow(ax, x0, y0, x1, y1, color="#5B6470", lw=1.6, ms=13, style="-|>", z=3,
          ls="-", rad=0.0):
    cs = f"arc3,rad={rad}" if rad else None
    a = FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style, color=color,
                        lw=lw, mutation_scale=ms, zorder=z, linestyle=ls,
                        connectionstyle=cs, shrinkA=2, shrinkB=2)
    ax.add_patch(a)


def header(ax, x0, x1, y0, text, color):
    rbox(ax, x0, y0, x1, y0 + 4.4, fc=color, ec=color, r=1.2, z=4)
    ax.text((x0 + x1) / 2, y0 + 2.2, text, ha="center", va="center",
            fontsize=8.5, fontweight="bold", color="white", zorder=5)


def mouse_icon(ax, cx, cy, s=1.0, color="#8A93A3"):
    ax.add_patch(Ellipse((cx, cy), 3.4 * s, 2.0 * s, fc=color, ec="none", zorder=6))
    ax.add_patch(Circle((cx + 1.9 * s, cy + 0.55 * s), 0.95 * s, fc=color,
                        ec="none", zorder=6))
    ax.add_patch(Circle((cx + 1.6 * s, cy + 1.35 * s), 0.42 * s, fc=color,
                        ec="none", zorder=6))
    ax.add_patch(Circle((cx + 2.35 * s, cy + 0.75 * s), 0.09 * s, fc="white",
                        ec="none", zorder=7))
    ax.add_patch(Circle((cx + 2.8 * s, cy + 0.42 * s), 0.12 * s, fc="#41464e",
                        ec="none", zorder=7))
    tail = FancyArrowPatch((cx - 1.7 * s, cy - 0.2 * s), (cx - 3.2 * s, cy + 1.1 * s),
                           arrowstyle="-", color=color, lw=1.4 * s,
                           connectionstyle="arc3,rad=0.5", zorder=6)
    ax.add_patch(tail)


def helix_icon(ax, cx, cy, s=1.0, color=C_BLUE):
    t = [i / 20 for i in range(21)]
    xs = [cx + (ti - 0.5) * 4.4 * s for ti in t]
    y1 = [cy + 0.9 * s * math.sin(2 * math.pi * ti * 1.5) for ti in t]
    y2 = [cy - 0.9 * s * math.sin(2 * math.pi * ti * 1.5) for ti in t]
    ax.plot(xs, y1, color=color, lw=1.3, zorder=6)
    ax.plot(xs, y2, color=color, lw=1.3, zorder=6)
    for i in range(0, 21, 4):
        ax.plot([xs[i], xs[i]], [y1[i], y2[i]], color=color, lw=0.9,
                alpha=0.7, zorder=6)


def methylation_icon(ax, cx, cy, s=1.0, color=C_TEAL):
    ax.plot([cx - 2.2 * s, cx + 2.2 * s], [cy, cy], color=color, lw=1.6, zorder=6)
    for dx in (-1.3, 0.0, 1.3):
        ax.plot([cx + dx * s, cx + dx * s], [cy, cy + 0.85 * s], color=color,
                lw=1.1, zorder=6)
        ax.add_patch(Circle((cx + dx * s, cy + 1.05 * s), 0.38 * s, fc=color,
                            ec="none", zorder=6))


def hourglass_icon(ax, cx, cy, s=1.0, color=C_PURPLE):
    ax.plot([cx - 1.6 * s, cx + 1.6 * s], [cy + 1.0 * s, cy - 1.0 * s],
            color=color, lw=1.4, zorder=6)
    ax.plot([cx - 1.6 * s, cx + 1.6 * s], [cy - 1.0 * s, cy + 1.0 * s],
            color=color, lw=1.4, zorder=6)
    ax.plot([cx - 1.6 * s, cx + 1.6 * s], [cy + 1.0 * s, cy + 1.0 * s],
            color=color, lw=1.4, zorder=6)
    ax.plot([cx - 1.6 * s, cx + 1.6 * s], [cy - 1.0 * s, cy - 1.0 * s],
            color=color, lw=1.4, zorder=6)


def db_icon(ax, cx, cy, s=1.0, color="#5B6470"):
    w, h = 1.8 * s, 2.4 * s
    ax.add_patch(FancyBboxPatch((cx - w / 2, cy - h / 2), w, h,
                                boxstyle="round,pad=0,rounding_size=0.3",
                                fc="white", ec=color, lw=1.1, zorder=6))
    ax.plot([cx - w / 2, cx + w / 2], [cy, cy], color=color, lw=1.0, zorder=6)


def magnifier_icon(ax, cx, cy, s=1.0, color=C_ORANGE):
    ax.add_patch(Circle((cx - 0.4 * s, cy + 0.4 * s), 1.05 * s, fc="none",
                        ec=color, lw=1.6, zorder=6))
    ax.plot([cx + 0.4 * s, cx + 1.6 * s], [cy - 0.4 * s, cy - 1.6 * s],
            color=color, lw=2.0, solid_capstyle="round", zorder=6)


def network_icon(ax, cx, cy, s=1.0):
    nodes = {
        "Tnf": (cx - 0.4 * s, cy + 0.6 * s, C_ORANGE),
        "Il6": (cx + 1.1 * s, cy + 0.9 * s, C_ORANGE),
        "Il1b": (cx + 1.6 * s, cy - 0.3 * s, C_ORANGE),
        "Nlrp3": (cx + 0.6 * s, cy + 1.9 * s, C_ORANGE),
        "Bche": (cx - 2.4 * s, cy + 1.4 * s, C_TEAL),
        "Mbp": (cx - 2.2 * s, cy - 1.2 * s, C_PURPLE),
        "Mog": (cx - 0.9 * s, cy - 1.9 * s, C_PURPLE),
        "Nqo1": (cx + 2.6 * s, cy - 1.4 * s, C_GREEN),
        "Pon1": (cx + 1.7 * s, cy - 2.3 * s, C_GREEN),
    }
    edges = [("Tnf", "Il6"), ("Tnf", "Il1b"), ("Il6", "Nlrp3"), ("Il6", "Bche"),
             ("Tnf", "Mbp"), ("Mbp", "Mog"), ("Il1b", "Nqo1"), ("Nqo1", "Pon1"),
             ("Tnf", "Nlrp3")]
    for a, b in edges:
        ax.plot([nodes[a][0], nodes[b][0]], [nodes[a][1], nodes[b][1]],
                color="#B0B5BD", lw=0.9, zorder=5)
    for n, (x, y, c) in nodes.items():
        ax.add_patch(Circle((x, y), 0.42 * s, fc=c, ec="white", lw=0.6, zorder=6))


def build(output_path):
    fig = plt.figure(figsize=(14, 9.2))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.set_aspect("equal")
    ax.axis("off")

    # ---------------------------------------------------------------- title
    ax.text(W / 2, 87.8, "Agentic LLM integration of multi-omics data resolves\na stress-to-demyelination network in Gulf War Illness",
            ha="center", va="center", fontsize=14.5, fontweight="bold",
            color=INK, linespacing=1.35)

    # ============================================================ column 1
    header(ax, 2, 30, 76, "1 · GWI model & datasets", C_BLUE)
    rbox(ax, 2.6, 15.5, 29.4, 74, fc="none", ec="#8A93A3", lw=0.9, r=1.2, z=1)

    rbox(ax, 3, 61, 29, 74, fc=C_LIGHT, ec="#8A93A3", lw=0.9)
    mouse_icon(ax, 8.5, 67.5, s=1.15)
    ax.text(19.5, 69.3, "BXD mice", ha="center", fontsize=8.5,
            fontweight="bold", color=INK)
    ax.text(19.5, 65.6, "CORT + DFP\n(GWI-emulants)", ha="center",
            fontsize=7.2, color=INK, linespacing=1.3)

    chips = [
        ("RNA-seq · 8 days", helix_icon, C_BLUE, "dataset 1"),
        ("320k methylation array · new", methylation_icon, C_TEAL, "dataset 2"),
        ("MBD-seq · 12 weeks", methylation_icon, C_TEAL, "dataset 3"),
        ("MBD-seq · 43 weeks", hourglass_icon, C_PURPLE, "dataset 4"),
    ]
    y_top = 56.5
    ch_h, gap = 6.6, 1.6
    for i, (label, icon, col, tag) in enumerate(chips):
        y0 = y_top - (i + 1) * ch_h - i * gap
        rbox(ax, 3, y0, 29, y0 + ch_h, fc="white", ec="#8A93A3", lw=0.9, r=1.0)
        icon(ax, 7.2, y0 + ch_h / 2, s=0.72, color=col)
        ax.text(16.2, y0 + ch_h * 0.60, label, ha="center", fontsize=7.0,
                color=INK)
        ax.text(16.2, y0 + ch_h * 0.22, tag, ha="center", fontsize=6.0,
                color=C_GREY, fontstyle="italic")

    arrow(ax, 1.6, 56.5, 1.6, 25.3, color=C_GREY_L, lw=1.0, ms=9)
    ax.text(0.5, 41, "exposure time-course", rotation=90, ha="center",
            va="center", fontsize=6.0, color=C_GREY)

    rbox(ax, 3, 17, 29, 23.5, fc="white", ec="#8A93A3", lw=0.9, r=1.0)
    ax.text(16, 21.4, "114 high-effect DMPs · 2,408 DE genes",
            ha="center", fontsize=6.4, color=INK)
    ax.text(16, 19.0, "≈ 3,000 methylation-QTL markers (LOD ≥ 4)",
            ha="center", fontsize=6.4, color=INK)

    # ============================================================ column 2
    header(ax, 36, 66, 76, "2 · Pinpoint — agentic LLM", C_ORANGE)
    rbox(ax, 36.6, 16, 65.4, 74, fc="#FFF7F2", ec=C_ORANGE, lw=1.1)
    ax.text(51, 71.5, "Pinpoint · Claude Opus 4.8", ha="center", fontsize=8.6,
            fontweight="bold", color=C_ORANGE)

    def agent(x0, x1, y0, y1, name, desc, color):
        rbox(ax, x0, y0, x1, y1, fc="white", ec=color, lw=1.1, r=1.0, z=5)
        ax.text((x0 + x1) / 2, (y0 + y1) / 2 + 1.15, name, ha="center",
                va="center", fontsize=7.4, fontweight="bold", color=color, zorder=6)
        ax.text((x0 + x1) / 2, (y0 + y1) / 2 - 1.15, desc, ha="center",
                va="center", fontsize=6.0, color="#4B5563", zorder=6)

    agent(43.5, 58.5, 64.5, 69.0, "Planner", "query → ordered steps", C_ORANGE)
    agent(38.5, 48.5, 50.5, 58.5, "Researcher", "GeneNetwork\nSPARQL", C_BLUE)
    agent(53.5, 63.5, 50.5, 58.5, "Expert", "NCBI\nEntrez", C_BLUE)
    agent(43.5, 58.5, 38.5, 44.5, "Reflector", "critiques · redirects", C_GREEN)
    agent(43.5, 58.5, 26.0, 32.0, "Manager", "routes · terminates", C_PURPLE)

    arrow(ax, 51, 64.5, 43.5, 58.5, color=C_GREY_L, lw=1.1, ms=9)
    arrow(ax, 51, 64.5, 58.5, 58.5, color=C_GREY_L, lw=1.1, ms=9)
    arrow(ax, 43.5, 50.5, 49, 44.5, color=C_GREY_L, lw=1.1, ms=9)
    arrow(ax, 58.5, 50.5, 53, 44.5, color=C_GREY_L, lw=1.1, ms=9)
    arrow(ax, 51, 38.5, 51, 32.0, color=C_GREY_L, lw=1.1, ms=9)
    # manager routing loop back to researcher / expert
    arrow(ax, 45.0, 29.0, 39.0, 29.0, color=C_GREY_L, lw=1.0, ms=8, ls=(0, (3, 2)))
    arrow(ax, 39.0, 29.0, 39.0, 54.5, color=C_GREY_L, lw=1.0, ms=8, ls=(0, (3, 2)))
    arrow(ax, 39.0, 54.5, 38.5, 54.5, color=C_GREY_L, lw=1.0, ms=8, ls=(0, (3, 2)))
    arrow(ax, 57.0, 29.0, 63.0, 29.0, color=C_GREY_L, lw=1.0, ms=8, ls=(0, (3, 2)))
    arrow(ax, 63.0, 29.0, 63.0, 54.5, color=C_GREY_L, lw=1.0, ms=8, ls=(0, (3, 2)))
    arrow(ax, 63.0, 54.5, 63.5, 54.5, color=C_GREY_L, lw=1.0, ms=8, ls=(0, (3, 2)))

    ax.text(51, 20.3, "3 iterations → consensus answer\nrepeated sampling ×10 → robust signals",
            ha="center", fontsize=6.6, color="#6B4326", linespacing=1.4)

    # arrow col1 -> col2
    arrow(ax, 30.0, 45, 36.4, 45, color="#5B6470", lw=1.8, ms=14)
    ax.text(33.2, 48.6, "gene lists ·\nQTL markers", ha="center", fontsize=6.4,
            color="#4B5563", linespacing=1.3)

    # ============================================================ column 3
    header(ax, 72, 100, 76, "3 · Genomic tasks", C_GREEN)

    # ---- task 1: fine-mapping
    rbox(ax, 73, 56, 99, 74, fc="white", ec="#8A93A3", lw=0.9)
    magnifier_icon(ax, 77.2, 67.8, s=0.95)
    ax.text(89, 71.5, "QTL fine-mapping", ha="center", fontsize=8.4,
            fontweight="bold", color=INK)
    ax.text(89, 67.5, "GEMMA mQTL (LOD ≥ 4)\n→ candidate genes in QTL regions",
            ha="center", fontsize=6.6, color="#4B5563", linespacing=1.4)
    rbox(ax, 75, 58, 97, 63.5, fc="#EAF3EC", ec=C_GREEN, lw=0.9, r=1.0)
    ax.text(86, 60.75, "×10 repeats → 78 prioritized genes", ha="center",
            fontsize=6.8, color="#2E5C26", fontweight="bold")

    # ---- task 2: network analysis
    rbox(ax, 73, 36.5, 99, 52.5, fc="white", ec="#8A93A3", lw=0.9)
    network_icon(ax, 77.8, 44.8, s=0.8)
    ax.text(88.5, 50.4, "Network analysis", ha="center", fontsize=8.4,
            fontweight="bold", color=INK)
    ax.text(90, 45.6, "genes across the 4 studies\n→ weighted gene–gene links\n(# shared annotations)",
            ha="center", fontsize=6.6, color="#4B5563", linespacing=1.4)
    rbox(ax, 75.5, 37, 98, 42, fc="#EAF3EC", ec=C_GREEN, lw=0.9, r=1.0)
    ax.text(86.75, 39.5, "×10 repeats → 53 genes · 47 edges · 5 modules",
            ha="center", fontsize=6.6, color="#2E5C26", fontweight="bold")

    # ---- task 3: PheWAS (classical QTL, linked to the network genes)
    rbox(ax, 73, 12, 99, 30.5, fc="white", ec=C_BLUE, lw=1.0, ls=(0, (4, 2.5)))
    db_icon(ax, 76.8, 27.4, s=0.62)
    ax.text(87.5, 28.2, "PheWAS · BXD phenome", ha="center", fontsize=8.0,
            fontweight="bold", color=INK)
    ax.text(87.5, 25.4, "GeneNetwork · 13,700 traits · precomputed GEMMA",
            ha="center", fontsize=5.9, color=C_BLUE)
    ax.text(87.5, 22.4, "traits with a QTL signal at\nnetwork-gene loci (± 1 Mbp)",
            ha="center", fontsize=6.4, color="#4B5563", linespacing=1.35)
    rbox(ax, 75, 13.6, 97, 19.2, fc="#EAF2FA", ec=C_BLUE, lw=0.9, r=1.0)
    ax.text(86, 16.4, "genes ↔ neurological & endocrine\nGWI-relevant traits",
            ha="center", fontsize=6.6, color="#274E6D", fontweight="bold",
            linespacing=1.35)

    # arrows col2 -> col3
    arrow(ax, 65.6, 65, 72.5, 65, color="#5B6470", lw=1.6, ms=12)
    arrow(ax, 65.6, 44.5, 72.5, 44.5, color="#5B6470", lw=1.6, ms=12)
    # arrow network analysis -> PheWAS
    arrow(ax, 86, 36.3, 86, 30.9, color=C_BLUE, lw=1.4, ms=10)
    ax.text(87.6, 33.6, "network genes", ha="left", fontsize=5.8, color=C_BLUE)

    # ============================================================ column 4
    header(ax, 106, 138, 76, "4 · Disease insights", C_PURPLE)

    rbox(ax, 107, 70.5, 137, 74.5, fc="#FFF3EA", ec=C_ORANGE, lw=1.0, r=1.2)
    ax.text(122, 72.5, "5 gene modules converge on IL6–TNF signalling",
            ha="center", va="center", fontsize=7.4, fontweight="bold",
            color="#9A3E00")

    def cascade_step(x0, x1, y0, y1, color, title, genes):
        rbox(ax, x0, y0, x1, y1, fc="white", ec=color, lw=1.3, r=1.2, z=4)
        ax.text((x0 + x1) / 2, (y0 + y1) / 2 + 1.55, title, ha="center",
                va="center", fontsize=7.3, fontweight="bold", color=color, zorder=5)
        ax.text((x0 + x1) / 2, (y0 + y1) / 2 - 1.65, genes, ha="center",
                va="center", fontsize=6.1, color="#4B5563", zorder=5)

    steps = [
        ("Exposure", "CORT + DFP (GWI-emulants)", C_GREY),
        ("Stress & cholinergic response", "Fkbp5 · Hspa1a/b · Bche · Chrm3", C_TEAL),
        ("Neuro-inflammation — IL6/TNF hub",
         "Tlr2 · Cd14 · Tnf · Il6 · Il1b · Nlrp3 · Casp4", C_ORANGE),
        ("Oxidative stress & detoxification",
         "Nqo1 · Hmox1 · Gsr · Txnrd1 · Pon1–3", C_GREEN),
        ("Demyelination", "Mbp · Mog · Plp1 ↓ → cognitive decline", C_PURPLE),
    ]
    y = 68.5
    bh, bgap = 9.2, 1.9
    for i, (title, genes, col) in enumerate(steps):
        y0 = y - (i + 1) * bh - i * bgap
        cascade_step(107, 137, y0, y0 + bh, col, title, genes)
        if i > 0:
            arrow(ax, 122, y0 + bh + bgap - 0.1, 122, y0 + bh + 0.1,
                  color=col, lw=1.6, ms=11)

    # arrow col3 -> col4 (network insight feeds the cascade)
    arrow(ax, 99.4, 47.2, 106.4, 47.2, color="#5B6470", lw=1.8, ms=14)

    save_figure(fig, output_path)
    plt.close(fig)
    print(f"saved {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-path",
        help="Path to save the graphical abstract at (PDF recommended)",
        default="graphical_abstract.pdf",
    )
    args = parser.parse_args()

    build(args.output_path)
