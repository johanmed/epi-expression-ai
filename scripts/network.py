"""Network visualisation for the GWI gene-connection analysis.

Draws a STRING-style circular network figure: nodes are placed on a circle
(grouped by connected component and ordered to minimise chord crossings),
node colours encode curated functional groups, node size encodes total
connection weight, and edge width/shade encode the connection weight.
Outputs a raster figure at the requested path plus a vector PDF twin
(recommended for LaTeX inclusion: small file, crisp at any scale, fast
compilation).
"""

import argparse
import json
import math
import os

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import networkx as nx
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import PathPatch
from matplotlib.path import Path

# Publishable TrueType fonts embedded in the PDF (editable text, small files)
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

DPI = 300  # raster output resolution (use the PDF twin for LaTeX)

GROUPS = {
    "Fkbp5": "Stress & inflammatory signalling", "Hspa1a": "Stress & inflammatory signalling",
    "Hspa1b": "Stress & inflammatory signalling", "Tnf": "Stress & inflammatory signalling",
    "Il1b": "Stress & inflammatory signalling", "Il6": "Stress & inflammatory signalling",
    "Nlrp3": "Stress & inflammatory signalling", "Nqo1": "Stress & inflammatory signalling",
    "Hmox1": "Stress & inflammatory signalling", "Gsr": "Stress & inflammatory signalling",
    "Txnrd1": "Stress & inflammatory signalling", "Srxn1": "Stress & inflammatory signalling",
    "Hsd11b1": "Stress & inflammatory signalling", "Casp4": "Stress & inflammatory signalling",
    "Casp12": "Stress & inflammatory signalling", "Cav1": "Stress & inflammatory signalling",
    "Cav2": "Stress & inflammatory signalling",
    "Bche": "Cholinergic signalling", "Chrm3": "Cholinergic signalling",
    "Chrna10": "Cholinergic signalling", "Chrnb2": "Cholinergic signalling",
    "Colq": "Cholinergic signalling",
    "Trem2": "Microglial activation", "Aif1": "Microglial activation",
    "Cd68": "Microglial activation", "P2ry12": "Microglial activation",
    "Tmem119": "Microglial activation",
    "Apoa1": "Lipoprotein metabolism", "Apoc1": "Lipoprotein metabolism",
    "Apoc3": "Lipoprotein metabolism", "Apoe": "Lipoprotein metabolism",
    "Ldlr": "Lipoprotein metabolism",
    "Hmgcr": "Cholesterol biosynthesis", "Sqle": "Cholesterol biosynthesis",
    "Fdft1": "Cholesterol biosynthesis", "Insig1": "Cholesterol biosynthesis",
    "Aqp4": "Astrocyte glutamate transport", "Gfap": "Astrocyte glutamate transport",
    "Slc1a2": "Astrocyte glutamate transport", "Slc1a3": "Astrocyte glutamate transport",
    "Pon1": "Paraoxonase detoxification", "Pon2": "Paraoxonase detoxification",
    "Pon3": "Paraoxonase detoxification",
    "Mbp": "Myelination", "Mog": "Myelination", "Plp1": "Myelination",
    "Ccl2": "Chemokine signalling", "Ccl3": "Chemokine signalling",
    "Cxcl10": "Chemokine signalling",
    "Nfatc4": "Calcineurin\u2013NFAT signalling", "Rcan1": "Calcineurin\u2013NFAT signalling",
    "Cd14": "Toll-like receptor pathway", "Tlr2": "Toll-like receptor pathway",
}

GROUP_ORDER = [
    "Stress & inflammatory signalling", "Cholinergic signalling",
    "Microglial activation", "Lipoprotein metabolism", "Cholesterol biosynthesis",
    "Astrocyte glutamate transport", "Paraoxonase detoxification", "Myelination",
    "Chemokine signalling", "Calcineurin\u2013NFAT signalling", "Toll-like receptor pathway",
    "Other",
]

PALETTE = {
    "Stress & inflammatory signalling": "#D55C00",
    "Cholinergic signalling": "#008B8B",
    "Microglial activation": "#8B6914",
    "Lipoprotein metabolism": "#3B6BA5",
    "Cholesterol biosynthesis": "#7B4FA3",
    "Astrocyte glutamate transport": "#4A9C3D",
    "Paraoxonase detoxification": "#C2316B",
    "Myelination": "#666666",
    "Chemokine signalling": "#B07A00",
    "Calcineurin\u2013NFAT signalling": "#2E9CB8",
    "Toll-like receptor pathway": "#946743",
    "Other": "#999999",
}

NODE_EDGE = "#303030"
NODE_SIZE_BASE = 60    # scatter size for a node with zero strength
NODE_SIZE_SCALE = 11   # scatter size added per unit of total connection weight
# Representative "total connection weight" values shown in the node-size legend
SIZE_LEGEND_VALUES = (5, 15, 30, 45)

# edge shade/width interpolation range (data-driven, see draw_network)
EDGE_SHADE_MIN, EDGE_SHADE_MAX = 0.38, 0.83
EDGE_WIDTH_MIN, EDGE_WIDTH_MAX = 0.55, 2.55

LABEL_RADIUS = 1.06
LABEL_HALO = [pe.withStroke(linewidth=2.5, foreground="white")]


def node_group(node):
    return GROUPS.get(node, "Other")


def node_color(node):
    return PALETTE.get(node_group(node), PALETTE["Other"])


def save_figure(fig, path, dpi=DPI, tight=True):
    """Save to `path`; if not already a PDF, also write a vector PDF twin."""
    fig.savefig(path, dpi=dpi, bbox_inches="tight" if tight else "",
                facecolor="white")
    root, ext = os.path.splitext(path)
    if ext.lower() != ".pdf":
        fig.savefig(root + ".pdf", bbox_inches="tight" if tight else "",
                    facecolor="white")


def parse_edge_list(input_path):
    """Parse the network-analysis JSON into an edge DataFrame (weight >= 5)."""
    with open(input_path) as i:
        dic = json.loads(i.read())

    collection = []
    for key in dic:
        gene1, gene2 = key.split("->")
        value = dic[key]
        weight = int(value.split("=")[1].split(",")[0].strip())
        if weight >= 5:
            collection.append([gene1.strip(), gene2.strip(), weight])
        else:
            continue
    df = pd.DataFrame(collection, columns=["source", "target", "weight"])
    df["norm_weight"] = df["weight"] / df["weight"].max()
    return df


def build_graph(df):
    """Undirected display graph; reciprocal/duplicate edges keep max weight."""
    G = nx.Graph()
    for _, row in df.iterrows():
        u, v, w = row["source"], row["target"], row["weight"]
        if G.has_edge(u, v):
            G[u][v]["weight"] = max(G[u][v]["weight"], w)
        else:
            G.add_edge(u, v, weight=w)
    return G


def circular_positions(G):
    """Nodes on a unit circle: components contiguous (largest first), nodes
    inside each component ordered by a spring-embedding angle to minimise
    chord crossings."""
    components = sorted(nx.connected_components(G),
                         key=lambda c: (-len(c), sorted(c)[0]))
    nodes = []
    for comp in components:
        if len(comp) <= 2:
            nodes.extend(sorted(comp))
            continue
        loc = nx.spring_layout(G.subgraph(comp), seed=7)
        nodes.extend(sorted(comp, key=lambda n: math.atan2(loc[n][1], loc[n][0])))

    n = len(nodes)
    pos = {}
    for i, node in enumerate(nodes):
        theta = 2 * math.pi * i / n
        pos[node] = (math.cos(theta), math.sin(theta))
    return pos


def draw_network(G, pos):
    """Render the STRING-style circular network with three legends."""
    strength = {n: sum(d["weight"] for _, _, d in G.edges(n, data=True))
                for n in G.nodes()}
    weights = [d["weight"] for _, _, d in G.edges(data=True)]
    wmin, wmax = min(weights), max(weights)

    def edge_width(w):
        return EDGE_WIDTH_MIN + (EDGE_WIDTH_MAX - EDGE_WIDTH_MIN) * (w - wmin) / (wmax - wmin)

    def edge_shade(w):
        return EDGE_SHADE_MIN + (EDGE_SHADE_MAX - EDGE_SHADE_MIN) * (w - wmin) / (wmax - wmin)

    fig, ax = plt.subplots(figsize=(10, 10))

    # curved edges, grey shade and width scaled by connection weight
    for u, v, d in G.edges(data=True):
        w = d["weight"]
        col = plt.cm.Greys(edge_shade(w))
        mx, my = (pos[u][0] + pos[v][0]) / 2, (pos[u][1] + pos[v][1]) / 2
        cx, cy = 0.22 * mx, 0.22 * my
        bx = 0.25 * pos[u][0] + 0.5 * cx + 0.25 * pos[v][0]
        by = 0.25 * pos[u][1] + 0.5 * cy + 0.25 * pos[v][1]
        patch = PathPatch(
            Path([(pos[u][0], pos[u][1]), (bx, by), (pos[v][0], pos[v][1])],
                 [Path.MOVETO, Path.CURVE3, Path.CURVE3]),
            fc="none", ec=col, lw=edge_width(w), alpha=0.85, zorder=1,
            capstyle="round")
        ax.add_patch(patch)

    # nodes: colour = functional group, size = total connection weight
    for node, (x, y) in pos.items():
        s = NODE_SIZE_BASE + NODE_SIZE_SCALE * strength[node]
        ax.scatter(x, y, s=s, c=node_color(node), edgecolors=NODE_EDGE,
                   linewidths=0.8, zorder=3)

    # radial gene labels (never upside-down), white halo for legibility
    for node, (x, y) in pos.items():
        deg = math.degrees(math.atan2(y, x))
        lx, ly = LABEL_RADIUS * x, LABEL_RADIUS * y
        if -90 <= deg <= 90:
            rot, ha = deg, "left"
        else:
            rot, ha = deg + 180, "right"
        ax.text(lx, ly, node, rotation=rot, rotation_mode="anchor", ha=ha,
                va="center", fontsize=8, fontstyle="italic", color="#1a1a1a",
                zorder=4, path_effects=LABEL_HALO)

    # legend 1: functional groups (node colours)
    present = {node_group(n) for n in G.nodes()}
    ordered = [g for g in GROUP_ORDER if g in present] + \
             sorted(present - set(GROUP_ORDER))
    group_handles = [
        Line2D([], [], marker="o", ls="", markersize=8,
               markerfacecolor=PALETTE[g], markeredgecolor=NODE_EDGE, label=g)
        for g in ordered
    ]
    legend_groups = ax.legend(
        handles=group_handles, loc="center left", bbox_to_anchor=(1.04, 0.5),
        frameon=False, fontsize=9, title="Functional group", title_fontsize=10)
    ax.add_artist(legend_groups)

    # legend 2: node size (total connection weight)
    size_handles = [
        Line2D([], [], marker="o", ls="",
               markersize=math.sqrt(NODE_SIZE_BASE + NODE_SIZE_SCALE * v),
               markerfacecolor="#B0B0B0", markeredgecolor=NODE_EDGE, label=str(v))
        for v in SIZE_LEGEND_VALUES
    ]
    legend_sizes = ax.legend(
        handles=size_handles, loc="upper right", bbox_to_anchor=(0.44, -0.01),
        frameon=False, fontsize=8, ncol=4, labelspacing=1.4, handletextpad=0.9,
        columnspacing=1.4, title="Node size: total connection weight",
        title_fontsize=9)
    ax.add_artist(legend_sizes)

    # legend 3: edge width (connection weight)
    weight_ticks = sorted({round(wmin + i * (wmax - wmin) / 3) for i in range(4)})
    weight_handles = [
        Line2D([], [], color=plt.cm.Greys(edge_shade(w)), lw=edge_width(w),
               label=str(w))
        for w in weight_ticks
    ]
    ax.legend(handles=weight_handles, loc="upper left", bbox_to_anchor=(0.56, -0.01),
              frameon=False, fontsize=8, ncol=4, title="Edge width: connection weight",
              title_fontsize=9)

    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.36, 1.36)
    ax.set_aspect("equal")
    ax.axis("off")
    return fig


def plot_connectivity_metrics(df, G_metrics):
    """Network connectivity metrics figure (bar plot per gene)."""
    degree_centrality = sorted(
        nx.degree_centrality(G_metrics).items(), key=lambda item: item[0]
    )
    betweenness_centrality = sorted(
        nx.betweenness_centrality(G_metrics).items(), key=lambda item: item[0]
    )
    eigenvector_centrality = sorted(
        nx.eigenvector_centrality(G_metrics, max_iter=1000).items(),
        key=lambda item: item[0]
    )
    degree_clustering = sorted(nx.clustering(G_metrics).items(),
                               key=lambda item: item[0])

    final_metrics = pd.DataFrame(
        {
            "Genes": list(dict(degree_centrality).keys()),
            "degree_centrality": list(dict(degree_centrality).values()),
            "betweenness_centrality": list(dict(betweenness_centrality).values()),
            "eigenvector_centrality": list(dict(eigenvector_centrality).values()),
            "degree_clustering": list(dict(degree_clustering).values()),
        }
    )
    fig = final_metrics.plot.bar(
        x="Genes", figsize=(10, 10), title="Network connectivity metrics"
    ).figure
    plt.ylabel("Scores")
    plt.tight_layout()
    return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-path", help="Path to file with network analysis metrics results"
    )
    parser.add_argument("--network-path", help="Path to save network figure at")
    parser.add_argument("--metric-path", help="Path to save network connectivity figure at")
    args = parser.parse_args()

    df = parse_edge_list(args.input_path)

    # ------------------------------------------------------------------ network
    G = build_graph(df)
    pos = circular_positions(G)
    fig = draw_network(G, pos)
    save_figure(fig, args.network_path)
    plt.show()

    # ------------------------------------------------- connectivity metrics
    # Metrics are computed on the directed edge list, as in the original
    # analysis (reciprocal connections kept as directed edges).
    G_metrics = nx.from_pandas_edgelist(
        df,
        source="source",
        target="target",
        edge_attr="norm_weight",
        create_using=nx.DiGraph,
    )
    fig_metrics = plot_connectivity_metrics(df, G_metrics)
    fig_metrics.savefig(args.metric_path, dpi=DPI, facecolor="white")
    root, ext = os.path.splitext(args.metric_path)
    if ext.lower() != ".pdf":
        fig_metrics.savefig(root + ".pdf", facecolor="white")
    plt.show()
