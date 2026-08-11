import argparse
import json

import matplotlib.pyplot as plt
import pandas as pd

import networkx as nx

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-path", help="Path to file with network analysis metrics results"
    )
    parser.add_argument("--network-path", help="Path to save network figure at")
    parser.add_argument("--metric-path", help="Path to save network connectivity figure at")
    args = parser.parse_args()

    with open(args.input_path) as i:
        content = i.read()
    dic = json.loads(content)

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
    df["norm_weight"] = df["weight"]/df["weight"].max()

    # Draw network
    G = nx.from_pandas_edgelist(
        df,
        source="source",
        target="target",
        edge_attr="norm_weight",
        create_using=nx.DiGraph(),
    )

    node_weights = {node: 0 for node in G.nodes()}
    for u, v, data in G.edges(data=True):
        w = data["norm_weight"]
        node_weights[u] += w
        node_weights[v] += w
    node_color_values = [node_weights[n] for n in G.nodes()]

    edge_weights = nx.get_edge_attributes(G, "norm_weight")
    min_w = min(edge_weights.values())
    max_w = max(edge_weights.values())
    width_map = {
        e: 0.5 + 3.0 * (w - min_w) / (max_w - min_w + 1e-9)
        for e, w in edge_weights.items()
    }

    pos = nx.spiral_layout(
        G.to_undirected(),
        scale=50,
        resolution=1,
        equidistant=True,
    )

    plt.figure(figsize=(20, 15))

    nx.draw_networkx_edges(
        G,
        pos,
        arrowstyle="-",
        arrowsize=30,
        edge_color="grey",
        alpha=0.95,
        connectionstyle="angle,rad=0.2",
        node_size=2000,
        width=[width_map[e] for e in G.edges()],
    )

    nodes = nx.draw_networkx_nodes(
        G,
        pos,
        node_color=node_color_values,
        cmap=plt.cm.YlOrRd,
        node_size=2000,
        alpha=0.95,
        linewidths=2,
    )

    plt.colorbar(nodes, shrink=0.5, label="Total connected weight")

    nx.draw_networkx_labels(
        G,
        pos,
        font_size=20,
        font_weight="bold",
    )

    plt.suptitle("Inferred Biological Connections in Gulf War Disease", fontsize=40)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(args.network_path, dpi=1000)
    plt.show()

    # Compute network connectivity metrics
    degree_centrality = sorted(
        nx.degree_centrality(G).items(), key=lambda item: item[0]
    )
    betweenness_centrality = sorted(
        nx.betweenness_centrality(G).items(), key=lambda item: item[0]
    )
    eigenvector_centrality = sorted(
        nx.eigenvector_centrality(G, max_iter=1000).items(), key=lambda item: item[0]
    )
    degree_clustering = sorted(nx.clustering(G).items(), key=lambda item: item[0])

    final_metrics = pd.DataFrame(
        {
            "Genes": list(dict(degree_centrality).keys()),
            "degree_centrality": list(dict(degree_centrality).values()),
            "betweenness_centrality": list(dict(betweenness_centrality).values()),
            "eigenvector_centrality": list(dict(eigenvector_centrality).values()),
            "degree_clustering": list(dict(degree_clustering).values()),
        }
    )
    final_metrics.plot.bar(
        x="Genes", figsize=(10, 10), title="Network connectivity metrics"
    )
    plt.ylabel("Scores")
    plt.tight_layout()
    plt.savefig(args.metric_path, dpi=1000)
    plt.show()
