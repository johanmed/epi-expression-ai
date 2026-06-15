import argparse
import json

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-path", help="Path to file with network analysis metrics results"
    )
    parser.add_argument("--network-path", help="Path to save network figure at")
    parser.add_argument("--metric-path", help="Path to save network figure at")
    args = parser.parse_args()

    with open(args.input_path) as i:
        content = i.read()
    dic = json.loads(content)

    collection = []
    for key in dic:
        gene1, gene2 = key.split("->")
        value = dic[key]
        weight = int(value.split("=")[1].split(",")[0].strip())
        if weight >= 3:
            collection.append([gene1.strip(), gene2.strip(), weight])
        else:
            continue
    df = pd.DataFrame(collection, columns=["source", "target", "weight"])

    G = nx.from_pandas_edgelist(
        df,
        source="source",
        target="target",
        edge_attr="weight",
        create_using=nx.DiGraph(),
    )
    
    pos = nx.spiral_layout(G, resolution=0.6, scale=20, equidistant=True)
    plt.figure(figsize=(10, 10))
    nx.draw_networkx_nodes(G, pos, node_color="lightblue", node_size=1000)
    nx.draw_networkx_labels(G, pos, font_size=10, font_weight="bold")
    nx.draw_networkx_edges(
        G,
        pos,
        arrowstyle="->",
        arrowsize=20,
        edge_color="gray",
        connectionstyle="angle,rad=0.2",
    )

    plt.title("Inferred Biological Connections in Gulf War Disease", fontsize=20)
    plt.axis("off")
    plt.savefig(args.network_path, dpi=1000)
    plt.show()

    # Compute network connectivity metrics
    degree_centrality = sorted(nx.degree_centrality(G).items(), key=lambda item: item[0])
    betweenness_centrality = sorted(nx.betweenness_centrality(G).items(), key=lambda item: item[0])
    eigenvector_centrality = sorted(nx.eigenvector_centrality(G, max_iter=1000).items(), key=lambda item: item[0])
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
        x="Genes",
        figsize=(10, 10),
        title="Network connectivity metrics"
    )
    plt.savefig(args.metric_path, dpi=1000)
    plt.show()


        
