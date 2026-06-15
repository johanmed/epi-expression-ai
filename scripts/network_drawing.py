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
    parser.add_argument("--output-path", help="Path to save network figure at")
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
    plt.savefig(args.output_path, dpi=1000)
    plt.show()
