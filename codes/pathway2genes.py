"""
Utility script to extract genes related to pathways identified from methylation analysis
"""

import requests


pathways = input("Please enter pathway identifier separated by comma:\n")
pathways = pathways.split(", ")

def fetch_process(pathway_id: str) -> list:
    response = requests.get(f"https://rest.kegg.jp/link/mmu/{pathway_id}")
    if response.status_code == 200:
        fetched = response.text.strip()
        genes = []
        results = fetched.split("\n")
        for result in results:
            gene = (result.split("\t")[1].split(":")[1]).strip()
            genes.append(gene)
        return genes
    else:
        raise ValueError("The pathway id is not valid")

new_list = []
for pathway in pathways:
    genes = fetch_process(pathway)
    new_list.extend(genes)
print(new_list)
