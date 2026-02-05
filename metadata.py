import json
from Bio import Phylo

# make metadata.json
file = f'nifH_final.newick'
tree = Phylo.read(file, "newick")

metadata = {}

for leaf in tree.get_terminals():

    id = leaf.name.replace("_", " ")
    genome = f"'{id.split('|')[0]}'"
    metadata[id] = {'species': 'test', 'genome': 'test'}

with open('metadata.json', 'w') as f:
    out = json.dump(metadata, f, indent=2)
