import argparse
import networkx as nx

from Bio.KEGG.KGML import KGML_parser
import pickle



parser = argparse.ArgumentParser(description='Create a DiGraph from KGML file.')
parser.add_argument('--keggpickles', nargs='+', type=argparse.FileType('r'), required=True)
parser.add_argument('--outpickle', type=argparse.FileType('w'), required=True)

args = parser.parse_args()


entries   = {}
relations = {}
for p in args.keggpickles:
    g = nx.gpickle.read_gpickle(p)
    for n in g.nodes():
        entries[n.graphics[0].name] = n

    for e in g.edges():
        relations[(e[0].graphics[0].name, e[1].graphics[0].name)] = g.get_edge_data(*e)


h = nx.DiGraph()
for e in entries:
    h.add_node(entries[e])

for r in relations:
    h.add_edge( entries[r[0]], entries[r[1]], relation=relations[r])

# write to pickle
nx.gpickle.write_gpickle(h, args.outpickle)
