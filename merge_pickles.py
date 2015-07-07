import argparse
import networkx as nx
from Bio.KEGG.KGML import KGML_parser


parser = argparse.ArgumentParser(description='Create a DiGraph from KGML file.')
parser.add_argument('--keggpickles', nargs='+', type=argparse.FileType('r'), required=True)
parser.add_argument('--outpickle', type=argparse.FileType('w'), required=True)

args = parser.parse_args()


# new graph
h = nx.DiGraph()

for p in args.keggpickles:
    g = nx.gpickle.read_gpickle(p)
    for e in g.edges():
        h.add_edge( *e, attr_dict=g.get_edge_data(*e))

# write to pickle
nx.gpickle.write_gpickle(h, args.outpickle)
