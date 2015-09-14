import argparse
import networkx as nx
from pxat import merge_pathway_subgraphs


parser = argparse.ArgumentParser(description='Merge pickles.')
parser.add_argument('--keggpickles', nargs='+', type=argparse.FileType('r'), required=True)
parser.add_argument('--outpickle', type=argparse.FileType('w'), required=True)

args = parser.parse_args()

h = merge_pathway_subgraphs(*[nx.gpickle.read_gpickle(p) for p in args.keggpickles])

nx.gpickle.write_gpickle(h, args.outpickle)
