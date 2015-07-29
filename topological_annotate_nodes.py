import argparse
import networkx as nx
from pxat import topological_annotate

parser = argparse.ArgumentParser(description='Topological classification of nodes in pathways.')
parser.add_argument('--inpickle', type=argparse.FileType('r'), required=True)
parser.add_argument('--outpickle', type=argparse.FileType('w'), required=True)

args = parser.parse_args()

g = nx.gpickle.read_gpickle(args.inpickle)

# write to pickle
nx.gpickle.write_gpickle(topological_annotate(g), args.outpickle)

