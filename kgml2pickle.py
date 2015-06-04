import argparse
import networkx as nx
from keggparser import parse_KGML
import pickle


parser = argparse.ArgumentParser(description='Create a DiGraph from KGML file.')
parser.add_argument('--kgml', type=argparse.FileType('r'), required=True)
parser.add_argument('--pickle', type=argparse.FileType('w'), required=True)

args = parser.parse_args()

g = parse_KGML.KGML2Graph(args.kgml)[1]

print g.nodes()

nx.write_gpickle(g,
                 "/tmp/hola.p")

#pickle.dump(g, args.pickle)

