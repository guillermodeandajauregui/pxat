import argparse
import networkx as nx
from pxat import kgml2graph, readKGML


parser = argparse.ArgumentParser(description='Create a DiGraph from KGML file.')
parser.add_argument('--kgml', type=argparse.FileType('r'), required=True)
parser.add_argument('--pickle', type=argparse.FileType('w'), required=True)

args = parser.parse_args()



# write to pickle
nx.gpickle.write_gpickle( kgml2graph( readKGML(args.kgml) ),
                          args.pickle)





