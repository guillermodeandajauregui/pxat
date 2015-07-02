import argparse
import networkx as nx

from Bio.KEGG.KGML import KGML_parser
import pickle



parser = argparse.ArgumentParser(description='Create a DiGraph from KGML file.')
parser.add_argument('--kgml', type=argparse.FileType('r'), required=True)
parser.add_argument('--pickle', type=argparse.FileType('w'), required=True)

args = parser.parse_args()


#then parse KGML file and return pathway object
pathway = KGML_parser.parse(args.kgml).next()


# make genes unique
uentries = {}
for e in pathway.entries:
    k = pathway.entries[e].graphics[0].name
    uentries[k] = e


g = nx.DiGraph()
# copy relations to edges
for relation in pathway.relations:
    g.add_edge( uentries[relation.entry1.graphics[0].name],
                uentries[relation.entry2.graphics[0].name], relation=relation)

# write to pickle
nx.gpickle.write_gpickle(g, args.pickle)
