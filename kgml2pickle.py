import argparse
import networkx as nx
from Bio.KEGG.KGML import KGML_parser



parser = argparse.ArgumentParser(description='Create a DiGraph from KGML file.')
parser.add_argument('--kgml', type=argparse.FileType('r'), required=True)
parser.add_argument('--pickle', type=argparse.FileType('w'), required=True)

args = parser.parse_args()


#then parse KGML file and return pathway object
pathway = KGML_parser.parse(args.kgml).next()

g = nx.DiGraph()
# copy relations to edges
for relation in pathway.relations:
    # no undefined or path entries, only hsa names
    if relation.entry1.name.startswith('hsa') and relation.entry2.name.startswith('hsa'):
        for e1 in relation.entry1.name.split():
            for e2 in relation.entry2.name.split():
                g.add_edge( e1,
                            e2, type=relation.type, subtypes=relation.subtypes)

# write to pickle
nx.gpickle.write_gpickle(g, args.pickle)
