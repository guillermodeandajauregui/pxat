import argparse
import networkx as nx
from Bio.KEGG.KGML import KGML_parser



parser = argparse.ArgumentParser(description='Create a DiGraph from KGML file.')
parser.add_argument('--kgml', type=argparse.FileType('r'), required=True)
parser.add_argument('--pickle', type=argparse.FileType('w'), required=True)

args = parser.parse_args()

# a shorter way
def readKGML(kgml):
    return KGML_parser.parse(kgml).next()


# then parse KGML file and return pathway object
def kgml2graph(pathway):
    g = nx.DiGraph()
    # copy relations to edges
    for relation in pathway.relations:
        # no undefined or path entries, only hsa names
        if relation.entry1.name.startswith('hsa') and relation.entry2.name.startswith('hsa'):
            for e1 in relation.entry1.name.split():
                for e2 in relation.entry2.name.split():
                    g.add_node( e1, pathways=set([pathway.title,]))
                    g.add_node( e2, pathways=set([pathway.title,]))
                    g.add_edge( e1,
                                e2, type=relation.type,
                                subtypes=relation.subtypes,
                                pathways=set([pathway.title,]))

    return g





# write to pickle
nx.gpickle.write_gpickle( kgml2graph( readKGML(args.kgml) ),
                          args.pickle)





