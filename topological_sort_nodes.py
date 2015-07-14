import argparse
import networkx as nx

parser = argparse.ArgumentParser(description='Topological classification of nodes in pathways.')
parser.add_argument('--inpickle', type=argparse.FileType('r'), required=True)
parser.add_argument('--outpickle', type=argparse.FileType('w'), required=True)

args = parser.parse_args()

g = nx.gpickle.read_gpickle(args.inpickle)
h = g.copy()

for n in g.nodes():
    if len(g.in_edges(n)) == 0:
        h.node[n]['signal'] = True
        for e in g.out_edges(n):
            h.node[e[1]]['receptor'] = True

for n in g.nodes():
    if len(g.out_edges(n)) == 0:
        h.node[n]['effector'] = True
        for e in g.in_edges(n):
            h.node[e[0]]['final_transducer'] = True        

for n in h.nodes():
    if not ('signal' in h.node[n] or 'receptor' in h.node[n] or 'effector' in h.node[n]):
        h.node[n]['transducer'] = True


# write to pickle
nx.gpickle.write_gpickle(h, args.outpickle)

