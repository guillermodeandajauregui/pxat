import argparse
import networkx as nx

parser = argparse.ArgumentParser(description='Topological classification of nodes in pathways.')
parser.add_argument('--input', type=argparse.FileType('r'), required=True)
parser.add_argument('--out', type=argparse.FileType('w'), required=True)
parser.add_argument('--format', choices=['edgelist', 'gml'], default='edgelist')

args = parser.parse_args()

g = nx.gpickle.read_gpickle(args.input)


if args.format == 'edgelist':
    nx.readwrite.edgelist.write_edgelist(g, args.out, data=True)

elif args.format == 'gml':
    nx.readwrite.gml.write_gml(g, args.out)
