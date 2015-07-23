import argparse
import networkx as nx

parser = argparse.ArgumentParser(description='Topological classification of nodes in pathways.')
parser.add_argument('--input', type=argparse.FileType('r'), required=True)
parser.add_argument('--out', type=argparse.FileType('w'), required=True)
parser.add_argument('--format', choices=['edgelist', 'gml', 'tsv'], default='edgelist')

args = parser.parse_args()

g = nx.gpickle.read_gpickle(args.input)


if args.format == 'edgelist':
    nx.readwrite.edgelist.write_edgelist(g, args.out, data=True)

elif args.format == 'gml':
    nx.readwrite.gml.write_gml(g, args.out)

elif args.format == 'tsv':
    for e in g.edges():
        edge_data = g.get_edge_data(*e)
        linea = "%s\t%s\t%s\n" % (e[0],e[1], edge_data)
        args.out.write(linea)
