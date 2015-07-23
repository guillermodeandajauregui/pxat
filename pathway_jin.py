import networkx as nx





def pathway_jin( p1, p2 ):
    g = nx.gpickle.read_gpickle('human_kegg_ann.pickle')
    pw1 = []
    pw2  = []

    for n in g.nodes():
        if p1 in g.node[n]['pathways']:
            pw1.append(n)
        if p2 in g.node[n]['pathways']:
            pw2.append(n)

    return float(len(set(pw1).intersection(set(pw2)))) / float(len(set(pw1).union(set(pw2))))


p1 = 'Apoptosis'
p2 = 'Estrogen signaling pathway'

print pathway_jin(p1, p2)
