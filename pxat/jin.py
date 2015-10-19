def pathway_jin( p1, p2, g ):
    """
    returns Jaccard index between two canonical pathways in the network
    """
    pw1  = []
    pw2  = []

    for n in g.nodes():
        if p1 in g.node[n]['pathways']:
            pw1.append(n)
        if p2 in g.node[n]['pathways']:
            pw2.append(n)

    return float(len(set(pw1).intersection(set(pw2)))) / float(len(set(pw1).union(set(pw2))))





def graph_jin(g1, g2):
    """
    returns Jaccard index between two graphs
    
    """
    return float(len(set(g1.node).intersection(set(g2.node)))) / float(len(set(g1.node).union(set(g2.node))))







def subgraph_jin(pw1, pw2, g):
    """
    returns Jaccard index between two subgraphs from two nbunches 
    
    """
    p1 = disubgraph(g, pw1).node
    p2 = disubgraph(g, pw2).node
    return float(len(set(p1).intersection(set(p2)))) / float(len(set(p1).union(set(p2))))
