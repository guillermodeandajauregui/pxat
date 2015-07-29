import networkx as nx


def pathway_jin( p1, p2, g ):

    pw1  = []
    pw2  = []

    for n in g.nodes():
        if p1 in g.node[n]['pathways']:
            pw1.append(n)
        if p2 in g.node[n]['pathways']:
            pw2.append(n)

    return float(len(set(pw1).intersection(set(pw2)))) / float(len(set(pw1).union(set(pw2))))




def pathway_xtalk( p, g):
    """
    returns pathways that crosstalk with given pathway, with their jin
    """



# class KEGGPathway(DiGraph):
#     def __init__(g):
#         self.g = g


#     def paths(self, s, t):
#         pass


# k = KEGGPathway(g)

# k.paths(s, t)
    
def nodes_from_pathway(p, g):
    nodes = []
    for n in g.nodes():
        if p in g.node[n]['pathways']:
            nodes.append(n)
    return nodes


def subgraph_from_pathways(pathways, g):
    nodes = []
    if type(pathways) == list:
        for p in pathways:
            nodes += nodes_from_pathway(p, g)
    elif type(pathways) == str:
        nodes += nodes_from_pathway(pathways, g)
        
    return nx.subgraph(g, nodes)


def paths(source, target, g):
    return [p for p in nx.all_simple_paths(g, source, target)]



def shortest_path(s, t, g):
    pass
    
def effector_from_signal(signal, g):
    signals   = set()
    effectors = set()

    for n in g.nodes():
        if 'signal' in g.node[n]:
            signals.add(n)

        if 'effector' in g.node[n]:
            effectors.add(n)

        for s in signals:
            for e in effectors:
                print s, e, [p for p in nx.all_simple_paths(g, s, e)]



def signals_from_effector(effector, g):
    pass



def effectors_from_receptor(receptor, g):
    pass



def receptors_from_effectors(effector, g):
    pass




def signals_from_pathway(p, g):
    pass



def effectors_from_pathway(p, g):
    pass


def receptors_from_pathway(p, g):
    pass


def transducers_from_pathway(p, g):
    pass


def final_transducers_from_pathway(p, g):
    pass
