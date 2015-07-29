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


def topological_annotate(g):
    # clear previous annotation if there be
    for n in g.nodes():
        if 'transducer' in g.node[n]:
            del(g.node[n]['transducer'])
        if 'signal' in g.node[n]:
            del(g.node[n]['signal'])
        if 'receptor' in g.node[n]:
            del(g.node[n]['receptor'])
        if 'effector' in g.node[n]:
            del(g.node[n]['effector'])
        if 'final_transducer' in g.node[n]:
            del(g.node[n]['final_transducer'])
                
    
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
        if not ('signal' in h.node[n] or 'effector' in h.node[n]):
            h.node[n]['transducer'] = True
    
    # Biological relevance of molecules that are both receptor AND
    # transducers is different from "true" receptors 
    for n in h.nodes():
        if 'receptor' in h.node[n] and 'transducer' in h.node[n]:
            in_classes = set()
            for e in h.in_edges(n):
                in_classes.update(h.node[e[0]].keys())
    
            if not 'transducer' in in_classes:
                del(h.node[n]['transducer'])
    return h




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
        
    h = nx.subgraph(g, nodes)

    remove_edges = []
    if type(pathways) == list:
        for p in pathways:
            for e in h.edges():
                if not p in h.get_edge_data(*e)['pathways']:
                    remove_edges.append(e)
    elif type(pathways) == str:
        for e in h.edges():
            if not pathways in h.get_edge_data(*e)['pathways']:
                remove_edges.append(e)
                
    for r in remove_edges:
        h.remove_edge(*r)

    return topological_annotate(h)


def paths(source, target, g):
    return [p for p in nx.all_simple_paths(g, source, target)]


    
# def effector_from_signal(signal, g):
#     signals   = set()
#     effectors = set()

#     for n in g.nodes():
#         if 'signal' in g.node[n]:
#             signals.add(n)

#         if 'effector' in g.node[n]:
#             effectors.add(n)

#         for s in signals:
#             for e in effectors:
#                 print s, e, [p for p in nx.all_simple_paths(g, s, e)]



def signals_to_node(node, g):
    all_signals = set()
    for n in g.nodes():
        if 'signal' in g.node[n]:
            all_signals.add(n)

    signals = []
    for s in all_signals:
        p = paths(s, node, g)
        if p:
            signals.append(s)
        
    return signals


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
