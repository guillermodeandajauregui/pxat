import networkx as nx
from Bio.KEGG.KGML import KGML_parser
from itertools import combinations
from pprint import pprint


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


    

def targets_from_node(g, node, topological_type):
    all_ttypes = set()
    for n in g.nodes():
        if topological_type in g.node[n]:
            all_ttypes.add(n)

    targets = []
    for s in all_ttypes:
        p = paths(node, s, g)
        if p:
            targets.append(s)
        
    return targets



def sources_to_node(g, node, topological_type):
    all_ttypes = set()
    for n in g.nodes():
        if topological_type in g.node[n]:
            all_ttypes.add(n)

    sources = []
    for s in all_ttypes:
        p = paths(s, node, g)
        if p:
            sources.append(s)
        
    return sources









# a shorter way
def readKGML(kgml):
    return KGML_parser.parse(kgml).next()









# parse KGML file and return pathway object
def kgml2graph(pathway):
    """
    pathway: a KGML pathway, as parsed by BioPython's KGML_parser
    
    returns: a NetworkX direted graph

    """
    g = nx.DiGraph()

    # dismember group entries
    groups = set()
    for entry in pathway.entries.values():
        if entry.type=='group':
            for c in list(entry.components):
                # tag their nodes
                for node in pathway.entries[c.id].name.split():
                    if node.startswith('hsa'):
                        group = (pathway.title, entry.id)
                        groups.add(group)
                        if node in g.nodes():
                            g.node[node]['groups'].append(group)
                        else:
                            g.add_node( node,
                                        pathways=set([pathway.title,]),
                                        groups=[group, ])

    # connect all nodes of a group in a cluster
    for group in groups:
        # get nodes for each group
        nodes = set()
        for n in g.nodes():
            if group in g.node[n]['groups']:
                nodes.add(n)
        # pair them, connect them
        for pair in combinations(nodes, 2):
            # must check for forward edge
            if g.get_edge_data(pair[0],pair[1]):
                f_edge_groups = g.get_edge_data(pair[0],pair[1])['groups']
                f_edge_groups.add(group)
            else:
                f_edge_groups = set([group, ])
            # and reverse edge
            if g.get_edge_data(pair[1],pair[0]):
                r_edge_groups = g.get_edge_data(pair[1],pair[0])['groups']
                r_edge_groups.add(group)
            else:
                r_edge_groups = set([group, ])
                
            # both are in the same group so bunch them together
            edge_groups = f_edge_groups.union(r_edge_groups)

            # upsert in and out edges
            g.add_edge( pair[0],
                        pair[1],
                        groups=edge_groups,
                        subtypes=['group', ],
                        pathways=set([pathway.title,]))                        

            g.add_edge( pair[1],
                        pair[0],
                        groups=edge_groups,
                        subtypes=['group', ],
                        pathways=set([pathway.title,]))                        
            


            
    # copy relations to edges
    for relation in pathway.relations:
        # no undefined or path entries, only hsa names
        if relation.entry1.name.startswith('hsa') and relation.entry2.name.startswith('hsa'):
            for e1 in relation.entry1.name.split():
                for e2 in relation.entry2.name.split():
                    g.add_node( e1, pathways=set([pathway.title,]))
                    g.add_node( e2, pathways=set([pathway.title,]))

                    edge = g.get_edge_data(e1, e2)
                    if edge:
                        if edge['subtypes'] == None:
                            subtypes = list()
                        else:
                            subtypes = edge['subtypes'].append(relation.subtypes)
                    else:
                        if relation.subtypes:
                            subtypes = relation.subtypes
                        else:
                            subtypes = list()
                        
                    g.add_edge( e1,
                                e2,
                                type     = relation.type,
                                subtypes = subtypes,
                                pathways = set([pathway.title,]))
    return g






def write_adj_matrix( path, g ):
    with open(path, 'w') as out:
        a = nx.adjacency_matrix(g).toarray()

        nodes = g.nodes()
        
        out.write("\t".join(['x',] + [n for n in nodes] ) + "\n")
    
        for row in a:
            line = [nodes.pop(0),] + list([str(int(i)) for i in row])
            out.write("\t".join(line) + "\n")
            
