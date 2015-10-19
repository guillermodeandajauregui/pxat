import networkx as nx
from KEGG import *

from pprint import pprint



def merge_pathway_subgraphs(*subgraphs):
    """
    Merges subgraphs. Merges 'pathways' fields of nodes and edges.
    """
    
    # new graph
    h = nx.DiGraph()

    for g in subgraphs:
        for n in g.nodes():
            if n in h.node:
                h.node[n]['pathways'].update(g.node[n]['pathways'])
            else:
                h.add_node(n, pathways=g.node[n]['pathways'])
    
        for e in g.edges():
            edge_data = h.get_edge_data(*e)
            if edge_data:
                edge_data['pathways'].update(g.get_edge_data(*e)['pathways'])
                h.add_edge( *e, attr_dict=edge_data)
            else:
                h.add_edge( *e, attr_dict=g.get_edge_data(*e))

    return h








def topological_annotate(g):
    """
    Will add topological tags to each node, depending on where it occurs in a graph.
    - signal
    - receptor: recieve signals
    - transducer: relays signals
    - final transducer: one hop before an effector
    - effector: exit nodes
    """
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






def topo_subgrapher(g, nbunch):
    """
    Returns a topologically re-annotated subgraph from graph G containing
    nbunch.
    """
    return topological_annotate(nx.DiGraph(g.subgraph(nbunch)))




def subgraph_from_pathways(pathways, g):
    """
    returns a topological annotated graph containing given pathways
    """
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




    
def nodes_by_subtype(g, type):
    """
    return nbunch of nodes of a given topological type
    """
    types = []
    for n in g.nodes():
        if type in g.node[n]:
            types.append(n)
    return types





def disubgraph(g, nbunch):
    """
    Returns a subgraph of nbunch in g
    """
    return nx.DiGraph(g.subgraph(nbunch))






def crosstalk_nodes(A, B):
    """
    returns a list of nodes shared by pathways A and B
    A and B are graphs
    """
    return set(A.node).intersection(set(B.node))




def pathway_xtalk( g, ):
    """
    returns pathways that crosstalk with given pathway
    """
    


def trajectories_in_graph(g):
    """
    all trajectories from signal to effector in graph
    """
    signals = nodes_by_subtype(g, 'signal')
    effectors = nodes_by_subtype(g, 'effector')
    trajectories = []

    for n in signals:
        for m in effectors:
            for t in nx.all_simple_paths(g, n, m):
                trajectories.append( t )

    return trajectories




def trajectories_from_nbunch(g, nbunch):
    """
    all trajectories in a subgraph containing nodes from nbunch 
    """
    subnetwork = nx.DiGraph(g.subgraph(nbunch))
    subnetwork = topological_annotate(subnetwork)
    return trajectories_in_graph(subnetwork)




def is_single_component(g, nbunch):
    """
    TRUE/FALSE, whether the subgraph of nbunch is connected in a single components
    """
    subnetwork = disubgraph(g, nbunch)
    subnetwork = subnetwork.to_undirected()
    if len(list(nx.connected_components(subnetwork))) > 1 :
      return False
    else:
      return True





def reduce_crosstalk_nodes(g, pw1, pw2):
    """
    shared subtrajectories formed by crosstalk nodes between two pathways may be reduced to a "single reduced node" 
    """
    crosstalkers = crosstalk_nodes(pw1, pw2)
    g_prime = topo_subgrapher(g, crosstalkers)
    subtrajectories = trajectories_in_graph(g_prime)
    return subtrajectories
    # h = copy g
    #for every subtrajectory in subtrajectories:
    #create one reduced node, name = node_1 + node_2 + ... node_n
    #all inputs to node_1 are inputs to reduced node 
    #all outputs of node_n are outputs of reduced node
    #remove nodes 1:n
    #return h, a new graph with reduced crosstalk nodes 





def nodes_from_pathway(p, g):
    """
    return nbunch given a pathway name
    """
    nodes = []
    for n in g.nodes():
        if p in g.node[n]['pathways']:
            nodes.append(n)
    return nodes







def paths(source, target, g):
    """
    sort of an alias to all_simple_paths
    """
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











# utility function to export for cytoscape
def write_adj_matrix( path, g ):
    with open(path, 'w') as out:
        a = nx.adjacency_matrix(g).toarray()

        nodes = g.nodes()
        
        out.write("\t".join(['x',] + [n for n in nodes] ) + "\n")
    
        for row in a:
            line = [nodes.pop(0),] + list([str(int(i)) for i in row])
            out.write("\t".join(line) + "\n")
            
