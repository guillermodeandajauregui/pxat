import networkx as nx
from Bio.KEGG.KGML import KGML_parser
from itertools import combinations

def kgml_file_to_digraph( kgml_file):
    filepath = file(kgml_file)
    gi       = readKGML(filepath)
    kegg_pw  = kgml2graph(gi)
    return(kegg_pw)




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
