import networkx as nx
import pxat



# run tests like
# $ nosetests




# test KEGG submodule
def test_kgml_file_to_digraph():
    g = pxat.KEGG.kgml_file_to_digraph('hsa04210.xml')
    
    # hsa04210 is apoptosis

    # nodes should be genes in that pathway, should come from KGML

    # edges are entries in KGML file

    # pathway name should be apoptosis

    

def test_kgml2graph():
    pass




# test JIN submodule
def test_pathway_jin():
    pass

def test_graph_jin():
    pass

def test_subgraph_jin():
    pass




# test main module
def test_merge_pathway_subgraphs():
    pass

def test_topological_annotate():
    pass

def test_topo_subgrapher():
    pass

def test_subgraph_from_pathways():
    pass

def test_nodes_by_subtype():
    pass

def test_disubgraph():
    pass

def test_crosstalk_nodes():
    pass

def test_pathway_xtalk():
    pass

def test_trajectories_in_graph():
    pass

def test_trajectories_from_nbunch():
    pass

def test_is_single_component():
    pass

def test_reduce_crosstalk_nodes():
    pass

def test_nodes_from_pathway():
    pass

def test_paths():
    pass

def test_targets_from_node():
    pass

def test_sources_to_node():
    pass

def test_write_adj_matrix():
    pass








def test_b():
    assert 'b' == 'b'

   

g = nx.gpickle.read_gpickle('red_prueba.pickle') 








#
# seek trajectories in human kegg, prototest
#

# g = nx.gpickle.read_gpickle('human_kegg_ann.pickle')

# signals   = set()
# effectors = set()

# for n in g.nodes():
#     if 'signal' in g.node[n]:
#         signals.add(n)

#     if 'effector' in g.node[n]:
#         effectors.add(n)

# for s in signals:
#     for e in effectors:
#         print s, e, [p for p in nx.all_simple_paths(g, s, e)]
