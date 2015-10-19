import networkx as nx
import pxat


# run tests like
# $ nosetests


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
