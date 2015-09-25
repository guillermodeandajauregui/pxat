from pyveplot import *
import networkx as nx
import random
        
h = Hiveplot( 'example.svg')
signal_axis   = Axis( (700, 900), (700, 41), stroke="black", stroke_dasharray="4 2", stroke_width=1.5) 
trans_axis    = Axis( (700,900), (700 + 2938,900), stroke="yellowgreen", stroke_width=2, stroke_opacity=0.6 )
effector_axis = Axis( (700,900), (101,900), stroke="blue", stroke_opacity="0.3", stroke_width=2)
h.axes = [ signal_axis, trans_axis, effector_axis ]


g = nx.gpickle.read_gpickle('human_kegg_classified.pickle')


# create and populate sigbnal axis
signals = []
for n in g.nodes():
    if 'signal' in g.node[n]:
        signals.append(n)

delta = 1.0 / len(signals)
offset = 0.01
for n in signals:
    offset += delta
    signal_axis.add_node(Node(n), offset)


# create and populate effector axis
effectors = []
for n in g.nodes():
    if 'effector' in g.node[n]:
        effectors.append(n)

delta = 1.0 / len(effectors)
offset = 0.01
for n in effectors:
    offset += delta
    effector_axis.add_node(Node(n), offset)



trans = []
for n in g.nodes():
    if not 'effector' in g.node[n] and not 'signal' in g.node[n]:
        trans.append(n)

delta = 1.0 / len(trans)
offset = 0.01
for n in trans:
    offset += delta
    trans_axis.add_node(Node(n), offset)


# edges from axis0 to axis1
print "connecting edges"
for e in g.edges():
    if (e[0] in signal_axis.nodes) and (e[1] in trans_axis.nodes):
        h.connect(signal_axis, e[0],
                  45, # angle of invisible axis for source control points
                  trans_axis, e[1],
                  -45, # angle of invisible axis for target control points
                  stroke_width = 0.34, # pass any SVG attributes to an edge
                  stroke_opacity = 0.4,
                  stroke = 'purple',
              )
    if (e[1] in signal_axis.nodes) and (e[0] in trans_axis.nodes):
        h.connect(signal_axis, e[1], 45,
                  trans_axis, e[0], -45,
                  stroke_width = 0.34, # pass any SVG attributes to an edge
                  stroke_opacity = 0.4,
                  stroke = 'purple',
              )

        # edges from axis1 to axis2
    if (e[0] in trans_axis.nodes) and (e[1] in effector_axis.nodes):
        h.connect(trans_axis, e[0], 15,
                  effector_axis, e[1], -15,
                  stroke_width = 0.34,
                  stroke_opacity = 0.4,
                  stroke = 'green',
              )
    if (e[1] in trans_axis.nodes) and (e[0] in effector_axis.nodes):
        h.connect(trans_axis, e[1], 15,
                  effector_axis, e[0], -15,
                  stroke_width = 0.34,
                  stroke_opacity = 0.4,
                  stroke = 'green',
              )
        
    # edges from axis0 to axis2
    if (e[0] in signal_axis.nodes) and (e[1] in effector_axis.nodes):
        h.connect(signal_axis, e[0], -45,
                  effector_axis, e[1], 45,
                  stroke_width = 0.34,
                  stroke_opacity = 0.4,
                  stroke = 'red',
              )
    if (e[1] in signal_axis.nodes) and (e[0] in effector_axis.nodes):
        h.connect(signal_axis, e[1], -45,
                  effector_axis, e[0], 45,
                  stroke_width = 0.34,
                  stroke_opacity = 0.9,
                  stroke = 'red',
              )
h.save()

