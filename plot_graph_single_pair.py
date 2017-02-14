import subprocess
import pygraphviz
from Gen_graph import *
from Graph_simplify import *
from Graph_assemble import *

# plot graph


###########################################################################
des_file=sys.argv[1]
edge_file=sys.argv[2]
fa_file=sys.argv[3]

# read_map: map read_name to read_index
# read_node_map: map read_index to node
read_map_single, read_map_pair, read_single, read_pair, des_list, read_db = get_seq_from_fa_mix(fa_file, des_file) # read dictionary
#pdb.set_trace()

G, read_node_map = create_graph_mix(des_list, edge_file, len(read_map_single))
G = collapse_graph_2(G,read_db,read_node_map)
subgraphs=nx.weakly_connected_components(G)
print "Graph construction finished!"
print "Original whole graph, nodes number %d, edges number %d."%(len(G),len(G.edges()))
plot_graph(G, 'Overlap_graph.png')

read_node_map_pair = {}
for read in read_node_map.keys():
    if read in read_pair:
        read_node_map_pair[read] = read_node_map[read]
paired_end_edges,PE_G = create_paired_end_graph(read_node_map_pair)
print "Pair-end graph, nodes number: %d, edges number: %d"%(len(PE_G), len(PE_G.edges()))
plot_graph(PE_G, 'Pair-end_graph.png')

pdb.set_trace()
# remove tips and bubbles
old_edge_num=len(G.edges())
new_edge_num=0
old_node_num=len(G.nodes())
new_node_num=0
loop=0
while (old_edge_num!=new_edge_num) or (old_node_num!=new_node_num):
    old_edge_num=len(G.edges())
    old_node_num=len(G.nodes())
    remove_tips(G, read_db, 500)
    remove_super_tips(G, read_db)
    collapse_graph_2(G,read_db,read_node_map)
    new_edge_num=len(G.edges())
    new_node_num=len(G.nodes())
    loop+=1

#pdb.set_trace()
read_node_map_pair = {}
for read in read_node_map.keys():
    if read in read_pair:
        read_node_map_pair[read] = read_node_map[read]
paired_end_edges,PE_G = create_paired_end_graph(read_node_map_pair)
print "After removing tips and bubbles, nodes number: %d, number of edges: %d"%(len(G), len(G.edges()))
plot_graph(G, 'Overlap_graph_removed_tips.png')
print "Pair-end graph, nodes number: %d, edges number: %d"%(len(PE_G), len(PE_G.edges()))
plot_graph(PE_G, 'Pair-end_graph_removed_tips.png')


## delete low overlap edges
for this_edge in G.edges():
    if int(G.edge[this_edge[0]][this_edge[1]][0]['label'])<170:
        G.remove_edge(this_edge[0],this_edge[1])
G= collapse_graph_2(G, read_db, read_node_map)
print "Collapsed graph after deleting low overlap edges, nodes number: %d, edges number: %d."%(len(G),len(G.edges()))

fig_name = 'Overlap_graph_removed_tips_delete_low_overlaps.png'
plot_graph(G, fig_name)

'''
# parameters
read_len = 230
Fragment_len = 450
overlap_cutoff = 170
binning_overlap_cutoff = 180
tip_len_cutoff = 500

graph_idx = 0
for gg in subgraphs:
    subgraph_simple = G.subgraph(gg)
    print "Original graph, nodes number %d, edges number %d."%(len(subgraph_simple), len(subgraph_simple.edges()))
    subgraph_simple=collapse_graph_2(G.subgraph(gg),read_db,read_node_map)
    print "Collapsed graph, nodes number %d, edges number %d."%(len(subgraph_simple), len(subgraph_simple.edges()))
    
    ## delete low overlap edges
    for this_edge in subgraph_simple.edges():
        if int(subgraph_simple.edge[this_edge[0]][this_edge[1]][0]['label'])<overlap_cutoff:
            subgraph_simple.remove_edge(this_edge[0],this_edge[1])
    subgraph_simple = collapse_graph_2(subgraph_simple, read_db, read_node_map)
    print "Collapsed graph after deleting low overlap edges, nodes number: %d, edges number: %d."%(len(subgraph_simple),len(subgraph_simple.edges()))
    
    fig_name = 'Subgraph_'+str(graph_idx)+'.png'
    plot_graph(subgraph_simple, fig_name)
'''
