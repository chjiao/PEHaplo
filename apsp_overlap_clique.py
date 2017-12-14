import re,sys,pdb,subprocess
import networkx as nx
import pygraphviz
from datetime import datetime
from lib.Graph_simplify import *
from lib.Graph_assemble import *

def DFS_transitive_reduction(G):
    for u in G.nodes():
        for v in G.successors(u):
            visited=set()
            stack=[v]
            while stack:
                vertex=stack.pop()
                if vertex not in visited:
                    if vertex in G[u] and vertex!=v:
                        G.remove_edge(u,vertex)
                    visited.add(vertex)
                    stack.extend(set(G.successors(vertex))-visited)

def generate_sequence_file(fa_file, seq_file='sequences.txt'):
    """
    read_map, # key: read_name, value: read_index
    read_name_list, # read names
    seq_dict, # key: read_index, value, sequence
    """
    read_map={}
    read_name_list=[]
    seq_dict={}
    f_out=open(seq_file,'w')
    read_idx=0
    with open(fa_file,'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                title=line
            else:
                seq=line
                read_map[title]=str(read_idx)
                read_name_list.append(str(read_idx))
                seq_dict[str(read_idx)]=seq
                read_idx+=1
                f_out.write(line+'\n')
    f_out.close()
    return read_name_list, read_map, seq_dict

def create_graph_apsp(overlap_file):
    """ 
    create the initial graph, node name: read_index
    read_node_dict: store the corresponding node for each read, key: read_index, value: corresponding node
    Apsp only generates plus-plus overlaps
    """
    G = nx.MultiDiGraph()
    read_node_dict = {}
    with open(overlap_file,'r') as f:
        for line in f:
            reads, overlap_len=line.strip().split(' ---> ')
            read_1,read_2=reads.strip().split()
            
            if not read_1 in G:
                G.add_node(read_1,read_ids=[read_1])
            read_node_dict[read_1]=read_1
            if not read_2 in G[read_1] and (not read_1==read_2): # bug fixed, judge whether the edge already exist, remove the edge connecting to self
                G.add_node(read_2,read_ids=[read_2])
                G.add_edge(read_1,read_2,label=overlap_len)
                read_node_dict[read_2]=read_2
            elif read_2 in G[read_1]:
                print "Duplicate edge found!",line.strip()
                if int(overlap_len)>int(G[read_1][read_2][0]['label']):
                    G[read_1][read_2][0]['label']=overlap_len

    return G, read_node_dict 

def create_graph_apsp_undirected(overlap_file):
    G = nx.Graph()
    with open(overlap_file,'r') as f:
        for line in f:
            reads, overlap_len=line.strip().split(' ---> ')
            read_1,read_2=reads.strip().split()
            
            if not read_1 in G:
                G.add_node(read_1)
            if not read_2 in G[read_1] and (not read_1==read_2): # bug fixed, judge whether the edge already exist, remove the edge connecting to self
                G.add_node(read_2)
                G.add_edge(read_1,read_2,label=overlap_len)
            elif read_2 in G[read_1]:
                print "Duplicate edge found!",line.strip()
                if int(overlap_len)>int(G[read_1][read_2]['label']):
                    G[read_1][read_2]['label']=overlap_len
    return G

def read_pair_file(pair_file, read_map):
    pair_dict = {}
    with open(pair_file,'r') as f:
        for line in f:
            line=line.strip()
            pair1, pair2 = line.split()
            if (not pair1 in read_map) or (not pair2 in read_map):
                continue
            idx1,idx2 = read_map[pair1],read_map[pair2]
            pair_dict[(idx1,idx2)]=1
    return pair_dict

def DFS_collapse_graph(G, read_node_dict, read_db):
    # G: at least two nodes
    if len(G.nodes())<2:
        return G
    starting_nodes=[n for n in G.nodes() if G.in_degree(n)==0]
    visited = set()
    stack_dict = {} # make sure that a node will not be pushed in the stack twice
    for start_node in starting_nodes:
        stack = [start_node]
        stack_dict[start_node] = 1
        while stack:
            vertex=stack.pop()
            if vertex not in visited:
                if G.out_degree(vertex)==1 and G.in_degree(G.successors(vertex)[0])==1:  # collapse
                    succ_node = G.successors(vertex)[0]
                    pred_nodes=G.predecessors(vertex)
                    succ_succ_nodes=G.successors(succ_node)
                    
                    ## update read_ids
                    combined_read_ids = G.node[vertex]['read_ids']
                    combined_read_ids.extend(G.node[succ_node]['read_ids'])

                    ## update node
                    combined_node=combined_read_ids[0]+'|'+str(len(combined_read_ids)-1)

                    ## update read_node_dict
                    for read in combined_read_ids:
                        read_node_dict[read] = combined_node

                    G.add_node(combined_node, read_ids=combined_read_ids)
                    for pred_node in pred_nodes:
                        o=G[pred_node][vertex][0]["label"]
                        G.add_edge(pred_node, combined_node, label=o)
                    for succ_succ_node in succ_succ_nodes:
                        if not 0 in G[succ_node][succ_succ_node]:
                            pdb.set_trace()
                        o=G[succ_node][succ_succ_node][0]["label"]
                        G.add_edge(combined_node, succ_succ_node, label=o)

                    ## update sequences
                    overlap_len = int(G[vertex][succ_node][0]["label"])
                    vertex_seq = read_db[vertex]
                    succ_seq = read_db[succ_node]
                    combined_seq = vertex_seq + succ_seq[overlap_len:]
                    read_db[combined_node] = combined_seq

                    ## clean up
                    G.remove_node(vertex)
                    G.remove_node(succ_node)

                    del read_db[vertex]
                    del read_db[succ_node]
                    stack.append(combined_node)
                    stack_dict[combined_node] = 1
                else:
                    visited.add(vertex)
                    if not vertex in G:
                        pdb.set_trace()
                    for succ_node in list(set(G.successors(vertex)) - visited):
                        if not succ_node in stack_dict:
                            stack_dict[succ_node] = 1
                            stack.append(succ_node)
                            #stack.extend(set(G.successors(vertex))-visited)

def linear_merge_linked_cliques(G, G_un, read_node_map, read_db, clique_cutoff):
    result = nx.MultiDiGraph()
    max_cliques = list(nx.find_cliques(G_un))
    max_cliques = [clique for clique in max_cliques if len(clique)>=clique_cutoff]

    print(len(max_cliques))
    for clique in max_cliques:
        #result.add_nodes_from(clique)
        for N in clique:
            result.add_node(N, read_ids = G.node[N]['read_ids'])
        G_clique = G.subgraph(clique)
        #result.add_edges_from(G_clique.edges())
        for E in G_clique.edges():
            if not result.has_edge(E[0], E[1]):
                o = G_clique[E[0]][E[1]][0]['label']
                result.add_edge(E[0], E[1], label = o) 
    
    # link the clique graph to the original graph
    starting_nodes = [x for x in result if result.in_degree(x)==0]
    ending_nodes = [x for x in result if result.out_degree(x)==0]
    nonClique_nodes = set(G.nodes()) - set(result.nodes())
    nonClique_nodes = nonClique_nodes | set(starting_nodes)
    nonClique_nodes = nonClique_nodes | set(ending_nodes)
    nonClique_subgraph = G.subgraph(list(nonClique_nodes))
    #result.add_edges_from(nonClique_subgraph.edges())
    for E in nonClique_subgraph.edges():
        if not result.has_edge(E[0],E[1]):
            if not E[0] in result:
                result.add_node(E[0], read_ids=G.node[E[0]]['read_ids'])
            if not E[1] in result:
                result.add_node(E[1], read_ids=G.node[E[1]]['read_ids'])
            o = nonClique_subgraph[E[0]][E[1]][0]['label']
            result.add_edge(E[0], E[1], label=o)

    idx = 0
    old_num = 0
    new_num = len(result.nodes())
    while old_num != new_num:
        idx+=1
        print idx
        old_num = len(result.nodes())
        #pdb.set_trace()
        result = linear_transitive_reduction(result)
        DFS_collapse_graph(result, read_node_map, read_db)
        #pdb.set_trace()
        new_num = len(result.nodes())
    #result = linear_transitive_reduction(result)
    
    it = 1
    G_un = result.to_undirected()
    max_cliques = list(nx.find_cliques(G_un))
    max_cliques = [clique for clique in max_cliques if len(clique)>=clique_cutoff]
    #pdb.set_trace()
    while(max_cliques):
        result2 = nx.MultiDiGraph()
        print("New round:", len(max_cliques))
        for clique in max_cliques:
            #result2.add_nodes_from(clique)
            for N in clique:
                result2.add_node(N, read_ids = G.node[N]['read_ids'])
            G_clique = result.subgraph(clique)
            #result2.add_edges_from(G_clique.edges())
            for E in G_clique.edges():
                if not result2.has_edge(E[0], E[1]):
                    o = G_clique[E[0]][E[1]][0]['label']
                    result.add_edge(E[0], E[1], label = o) 
        
        # link the clique graph to the original graph
        starting_nodes = [x for x in result2 if result2.in_degree(x)==0]
        ending_nodes = [x for x in result2 if result2.out_degree(x)==0]
        nonClique_nodes = set(result.nodes()) - set(result2.nodes())
        nonClique_nodes = nonClique_nodes | set(starting_nodes)
        nonClique_nodes = nonClique_nodes | set(ending_nodes)
        nonClique_subgraph = result.subgraph(list(nonClique_nodes))
        #result2.add_edges_from(nonClique_subgraph.edges())
        for E in nonClique_subgraph.edges():
            if not result2.has_edge(E[0],E[1]):
                if not E[0] in result2:
                    result2.add_node(E[0], read_ids=result.node[E[0]]['read_ids'])
                if not E[1] in result2:
                    result2.add_node(E[1], read_ids=result.node[E[1]]['read_ids'])
                o = nonClique_subgraph[E[0]][E[1]][0]['label']
                result2.add_edge(E[0], E[1], label=o)

        idx = 0
        old_num = 0
        new_num = len(result.nodes())
        while old_num != new_num:
            idx+=1
            print "New round:", idx
            old_num = len(result2.nodes())
            result2 = linear_transitive_reduction(result2)
            DFS_collapse_graph(result2, read_node_map, read_db)
            new_num = len(result2.nodes())

        result = result2 
        G_un = result2.to_undirected()
        max_cliques = list(nx.find_cliques(G_un))
        max_cliques = [clique for clique in max_cliques if len(clique)>=clique_cutoff]
        print "%dth iteration finished!" % it
        it += 1
    """
    G_un = result.to_undirected()
    max_cliques = list(nx.find_cliques(G_un))
    max_cliques = [clique for clique in max_cliques if len(clique)>=clique_cutoff]
    pdb.set_trace()
    """
    return result


def transitive_reduction(G):
    TR = nx.MultiDiGraph()
    TR.add_nodes_from(G.nodes())
    for u in G:
        u_edges = set(G[u])
        for v in G[u]:
            u_edges -= {y for x, y in nx.dfs_edges(G, v)}
        TR.add_edges_from((u,v) for v in u_edges)
    return TR

def linear_transitive_reduction(G):
    TR = nx.MultiDiGraph()
    for N in G.nodes():
        if not 'read_ids' in G.node[N]:
            pdb.set_trace()
        TR.add_node(N, read_ids = G.node[N]['read_ids'])
    #TR.add_nodes_from(G.nodes())
    for u in G:
        u_edges = set(G[u])
        for v in G[u]:
            u_edges -= set(G[v])
        #TR.add_edges_from((u,v) for v in u_edges)
        for v in u_edges:
            o = G[u][v][0]['label']
            TR.add_edge(u, v, label = o)
    #TR = transitive_reduction(TR)
    return TR

def plot_graph(G, figname):
    G_plot=nx.drawing.nx_agraph.to_agraph(G)
    G_plot.draw(figname,prog='dot')

#######################################################
fa_file=sys.argv[1]
pair_file=sys.argv[2]
overlap = sys.argv[3]

# parameters
read_len = int(sys.argv[4])
Fragment_len = int(sys.argv[5])
overlap_cutoff = int(sys.argv[6])
binning_overlap_cutoff = overlap_cutoff+10
tip_len_cutoff = 400

des_list, read_map, read_db = generate_sequence_file(fa_file, 'sequences.txt')
subprocess.call("Apsp sequences.txt -p 4 -m %s -o 2 >overlap_whole.txt" % overlap, shell=True)
overlap_file='overlap_whole.txt'
G, read_node_map = create_graph_apsp(overlap_file)
G_undirected = create_graph_apsp_undirected(overlap_file)
print "The nodes of the whole graph is: %d; the edges of the whole graph is: %d." % (len(G.nodes()), len(G.edges()))
pair_dict = read_pair_file(pair_file, read_map)
#pdb.set_trace()

idx=0
subgraphs=nx.weakly_connected_components(G)

# saved files
f_path1=open('Paths.txt','w')
f_contig1=open('Contigs.fa','w')
f_node0 = open('PEG_nodes_sequences_before_removing.fa', 'w')
f_node = open('PEG_nodes_sequences.fa','w')

for subgraph in subgraphs:
    start_time=datetime.now()
    G_subgraph=G.subgraph(subgraph)
    G_un_subgraph=G_undirected.subgraph(subgraph)
    print "The nodes of the subgraph is: %d; the edges of the subgraph is: %d." % (len(G_subgraph.nodes()), len(G_subgraph.edges()))

    ## merge the cliques
    #merge_isolated_cliques(G_subgraph, G_un_subgraph, read_node_map, read_db)
    #merge_linked_cliques_3(G_subgraph, G_un_subgraph, read_node_map, read_db, 4)
    #pdb.set_trace()
    G_subgraph = linear_merge_linked_cliques(G_subgraph, G_un_subgraph, read_node_map, read_db, 4)
    print "The nodes of the subgraph now is: %d; edges is: %d." % (len(G_subgraph.nodes()), len(G_subgraph.edges())) 

    DFS_collapse_graph(G_subgraph, read_node_map, read_db)
    print "The nodes of the subgraph after merge cliques is: %d; edges is: %d." % (len(G_subgraph.nodes()),len(G_subgraph.edges()))
    print datetime.now()-start_time

    ## remove transitive edges
    idx+=1
    DFS_transitive_reduction(G_subgraph)
    DFS_collapse_graph(G_subgraph, read_node_map, read_db)
    print "The nodes of the subgraph after removing transitive edges is: %d; edges is: %d." % (len(G_subgraph.nodes()),len(G_subgraph.edges()))

    #pdb.set_trace()
    
    # save the graph
    #graph_save_name = 'graph_'+str(idx)+'.gpickle'
    #nx.write_gpickle(G_subgraph, graph_save_name)
    
    ## delete low overlap edges
    for this_edge in G_subgraph.edges():
        if not 0 in G_subgraph[this_edge[0]][this_edge[1]]:
            pdb.set_trace()
        if int(G_subgraph[this_edge[0]][this_edge[1]][0]['label'])<overlap_cutoff:
            G_subgraph.remove_edge(this_edge[0], this_edge[1])
    DFS_collapse_graph(G_subgraph, read_node_map, read_db)
    print "Collapsed graph after deleting low overlap edges, nodes number: %d, edges number: %d." % (len(G_subgraph),len(G_subgraph.edges()))
    
    for N in G_subgraph.nodes():
        f_node0.write('>'+N+'\n'+read_db[N]+'\n')
    
    """
    paired_end_edges, PE_G = create_paired_end_graph_with_pairs(read_node_map, pair_dict)
    for N in G_subgraph.nodes():
        if not N in PE_G:
            PE_G.add_node(N)
    """
    #"""
    ## binning graph
    old_edge_num=len(G_subgraph.edges())
    new_edge_num=0
    old_node_num=len(G_subgraph.nodes())
    new_node_num=0
    loop=0
    while old_edge_num!=new_edge_num or old_node_num!=new_node_num:
        loop+=1
        old_edge_num=len(G_subgraph.edges())
        old_node_num=len(G_subgraph.nodes())

        paired_end_edges, PE_G = create_paired_end_graph_with_pairs(read_node_map, pair_dict)
        for N in G_subgraph.nodes():
            if not N in PE_G:
                PE_G.add_node(N)

        pair_end_binning(G_subgraph, PE_G, binning_overlap_cutoff)  # delete edges with no pair-end supporting
        #G_subgraph=collapse_graph_2(G_subgraph, read_db, read_node_map)
        #DFS_collapse_graph(G_subgraph, read_node_map, read_db)
        new_edge_num=len(G_subgraph.edges())
        new_node_num=len(G_subgraph.nodes())
        print "After binning, the number of nodes: %d, number of edges: %d."%(new_node_num, new_edge_num)
    print 'Loop:',loop
    print "After binning, the number of nodes: %d, number of edges: %d."%(new_node_num, new_edge_num)
    
    # remove tips and bubbles
    old_edge_num=len(G_subgraph.edges())
    new_edge_num=0
    old_node_num=len(G_subgraph.nodes())
    new_node_num=0
    loop=0
    while (old_edge_num!=new_edge_num) or (old_node_num!=new_node_num):
        old_edge_num=len(G_subgraph.edges())
        old_node_num=len(G_subgraph.nodes())
        remove_tips(G_subgraph, read_db, tip_len_cutoff)
        remove_super_tips(G_subgraph, read_db)
        collapse_graph_2(G_subgraph,read_db,read_node_map)
        new_edge_num=len(G_subgraph.edges())
        new_node_num=len(G_subgraph.nodes())
        loop+=1

    paired_end_edges,PE_G = create_paired_end_graph_with_pairs(read_node_map, pair_dict)
    for N in G_subgraph.nodes():
        if not N in PE_G:
            PE_G.add_node(N)
    print "After removing tips and bubbles, nodes number: %d, number of edges: %d" % (len(G_subgraph), len(G_subgraph.edges()))
    #plot_graph(G_subgraph, "HIV_graph.png") 
    #pdb.set_trace()
    #"""

    for N in G_subgraph.nodes():
        f_node.write('>'+N+'\n'+read_db[N]+'\n')

    ## find path and assembly
    #pdb.set_trace()
    all_paths=[]
    assembled_nodes=set([])
    starting_nodes=[n for n in G_subgraph.nodes() if G_subgraph.in_degree(n)==0]
    for start_node in starting_nodes:
        print "Begin a new start node:",start_node
        paths=list(DFS_paths_single_pair_end(G_subgraph, start_node, PE_G, read_db, Fragment_len))
        all_paths.extend(paths)
        print len(all_paths)
    
    for path in all_paths:
        out_path="--".join(path)
        f_path1.write(out_path+'\n')
        f_path1.flush()
        assembled_nodes=assembled_nodes.union(set(path))

    contig_index=0
    contigs=get_assemblie(G_subgraph, all_paths, read_db)  # a dictionary, key: path information, value: assembled sequence
    for contig_key in contigs:
        if len(contigs[contig_key])>=Fragment_len:
            title='>contig'+'_'+str(len(contigs[contig_key]))+'_'+str(contig_index)
            f_contig1.write(title+'\n'+contigs[contig_key]+'\n')
            contig_index+=1
            f_contig1.flush()
   
    #pdb.set_trace()
    paths_unassembled=[]
    ## Handle the unassembled nodes
    unassembled_nodes=set(G_subgraph.nodes())-assembled_nodes
    print "Unassembled nodes!", len(unassembled_nodes)
    old_un_num = len(unassembled_nodes)
    new_un_num = 0
    loop = 0
    while (len(unassembled_nodes)>0 and old_un_num!=new_un_num):
        old_un_num = len(unassembled_nodes)
        sub_sub_graph=G_subgraph.subgraph(list(unassembled_nodes))

        new_starting_nodes=[n for n in sub_sub_graph.nodes() if sub_sub_graph.in_degree(n)==0]
        for start_node in new_starting_nodes:
            print start_node
            paths=list(DFS_paths_single_pair_end_unassembled(G_subgraph, start_node, PE_G, read_db, Fragment_len))
            paths_unassembled.extend(paths)
            for path in paths:
                out_path="--".join(path)
                f_path1.write(out_path+'\n')
                f_path1.flush()
                assembled_nodes=assembled_nodes.union(set(path))
        unassembled_nodes=set(G_subgraph.nodes())-assembled_nodes
        new_un_num = len(unassembled_nodes)
        loop+=1
        
        #pdb.set_trace()
        print "Loop: %d." % loop
        print "Unassembled nodes!",len(set(G_subgraph.nodes())-assembled_nodes)

    ## assemble contigs from the subgraph
    contigs=get_assemblie(G_subgraph, paths_unassembled, read_db)  # a dictionary, key: path information, value: assembled sequence
    for contig_key in contigs:
        if len(contigs[contig_key])>Fragment_len:
            title='>contig'+'_'+str(len(contigs[contig_key]))+'_'+str(contig_index)
            f_contig1.write(title+'\n'+contigs[contig_key]+'\n')
            contig_index+=1

f_path1.close()
f_contig1.close()
f_node0.close()
f_node.close()
