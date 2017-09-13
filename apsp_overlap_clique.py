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
    return G

def get_isolated_cliques(cliques):
    clique_nodes = {}
    isolated_cliques = []
    for clique in cliques:
        flag = 0
        for N in clique:
            if N in clique_nodes:
                flag = 1
                break
        if not flag:
            isolated_cliques.append(clique)
            for N in clique:
                clique_nodes[N] = 1
    return isolated_cliques

def get_linked_cliques(cliques):
    clique_nodes = {}
    G=nx.Graph()
    for i in range(len(cliques)):
        clique = cliques[i]
        for N in clique:
            if not N in clique_nodes:
                clique_nodes[N] = [i]
            else:
                clique_nodes[N].append(i)
        G.add_node(i)

    edges = [x for x in clique_nodes.values() if len(x)>1]
    for edge in edges:
        for m in range(len(edge)-1):
            for n in range(m+1, len(edge)):
                G.add_edge(edge[m],edge[n])

    #pdb.set_trace()
    linked_cliques = []
    subgraphs = nx.connected_component_subgraphs(G)
    for subgraph in subgraphs:
        new_clique = []
        for clique_idx in subgraph:
            new_clique.extend(cliques[clique_idx])
        #if len(subgraph)==5:
        linked_cliques.append(list(set(new_clique)))
    return linked_cliques

def merge_linked_cliques(G, G_un, read_node_map, read_db, clique_cutoff):
    max_cliques = list(nx.find_cliques(G_un))
    max_cliques = [clique for clique in max_cliques if len(clique)>=clique_cutoff]

    linked_cliques = get_linked_cliques(max_cliques)
    while linked_cliques:
        #pdb.set_trace()
        for clique in linked_cliques:
            print len(clique)
            G_clique = G.subgraph(clique)
            DFS_transitive_reduction(G_clique)

            start_nodes = [x for x in G_clique if G_clique.in_degree(x)==0]
            end_nodes = [n for n in G_clique if G_clique.out_degree(n)==0]
            
            DFS_collapse_graph(G_clique, read_node_map, read_db) # After collapse the clique, only one node should be left
            # add the clique graph to the original graph
            for this_edge in G_clique.edges():
                n1, n2= this_edge[0], this_edge[1]
                if not G.has_edge(this_edge[0], this_edge[1]):
                    if not n1 in G:
                        N_read_ids = G_clique.node[n1]['read_ids']
                        G.add_node(n1, read_ids = N_read_ids)
                    if not n2 in G:
                        N_read_ids = G_clique.node[n2]['read_ids']
                        G.add_node(n2, read_ids = N_read_ids)
                    o = G_clique[this_edge[0]][this_edge[1]][0]['label']
                    G.add_edge(this_edge[0], this_edge[1], label = o)
            
            for start_node in start_nodes:
                start_read_id = start_node.split('|')[0]
                new_start_node = read_node_map[start_read_id] # collapsed node
                start_pred_nodes = G.predecessors(start_node) 
                new_read_ids = G_clique.node[new_start_node]['read_ids']
                G.add_node(new_start_node, read_ids=new_read_ids)
                G_un.add_node(new_start_node)
                for pred_node in start_pred_nodes:
                    o = G[pred_node][start_node][0]['label']
                    if not new_start_node in G[pred_node]:
                        G.add_edge(pred_node, new_start_node, label = o)
                        G_un.add_edge(pred_node, new_start_node, label = o)
            
            for end_node in end_nodes:
                end_read_id = end_node.split('|')[0]
                new_end_node = read_node_map[end_read_id] # collapsed node
                end_succ_nodes = G.successors(end_node)
                new_read_ids = G_clique.node[new_end_node]['read_ids']
                G.add_node(new_end_node, read_ids = new_read_ids)
                G_un.add_node(new_end_node)
                for succ_node in end_succ_nodes:
                    o = G[end_node][succ_node][0]['label']
                    if not succ_node in G[new_end_node]:
                        G.add_edge(new_end_node, succ_node, label = o)
                        G_un.add_edge(new_end_node, succ_node, label = o)
            
            ## clean up
            for clique_node in clique:
                if not clique_node in G_clique:
                    G.remove_node(clique_node)
                    G_un.remove_node(clique_node)

        #pdb.set_trace()
        max_cliques = list(nx.find_cliques(G_un))
        max_cliques = [clique for clique in max_cliques if len(clique)>clique_cutoff]
        linked_cliques = get_linked_cliques(max_cliques)

def merge_linked_cliques_2(G, G_un, read_node_map, read_db, clique_cutoff):
    max_cliques = list(nx.find_cliques(G_un))
    max_cliques = [clique for clique in max_cliques if len(clique)>clique_cutoff]
    
    while max_cliques:
        G_all_clique = nx.MultiDiGraph()
        for clique in max_cliques:
            print len(clique)
            GC_single = G.subgraph(clique)
            DFS_transitive_reduction(GC_single)
            for this_edge in GC_single.edges():
                n1, n2 = this_edge[0], this_edge[1]
                if not n1 in G_all_clique:
                    G_all_clique.add_node(n1, read_ids = G.node[n1]['read_ids'])
                if not n2 in G_all_clique:
                    G_all_clique.add_node(n2, read_ids = G.node[n2]['read_ids'])
                if not G_all_clique.has_edge(n1, n2):
                    overlap = G[n1][n2][0]['label']
                    G_all_clique.add_edge(n1, n2, label=overlap)
        #DFS_transitive_reduction(G_all_clique)
        print "Transitive reduction of all cliques finished!", len(G_all_clique), len(G_all_clique.edges())

        subgraphs = nx.weakly_connected_components(G_all_clique)
        for subgraph in subgraphs:
            print len(subgraph)
            G_clique = G_all_clique.subgraph(subgraph)

            start_nodes = [x for x in G_clique if G_clique.in_degree(x)==0]
            end_nodes = [n for n in G_clique if G_clique.out_degree(n)==0]

            DFS_collapse_graph(G_clique, read_node_map, read_db)
            DFS_transitive_reduction(G_clique)
            DFS_collapse_graph(G_clique, read_node_map, read_db)


            # add the clique graph to the original graph
            for this_edge in G_clique.edges():
                n1, n2= this_edge[0], this_edge[1]
                if not G.has_edge(this_edge[0], this_edge[1]):
                    if not n1 in G:
                        N_read_ids = G_clique.node[n1]['read_ids']
                        G.add_node(n1, read_ids = N_read_ids)
                    if not n2 in G:
                        N_read_ids = G_clique.node[n2]['read_ids']
                        G.add_node(n2, read_ids = N_read_ids)
                    o = G_clique[this_edge[0]][this_edge[1]][0]['label']
                    G.add_edge(this_edge[0], this_edge[1], label = o)
            
            for start_node in start_nodes:
                start_read_id = start_node.split('|')[0]
                new_start_node = read_node_map[start_read_id] # collapsed node
                start_pred_nodes = G.predecessors(start_node) 
                new_read_ids = G_clique.node[new_start_node]['read_ids']
                G.add_node(new_start_node, read_ids=new_read_ids)
                G_un.add_node(new_start_node)
                for pred_node in start_pred_nodes:
                    o = G[pred_node][start_node][0]['label']
                    if not new_start_node in G[pred_node]:
                        G.add_edge(pred_node, new_start_node, label = o)
                        G_un.add_edge(pred_node, new_start_node, label = o)
            
            for end_node in end_nodes:
                end_read_id = end_node.split('|')[0]
                new_end_node = read_node_map[end_read_id] # collapsed node
                end_succ_nodes = G.successors(end_node)
                new_read_ids = G_clique.node[new_end_node]['read_ids']
                G.add_node(new_end_node, read_ids = new_read_ids)
                G_un.add_node(new_end_node)
                for succ_node in end_succ_nodes:
                    o = G[end_node][succ_node][0]['label']
                    if not succ_node in G[new_end_node]:
                        G.add_edge(new_end_node, succ_node, label = o)
                        G_un.add_edge(new_end_node, succ_node, label = o)
            
            ## clean up
            for clique_node in subgraph:
                if not clique_node in G_clique:
                    G.remove_node(clique_node)
                    G_un.remove_node(clique_node)

        max_cliques = list(nx.find_cliques(G_un))
        max_cliques = [clique for clique in max_cliques if len(clique)>clique_cutoff]

def transitive_reduction(G):
    TR = nx.MultiDiGraph()
    TR.add_nodes_from(G.nodes())
    for u in G:
        u_edges = set(G[u])
        for v in G[u]:
            u_edges -= {y for x, y in nx.dfs_edges(G, v)}
        TR.add_edges_from((u,v) for v in u_edges)
    return TR

def merge_linked_cliques_3(G, G_un, read_node_map, read_db, clique_cutoff):
    max_cliques = list(nx.find_cliques(G_un))
    max_cliques = [clique for clique in max_cliques if len(clique)>clique_cutoff]
    
    while max_cliques:
        #pdb.set_trace()
        G_all_clique = nx.MultiDiGraph()
        idx = 0
        for clique in max_cliques:
            idx+=1
            print idx, len(clique)
            GC_single = G.subgraph(clique)
            #DFS_transitive_reduction(GC_single)
            GC_single = transitive_reduction(GC_single)
            for this_edge in GC_single.edges():
                n1, n2 = this_edge[0], this_edge[1]
                if not n1 in G_all_clique:
                    G_all_clique.add_node(n1, read_ids = G.node[n1]['read_ids'])
                if not n2 in G_all_clique:
                    G_all_clique.add_node(n2, read_ids = G.node[n2]['read_ids'])
                if not G_all_clique.has_edge(n1, n2):
                    overlap = G[n1][n2][0]['label']
                    G_all_clique.add_edge(n1, n2, label=overlap)
            for N in GC_single.nodes():
                if not N in G_all_clique:
                    G_all_clique.add_node(N, read_ids = G.node[N]['read_ids'])

        #DFS_transitive_reduction(G_all_clique)
        print "Transitive reduction of all cliques finished!", len(G_all_clique.nodes()), len(G_all_clique.edges())
        #pdb.set_trace()

        subgraphs = nx.weakly_connected_components(G_all_clique)
        idx = 0
        for subgraph in subgraphs:
            idx +=1
            print idx,len(subgraph)
            G_clique = G_all_clique.subgraph(subgraph)

            start_nodes = [x for x in G_clique if G_clique.in_degree(x)==0]
            end_nodes = [n for n in G_clique if G_clique.out_degree(n)==0]

            start_pred_nodes = {}
            start_overlaps = {}
            for start_node in start_nodes:
                start_pred_nodes[start_node] = G.predecessors(start_node)
                for pred_node in start_pred_nodes[start_node]:
                    if not pred_node in G or not start_node in G[pred_node]:
                        pdb.set_trace()
                    o = G[pred_node][start_node][0]['label']
                    start_overlaps[(pred_node, start_node)] = o
            #start_pred_nodes = list(set(start_pred_nodes))
            
            end_succ_nodes = {}
            end_overlaps = {}
            for end_node in end_nodes:
                end_succ_nodes[end_node] = G.successors(end_node)
                for succ_node in end_succ_nodes[end_node]:
                    o = G[end_node][succ_node][0]['label']
                    end_overlaps[(end_node, succ_node)] = o
            #end_succ_nodes = list(set(end_succ_nodes))

            # remove the original nodes in graph
            for N in G_clique:
                if N in G:
                    G.remove_node(N)
                    G_un.remove_node(N)
            
            # simplify the linked cliques graph
            DFS_collapse_graph(G_clique, read_node_map, read_db)
            DFS_transitive_reduction(G_clique)
            DFS_collapse_graph(G_clique, read_node_map, read_db)

            # add the clique graph to the original graph
            for this_edge in G_clique.edges():
                n1, n2= this_edge[0], this_edge[1]
                if not G.has_edge(this_edge[0], this_edge[1]):
                    if not n1 in G:
                        N_read_ids = G_clique.node[n1]['read_ids']
                        G.add_node(n1, read_ids = N_read_ids)
                    if not n2 in G:
                        N_read_ids = G_clique.node[n2]['read_ids']
                        G.add_node(n2, read_ids = N_read_ids)
                    o = G_clique[this_edge[0]][this_edge[1]][0]['label']
                    G.add_edge(this_edge[0], this_edge[1], label = o)
                    G_un.add_edge(this_edge[0], this_edge[1], label = o)

            for N in G_clique.nodes():
                if not N in G:
                    N_read_ids = G_clique.node[N]['read_ids']
                    G.add_node(N, read_ids = N_read_ids)
                    G_un.add_node(N)

            
            for start_node in start_nodes:
                start_read_id = start_node.split('|')[0]
                new_start_node = read_node_map[start_read_id] # collapsed node
                pred_nodes = start_pred_nodes[start_node]
                for pred_node in pred_nodes:
                    if pred_node in G:
                        if not G.has_edge(pred_node, start_node):
                            o = start_overlaps[(pred_node, start_node)]
                            G.add_edge(pred_node, new_start_node, label = o)
                            G_un.add_edge(pred_node, new_start_node, label = o)
            
            for end_node in end_nodes:
                end_read_id = end_node.split('|')[0]
                new_end_node = read_node_map[end_read_id] # collapsed node
                succ_nodes = end_succ_nodes[end_node]
                for succ_node in succ_nodes:
                    if succ_node in G:
                        if not G.has_edge(end_node, succ_node):
                            o = end_overlaps[(end_node, succ_node)]
                            G.add_edge(new_end_node, succ_node, label = o)
                            G_un.add_edge(new_end_node, succ_node, lable = o)
            
        print len(G_un.nodes()),len(G_un.edges())
    
        max_cliques = list(nx.find_cliques(G_un))
        max_cliques = [clique for clique in max_cliques if len(clique)>clique_cutoff]

def merge_isolated_cliques(G, G_un, read_node_map, read_db):
    max_cliques = list(nx.find_cliques(G_un))
    max_cliques = [clique for clique in max_cliques if len(clique)>3]
    isolated_cliques = get_isolated_cliques(max_cliques)
    while isolated_cliques:
        for clique in isolated_cliques:
            print len(clique)
            G_clique = G.subgraph(clique)
            DFS_transitive_reduction(G_clique)

            start_nodes = [x for x in G_clique if G_clique.in_degree(x)==0]
            start_node = start_nodes[0]
            end_nodes = [n for n in G_clique if G_clique.out_degree(n)==0]
            end_node = end_nodes[0]
            
            DFS_collapse_graph(G_clique, read_node_map, read_db) # After collapse the clique, only one node should be left
            if len(G_clique)>1:
                pdb.set_trace()
            
            combined_node = G_clique.nodes()[0]
            combined_read_ids = G_clique.node[combined_node]['read_ids']
            G.add_node(combined_node, read_ids = combined_read_ids)
            G_un.add_node(combined_node)

            start_pred_nodes = G.predecessors(start_node)
            for pred_node in start_pred_nodes:
                o = G[pred_node][start_node][0]['label']
                G.add_edge(pred_node, combined_node, label = o)
                G_un.add_edge(pred_node, combined_node, label = o)

            end_succ_nodes = G.successors(end_node)
            for succ_node in end_succ_nodes:
                o = G[end_node][succ_node][0]['label']
                G.add_edge(combined_node, succ_node, label = o)
                G_un.add_edge(combined_node, succ_node, label = o)
            
            ## clean up
            for clique_node in clique:
                G.remove_node(clique_node)
                G_un.remove_node(clique_node)

        max_cliques = list(nx.find_cliques(G_un))
        max_cliques = [clique for clique in max_cliques if len(clique)>3]
        isolated_cliques = get_isolated_cliques(max_cliques)

def merge_cliques(G, G_un, read_node_map, read_db):
    max_cliques = list(nx.find_cliques(G_un))
    max_cliques = [clique for clique in max_cliques if len(clique)>3]
    while max_cliques:
        clique = max_cliques[0]
        print len(clique)
        G_clique = G.subgraph(clique)
        DFS_transitive_reduction(G_clique)

        start_nodes = [x for x in G_clique if G_clique.in_degree(x)==0]
        start_node = start_nodes[0]
        end_nodes = [n for n in G_clique if G_clique.out_degree(n)==0]
        end_node = end_nodes[0]
        
        #pdb.set_trace()
        DFS_collapse_graph(G_clique, read_node_map, read_db) # After collapse the clique, only one node should be left
        if len(G_clique)>1:
            pdb.set_trace()
        
        combined_node = G_clique.nodes()[0]
        combined_read_ids = G_clique.node[combined_node]['read_ids']
        G.add_node(combined_node, read_ids = combined_read_ids)
        G_un.add_node(combined_node)

        start_pred_nodes = G.predecessors(start_node)
        for pred_node in start_pred_nodes:
            o = G[pred_node][start_node][0]['label']
            G.add_edge(pred_node, combined_node, label = o)
            G_un.add_edge(pred_node, combined_node, label = o)

        end_succ_nodes = G.successors(end_node)
        for succ_node in end_succ_nodes:
            o = G[end_node][succ_node][0]['label']
            G.add_edge(combined_node, succ_node, label = o)
            G_un.add_edge(combined_node, succ_node, label = o)
        
        ## clean up
        for clique_node in clique:
            G.remove_node(clique_node)
            G_un.remove_node(clique_node)

        max_cliques = list(nx.find_cliques(G_un))
        max_cliques = [clique for clique in max_cliques if len(clique)>3]

def assemble_path_contig(G, PE_G, read_db, Fragment_len, path_file='HIV_paths.txt', contig_file='HIV_contigs.fa'):
    f_path = open(path_file,'w')
    f_contig = open(contig_file,'w')

    all_paths=[]
    assembled_nodes=set([])
    starting_nodes=[n for n in G.nodes() if G.in_degree(n)==0]
    for start_node in starting_nodes:
        print "Begin a new start node:",start_node
        paths=list(DFS_paths_single_pair_end(G, start_node, PE_G, read_db, Fragment_len))
        all_paths.extend(paths)
        print len(all_paths)
    
    for path in all_paths:
        out_path="--".join(path)
        f_path.write(out_path+'\n')
        assembled_nodes=assembled_nodes.union(set(path))

    contig_index=0
    contigs=get_assemblie3(G, all_paths, read_db)  # a dictionary, key: path information, value: assembled sequence
    for contig_key in contigs:
        title='>contig'+'_'+str(len(contigs[contig_key]))+'_'+str(contig_index)
        f_contig.write(title+'\n'+contigs[contig_key]+'\n')
        contig_index+=1
   
    ## Handle the unassembled nodes
    paths_unassembled=[]
    unassembled_nodes=set(G.nodes())-assembled_nodes
    print "Unassembled nodes!", len(unassembled_nodes)
    old_un_num = len(unassembled_nodes)
    new_un_num = 0
    loop = 0
    while (len(unassembled_nodes)>0 and old_un_num!=new_un_num):
        old_un_num = len(unassembled_nodes)
        sub_sub_graph=G.subgraph(list(unassembled_nodes))

        new_starting_nodes=[n for n in sub_sub_graph.nodes() if sub_sub_graph.in_degree(n)==0]
        for start_node in new_starting_nodes:
            print start_node
            paths=list(DFS_paths_single_pair_end_unassembled(G, start_node, PE_G, read_db, Fragment_len))
            paths_unassembled.extend(paths)
            for path in paths:
                out_path="--".join(path)
                f_path.write(out_path+'\n')
                assembled_nodes=assembled_nodes.union(set(path))
        unassembled_nodes=set(G.nodes())-assembled_nodes
        new_un_num = len(unassembled_nodes)
        loop+=1
        
        #pdb.set_trace()
        print "Loop: %d." % loop
        print "Unassembled nodes!",len(set(G.nodes())-assembled_nodes)

    ## assemble contigs from the subgraph
    contigs=get_assemblie3(G, paths_unassembled, read_db)  # a dictionary, key: path information, value: assembled sequence
    for contig_key in contigs:
        title='>contig'+'_'+str(len(contigs[contig_key]))+'_'+str(contig_index)
        f_contig.write(title+'\n'+contigs[contig_key]+'\n')
        contig_index+=1
    f_path.close()
    f_contig.close()


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

# parameters
#read_len = 200
#Fragment_len = 450
#overlap_cutoff = 180
#binning_overlap_cutoff = 190
#tip_len_cutoff = 400

idx=0
subgraphs=nx.weakly_connected_components(G)

# saved files
f_path1=open('Paths.txt','w')
f_contig1=open('Contigs.fa','w')
f_node = open('PEG_nodes_sequences.fa','w')

for subgraph in subgraphs:
    start_time=datetime.now()
    G_subgraph=G.subgraph(subgraph)
    G_un_subgraph=G_undirected.subgraph(subgraph)
    print "The nodes of the subgraph is: %d; the edges of the subgraph is: %d." % (len(G_subgraph.nodes()), len(G_subgraph.edges()))

    ## merge the cliques
    #merge_isolated_cliques(G_subgraph, G_un_subgraph, read_node_map, read_db)
    merge_linked_cliques_3(G_subgraph, G_un_subgraph, read_node_map, read_db, 4)

    #pdb.set_trace()
    DFS_collapse_graph(G_subgraph, read_node_map, read_db)
    print "The nodes of the subgraph after merge cliques is: %d; edges is: %d." % (len(G_subgraph.nodes()),len(G_subgraph.edges()))
    print datetime.now()-start_time

    #pdb.set_trace()
    ## remove transitive edges
    idx+=1
    DFS_transitive_reduction(G_subgraph)
    DFS_collapse_graph(G_subgraph, read_node_map, read_db)
    print "The nodes of the subgraph after removing transitive edges is: %d; edges is: %d." % (len(G_subgraph.nodes()),len(G_subgraph.edges()))

    #plot_graph(G_subgraph, 'G_subgraph.png')
    #pdb.set_trace()
    
    ## delete low overlap edges
    for this_edge in G_subgraph.edges():
        if not 0 in G_subgraph[this_edge[0]][this_edge[1]]:
            pdb.set_trace()
        if int(G_subgraph[this_edge[0]][this_edge[1]][0]['label'])<overlap_cutoff:
            G_subgraph.remove_edge(this_edge[0], this_edge[1])
    DFS_collapse_graph(G_subgraph, read_node_map, read_db)
    print "Collapsed graph after deleting low overlap edges, nodes number: %d, edges number: %d." % (len(G_subgraph),len(G_subgraph.edges()))

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
        G_subgraph=collapse_graph_2(G_subgraph, read_db, read_node_map)
        new_edge_num=len(G_subgraph.edges())
        new_node_num=len(G_subgraph.nodes())
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
    #pdb.set_trace()
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
    contigs=get_assemblie3(G_subgraph, all_paths, read_db)  # a dictionary, key: path information, value: assembled sequence
    for contig_key in contigs:
        #title='>'+contig_key+'_'+str(len(contigs[contig_key]))+'_'+str(contig_index)
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
    contigs=get_assemblie3(G_subgraph, paths_unassembled, read_db)  # a dictionary, key: path information, value: assembled sequence
    for contig_key in contigs:
        #title='>'+contig_key+'_'+str(len(contigs[contig_key]))+'_'+str(contig_index)
        if len(contigs[contig_key])>Fragment_len:
            title='>contig'+'_'+str(len(contigs[contig_key]))+'_'+str(contig_index)
            f_contig1.write(title+'\n'+contigs[contig_key]+'\n')
            contig_index+=1

f_path1.close()
f_contig1.close()
f_node.close()
