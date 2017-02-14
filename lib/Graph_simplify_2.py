import re,sys,pdb
import networkx as nx

def collapse_graph(G, candidates,read_db, read_node_dict):
    # node collapsed: combined node
    # read_node_dict: store the corresponding node for each read, key: read_base, value: [.1 node, .2 node]
    while True:
        nodes_to_combine = []
        if not candidates:
            all_node = G.nodes()
        else:
            all_node = candidates

        for node in all_node:
            #pdb.set_trace()
            if G.in_degree(node) == 1 and G.out_degree(G.predecessors(node)[0]) == 1:
                predecessor = G.predecessors(node)[0]
                if compare_list(G.node[node]['species'],G.node[predecessor]['species']):
                    nodes_to_combine.append(node)
                    if candidates:
                        candidates.remove(node)

        if not nodes_to_combine:
            break

        for node_to_combine in nodes_to_combine:
            if G.in_degree(node_to_combine)==0:
                pdb.set_trace()

            predecessor = G.predecessors(node_to_combine)[0]
            predecessors_predecessors = G.predecessors(predecessor)
            successors = G.successors(node_to_combine)

            # update graph
            ## update species
            pre_species=G.node[predecessor]['species']
            node_to_combine_species=G.node[node_to_combine]['species']
            
            if compare_list(pre_species,node_to_combine_species):
                combined_species=join_species(pre_species,node_to_combine_species)
            else:
                continue

            ## update read id
            combined_read_ids=G.node[predecessor]['read_ids']
            #if pre_read_ids==None:
            #    pdb.set_trace()
            combined_read_ids.extend(G.node[node_to_combine]['read_ids'])

            ## update node
            if node_to_combine.find('|')>-1 and predecessor.find('|')==-1:
                combined_node = predecessor+'|'+str(int(node_to_combine.split('|')[1])+1)

            elif node_to_combine.find('|')>-1 and predecessor.find('|')>-1:
                combined_node = predecessor.split('|')[0]+'|'+str(int(predecessor.split('|')[1])+int(node_to_combine.split('|')[1])+1)
            elif predecessor.find('|')>-1 and node_to_combine.find('|')==-1:
                combined_node = predecessor.split('|')[0]+'|'+str(int(predecessor.split('|')[1])+1)
            else:
                combined_node = predecessor +'|'+str(1)
            #pdb.set_trace()
           
            ## update read_node_dict
            for read in combined_read_ids:
                if read.endswith('1'):
                    read_node_dict[read.split('/')[0]][0]=combined_node
                else:
                    read_node_dict[read.split('/')[0]][1]=combined_node

            if len(combined_species)==1:
                color_node=color_dict[combined_species[0]]
            else:
                color_node='black'
            G.add_node(combined_node,species=combined_species,read_ids=combined_read_ids,color=color_node,penwidth='2.0')
            for predecessors_predecessor in predecessors_predecessors:
                o = G[predecessors_predecessor][predecessor][0]["label"]
                G.add_edge(predecessors_predecessor, combined_node, label= o, color='blue', penwidth='2.0')

            for successor in successors:
                o = G[node_to_combine][successor][0]["label"]
                G.add_edge(combined_node, successor, label= o, color='blue', penwidth='2.0')

            # update sequences
            overlap_to_predecessor = int(G[predecessor][node_to_combine][0]["label"])
            offset = len(read_db[predecessor]) - overlap_to_predecessor

            pred_seq = read_db[predecessor]
            node_seq = read_db[node_to_combine]
            combined_seq = pred_seq + node_seq[overlap_to_predecessor:]

            read_db[combined_node] = combined_seq

            # clean up
            G.remove_node(node_to_combine)
            G.remove_node(predecessor)

            del read_db[node_to_combine]
            del read_db[predecessor]

            if node_to_combine in nodes_to_combine:
                nodes_to_combine.remove(node_to_combine)
            if predecessor in nodes_to_combine:
                nodes_to_combine.remove(predecessor)
    return G

def collapse_graph_2(G, read_db, read_node_dict):
    # node collapsed: combined node
    # read_node_dict: store the corresponding node for each read, key: read_base, value: [.1 node, .2 node]
    while True:
        nodes_to_combine = []
        all_node = G.nodes()

        for node in all_node:
            #pdb.set_trace()
            if G.in_degree(node) == 1 and G.out_degree(G.predecessors(node)[0]) == 1:
                predecessor = G.predecessors(node)[0]
                nodes_to_combine.append(node)

        if not nodes_to_combine:
            break

        for node_to_combine in nodes_to_combine:
            predecessor = G.predecessors(node_to_combine)[0]
            predecessors_predecessors = G.predecessors(predecessor)
            successors = G.successors(node_to_combine)
            
            ## update read_ds
            combined_read_ids=G.node[predecessor]['read_ids'] 
            combined_read_ids.extend(G.node[node_to_combine]['read_ids'])

            ## update node
            combined_node=combined_read_ids[0]+'|'+str(len(combined_read_ids)-1)
           
            ## update read_node_dict
            for read in combined_read_ids:
                read_node_dict[read]=combined_node

            G.add_node(combined_node,read_ids=combined_read_ids,penwidth='2.0')
            for predecessors_predecessor in predecessors_predecessors:
                o = G[predecessors_predecessor][predecessor][0]["label"]
                G.add_edge(predecessors_predecessor, combined_node, label= o, color='blue', penwidth='2.0')

            for successor in successors:
                o = G[node_to_combine][successor][0]["label"]
                G.add_edge(combined_node, successor, label= o, color='blue', penwidth='2.0')

            # update sequences
            overlap_to_predecessor = int(G[predecessor][node_to_combine][0]["label"])
            offset = len(read_db[predecessor]) - overlap_to_predecessor

            pred_seq = read_db[predecessor]
            node_seq = read_db[node_to_combine]
            combined_seq = pred_seq + node_seq[overlap_to_predecessor:]

            read_db[combined_node] = combined_seq

            # clean up
            G.remove_node(node_to_combine)
            G.remove_node(predecessor)

            del read_db[node_to_combine]
            del read_db[predecessor]

            if predecessor in nodes_to_combine:
                nodes_to_combine.remove(predecessor)
    return G

def collapse_graph_3(G, read_db, PE_G):
    # node collapsed: combined node
    # PE_G: the pair-end graph
    while True:
        nodes_to_combine = []
        all_node = G.nodes()

        for node in all_node:
            #pdb.set_trace()
            if G.in_degree(node) == 1 and G.out_degree(G.predecessors(node)[0]) == 1:
                predecessor = G.predecessors(node)[0]
                nodes_to_combine.append(node)

        if not nodes_to_combine:
            break

        for node_to_combine in nodes_to_combine:
            if G.in_degree(node_to_combine)==0:
                pdb.set_trace()

            predecessor = G.predecessors(node_to_combine)[0]
            predecessors_predecessors = G.predecessors(predecessor)
            successors = G.successors(node_to_combine)

            # update graph
            ## update species
            pre_species=G.node[predecessor]['species']
            node_to_combine_species=G.node[node_to_combine]['species']
            
            combined_species=join_species(pre_species,node_to_combine_species)

            ## update read id
            combined_read_ids=G.node[predecessor]['read_ids']
            #if pre_read_ids==None:
            #    pdb.set_trace()
            combined_read_ids.extend(G.node[node_to_combine]['read_ids'])

            ## update node
            #combined_node=str(int(combined_read_ids[0].split('/')[0])-1)+'|'+str(len(combined_read_ids)-1)
            if node_to_combine.find('|')>-1 and predecessor.find('|')==-1:
                combined_node = predecessor+'|'+str(int(node_to_combine.split('|')[1].split('-')[0])+1)

            elif node_to_combine.find('|')>-1 and predecessor.find('|')>-1:
                combined_node = predecessor.split('|')[0]+'|'+str(int(predecessor.split('|')[1].split('-')[0])+int(node_to_combine.split('|')[1].split('-')[0])+1)
            elif predecessor.find('|')>-1 and node_to_combine.find('|')==-1:
                combined_node = predecessor.split('|')[0]+'|'+str(int(predecessor.split('|')[1].split('-')[0])+1)
            else:
                combined_node = predecessor +'|'+str(1)
           
            ## update PE_G
            PE_G.add_node(combined_node)
            pair_pres=set(PE_G.predecessors(predecessor))
            pair_pres=pair_pres | set(PE_G.predecessors(node_to_combine))
            for pre in pair_pres:
                o=0
                if predecessor in PE_G[pre]:
                    o+=int(PE_G[pre][predecessor]['label'])
                if node_to_combine in PE_G[pre]:
                    o+=int(PE_G[pre][node_to_combine]['label'])
                PE_G.add_edge(pre,combined_node,label=str(o))

            pair_sucs=set(PE_G.successors(predecessor))
            pair_sucs=pair_sucs | set(PE_G.successors(node_to_combine))
            for suc in pair_sucs:
                o=0
                if suc in PE_G[predecessor]:
                    o+=int(PE_G[predecessor][suc]['label'])
                if suc in PE_G[node_to_combine]:
                    o+=int(PE_G[node_to_combine][suc]['label'])
                PE_G.add_edge(combined_node,suc,label=str(o))

            ## update node color
            if len(combined_species)==1:
                color_node=color_dict[combined_species[0]]
            else:
                color_node='black'
            G.add_node(combined_node,species=combined_species,read_ids=combined_read_ids,color=color_node,penwidth='2.0')
            for predecessors_predecessor in predecessors_predecessors:
                o = G[predecessors_predecessor][predecessor][0]["label"]
                G.add_edge(predecessors_predecessor, combined_node, label= o, color='blue', penwidth='2.0')

            for successor in successors:
                o = G[node_to_combine][successor][0]["label"]
                G.add_edge(combined_node, successor, label= o, color='blue', penwidth='2.0')

            # update sequences
            overlap_to_predecessor = int(G[predecessor][node_to_combine][0]["label"])
            offset = len(read_db[predecessor]) - overlap_to_predecessor

            pred_seq = read_db[predecessor]
            node_seq = read_db[node_to_combine]
            combined_seq = pred_seq + node_seq[overlap_to_predecessor:]

            read_db[combined_node] = combined_seq

            # clean up
            G.remove_node(node_to_combine)
            G.remove_node(predecessor)

            PE_G.remove_node(node_to_combine)
            PE_G.remove_node(predecessor)

            del read_db[node_to_combine]
            del read_db[predecessor]

            if predecessor in nodes_to_combine:
                nodes_to_combine.remove(predecessor)


def is_false_connection(G,node1,node2,paired_end_dict):
    # decide whether (node1,node2) is a false connection in G based on paired-end information
    flag = 1  # 1: is false connection, no paired-end evidence; 0: is not false connection
    predecessors=G.predecessors(node1)
    if (node1,node2) in paired_end_dict or (node2,node1) in paired_end_dict:
        flag=0
    for pre_node in predecessors:
        if (pre_node,node2) in paired_end_dict or (node2, pre_node) in paired_end_dict:
            flag=0
    return flag

def create_paired_end_connection(subgraph_paired_end_edges):
    PE_node_dict={}
    PE_group=[]
    PE_G=nx.Graph()
    idx=0
    for pair in subgraph_paired_end_edges:
        PE_G.add_edge(pair[0],pair[1])
    PE_subgraphs_nodes=nx.connected_components(PE_G)
    for PE_graph_nodes in PE_subgraphs_nodes:
        PE_group.append(PE_graph_nodes)
        for node in PE_graph_nodes:
            if not node in PE_node_dict:
                PE_node_dict[node]=idx
        idx+=1
    return PE_node_dict,PE_group


def create_paired_end_graph(read_node_dict):
    """
    create the pair-end graph from the read_node_dict
    pair_end_edges, key: nodes pair, value: pairs between these two nodes
    """
    paired_end_edges={}
    for read in read_node_dict.keys():
        if int(read)%2==0: # pair.1
            pair=str(int(read)+1)
            if pair in read_node_dict:
                node_1=read_node_dict[read]
                node_2=read_node_dict[pair]
                if not (node_1,node_2) in paired_end_edges:
                    paired_end_edges[(node_1,node_2)]=1
                else:
                    paired_end_edges[(node_1,node_2)]+=1
        elif int(read)%2==1: # pair.2
            pair=str(int(read)-1)
            if pair in read_node_dict:
                node_1=read_node_dict[read]
                node_2=read_node_dict[pair]
                if not (node_2,node_1) in paired_end_edges:
                    paired_end_edges[(node_2,node_1)]=1
                else:
                    paired_end_edges[(node_2,node_1)]+=1
    PE_G=nx.DiGraph()   
    for pair in paired_end_edges:
        PE_G.add_edge(pair[0],pair[1],color='red',label=str(paired_end_edges[pair]))
    return paired_end_edges,PE_G

def create_paired_end_graph_with_pairs(read_node_dict, pairs_dict):
    """
    create the pair-end graph from read_node_dict and pairs_dict
    pairs_dict, key: nodes pair, in the index way
    pair_end_edges, key: nodes pair, value: pairs between these two nodes
    """
    paired_end_edges={}
    for pair in pairs_dict.keys():
        if (not pair[0] in read_node_dict) or (not pair[1] in read_node_dict):
            continue
        pair_n1 = read_node_dict[pair[0]]
        pair_n2 = read_node_dict[pair[1]]
        if not (pair_n1, pair_n2) in paired_end_edges:
            paired_end_edges[(pair_n1, pair_n2)]=1
        else:
            paired_end_edges[(pair_n1, pair_n2)]+=1
    
    #PE_G=nx.DiGraph()
    #for pair in paired_end_edges:
    #    PE_G.add_edge(pair[0], pair[1], color='red', label=str(paired_end_edges[pair]))
    
    PE_G=nx.Graph()
    for pair in paired_end_edges:
        if PE_G.has_edge(pair[0], pair[1]):
            new_label = str(paired_end_edges[pair]+int(PE_G[pair[0]][pair[1]]['label']))
            PE_G.add_edge(pair[0], pair[1], color='red', label = new_label)
        else:
            PE_G.add_edge(pair[0], pair[1], color='red', label = str(paired_end_edges[pair]))

    return paired_end_edges, PE_G

def pair_end_binning(G, PE_G, overlap_cutoff):
    """
    Algorithm: for each edge in overlap graph (node1--node2), if no pair-end connection between the two edge nodes, no pair-end connection between predecessors of node1 and node2, in_degree of node2 >=2, delete the edge in overlap graph
    """

    for E in G.edges():
        if (not E[1] in PE_G[E[0]]) and G.in_degree(E[1])>1 and G.out_degree(E[0])>1:
            flag1,flag2=0,0
            for pre in G.predecessors(E[0]):
                if G.out_degree(pre)==1 and (E[1] in PE_G[pre]):
                    flag1+=1
            for suc in G.successors(E[1]):
                if G.in_degree(suc)==1 and suc in PE_G[E[0]]:
                    flag2+=1
            overlap=int(G[E[0]][E[1]][0]['label'])
            if not flag1 and not flag2:
                if overlap<overlap_cutoff:
                    G.remove_edge(E[0],E[1])

    for N in G.nodes():
        if G.in_degree(N)==1 and G.out_degree(N)==1:
            continue
        if G.in_degree(N)==1 and G.out_degree(N)>1:
            flag=0
            for suc in G.successors(N):
                if suc in PE_G[N]:
                    flag+=1
            if flag>0:
                for suc in G.successors(N):
                    if not suc in PE_G[N] and G.in_degree(suc)>1:
                        G.remove_edge(N,suc)
            else:
                flag=0
                for suc in G.successors(N):
                    if G.in_degree(suc)==1 and G.out_degree(suc)==1:
                        flag+=1
                if flag==1:
                    for suc in G.successors(N):
                        if G.in_degree(suc)>1:
                            G.remove_edge(N,suc)

        elif G.in_degree(N)>1 and G.out_degree(N)==1:
            flag=0
            for pre in G.predecessors(N):
                if N in PE_G[pre]:
                    flag+=1
            if flag>0:
                for pre in G.predecessors(N):
                    if not N in PE_G[pre] and G.out_degree(pre)>1:
                        G.remove_edge(pre,N)
            else:
                flag=0
                for pre in G.predecessors(N):
                    if G.out_degree(pre)==1 and G.in_degree(pre)==1:
                        flag+=1
                if flag==1:
                    for pre in G.predecessors(N):
                        if G.out_degree(pre)>1:
                            G.remove_edge(pre,N)


def get_part_graph(G, start_node, nodes_num):
    sub_graph_nodes=set([start_node])
    stack=[[start_node]]
    loop=0
    while (len(sub_graph_nodes)<=nodes_num) and stack:
        succ_nodes=stack.pop()
        next_succ_nodes=set([])
        for succ_node in succ_nodes:
            sub_graph_nodes=sub_graph_nodes.union(set(G.successors(succ_node)))
            sub_graph_nodes=sub_graph_nodes.union(set(G.predecessors(succ_node)))
            next_succ_nodes=next_succ_nodes.union(set(G.successors(succ_node)))
        stack.append(list(next_succ_nodes))
        loop+=1
        if loop>10*nodes_num:
            break
    return G.subgraph(sub_graph_nodes)

def plot_graph(G, figname):
    G_plot=nx.drawing.nx_agraph.to_agraph(G)
    G_plot.draw(figname,prog='dot')

def remove_tips(G, seq_db, node_len_cutoff):
    """
    Algorithm: for nodes with either multiple predecessors or multiple successors:
    The removed predecessors or successors are either starting nodes or ending nodes
    (1) if at least one predecessor or successor has more connections rather than the current node, remove other predecessors or successors that only have connections to the current node
    (2) if all predecessors or successors only have connections to the current node, keep the longest predecessor or successor 
    """

    branch_start_nodes=[n for n in G.nodes() if len(G.predecessors(n))>1]
    branch_end_nodes=[n for n in G.nodes() if len(G.successors(n))>1]
    
    for branch_node in branch_start_nodes:
        flag=0
        max_len=0
        max_len_node=''
        for pred_node in G.predecessors(branch_node):
            if max_len<len(seq_db[pred_node]):
                max_len=len(seq_db[pred_node])
                max_len_node=pred_node
            if G.in_degree(pred_node)>0 or G.out_degree(pred_node)>1:
                flag=1
                break
        if flag:
            for pred_node in G.predecessors(branch_node):
                if G.in_degree(pred_node)==0 and G.out_degree(pred_node)==1 and len(seq_db[pred_node])<node_len_cutoff:
                    G.remove_node(pred_node)
        else:
            for pred_node in G.predecessors(branch_node):
                if G.in_degree(pred_node)==0 and G.out_degree(pred_node)==1 and len(seq_db[pred_node])<max_len:
                    G.remove_node(pred_node)

    for branch_node in branch_end_nodes:
        flag=0
        max_len=0
        max_len_node=''
        for succ_node in G.successors(branch_node):
            if max_len<len(seq_db[succ_node]):
                max_len=len(seq_db[succ_node])
                max_len_node=succ_node
            if G.in_degree(succ_node)>1 or G.out_degree(succ_node)>0:
                flag=1
                break
        if flag:
            for succ_node in G.successors(branch_node):
                if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==0 and len(seq_db[succ_node])<node_len_cutoff:
                    G.remove_node(succ_node)
        else:
            for succ_node in G.successors(branch_node):
                if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==0 and len(seq_db[succ_node])<max_len:
                    G.remove_node(succ_node)

def remove_super_tips(G, seq_db):
    """
    Algorithm: for each staring node, if all of its successors(at least two) have more than one predecessor, remove the starting node
    for each ending node, if all of its predecessors(at least two) have more than one successor, remove the ending node
    """

    starting_nodes=[n for n in G.nodes() if G.in_degree(n)==0]
    ending_nodes=[n for n in G.nodes() if G.out_degree(n)==0]
    for start_node in starting_nodes:
        flag=0
        if G.out_degree(start_node)>1:
            for succ_node in G.successors(start_node):
                if G.in_degree(succ_node)==1:
                    flag=1
                    break
            if not flag:
                G.remove_node(start_node)
    
    for end_node in ending_nodes:
        flag=0
        if G.in_degree(end_node)>1:
            for pred_node in G.predecessors(end_node):
                if G.out_degree(pred_node)==1:
                    flag=1
                    break
            if not flag:
                G.remove_node(end_node)


