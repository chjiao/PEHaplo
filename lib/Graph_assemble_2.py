import re,sys,pdb
import networkx as nx

def SISO_clustering(G,PE_G):
    SISO_list=[]
    for N in G:
        if G.in_degree(N)<=1 and G.out_degree(N)<=1:
            SISO_list.append(N)
    
    cluster_group=[]
    SISO_dict={}
    for N in SISO_list:
        #pdb.set_trace()
        SISO_dict[N]=1
        stack=[N]
        group=set([N])
        while stack:
            #pdb.set_trace()
            SISO_N=stack.pop()
            for pre in PE_G.predecessors(SISO_N):
                group.add(pre)
                if G.in_degree(pre)<=1 and G.out_degree(pre)<=1 and not pre in SISO_dict:
                    stack.append(pre)
                    #if G.node[pre]['species']!=G.node[SISO_N]['species']:
                    #    pdb.set_trace()
                    SISO_dict[pre]=1
                    if pre in SISO_list:
                        SISO_list.remove(pre)
            for suc in PE_G.predecessors(SISO_N):
                group.add(suc)
                if G.in_degree(suc)<=1 and G.out_degree(suc)<=1 and not suc in SISO_dict:
                    stack.append(suc)
                    #if G.node[suc]['species']!=G.node[SISO_N]['species']:
                    #    pdb.set_trace()
                    SISO_dict[suc]=1
                    if suc in SISO_list:
                        SISO_list.remove(suc)
        cluster_group.append(group)
    return cluster_group

def SISO_clustering_2(G,PE_G):
    # return a graph containing the clustering results
    SISO_list=[]
    for N in G:
        if G.in_degree(N)<=1 and G.out_degree(N)<=1:
            SISO_list.append(N)
    
    #cluster_group=[]
    SISO_dict={}
    PE_group_G=nx.Graph()
    for N in SISO_list:
        #pdb.set_trace()
        SISO_dict[N]=1
        stack=[N]
        group=set([N])
        while stack:
            #pdb.set_trace()
            SISO_N=stack.pop()
            for pre in PE_G.predecessors(SISO_N):
                group.add(pre)
                if G.in_degree(pre)<=1 and G.out_degree(pre)<=1 and not pre in SISO_dict:
                    stack.append(pre)
                    SISO_dict[pre]=1
                    if pre in SISO_list:
                        SISO_list.remove(pre)
            for suc in PE_G.predecessors(SISO_N):
                group.add(suc)
                if G.in_degree(suc)<=1 and G.out_degree(suc)<=1 and not suc in SISO_dict:
                    stack.append(suc)
                    SISO_dict[suc]=1
                    if suc in SISO_list:
                        SISO_list.remove(suc)
        group=list(group)
        for i in range(len(group)-1):
            for j in range(i+1,len(group)):
                PE_group_G.add_edge(group[i],group[j])

    for N in G.nodes():
        if not N in PE_group_G:
            PE_group_G.add_node(N)
    
    return PE_group_G


def DFS_paths_interative(G, start_node, end_node):
    stack=[(start_node,[start_node])]
    while stack:
        (vertex, path)=stack.pop()
        for succ_node in set(G.successors(vertex)) - set(path):
            #pdb.set_trace()
            if succ_node==end_node:
                yield path+[succ_node]
            else:
                stack.append((succ_node,path+[succ_node]))

def get_paired_score(path,succ_node,paired_end_edges):
    paired_score1,paired_score2=0,0
    for this_node in path[0:-1]:
        if (this_node,succ_node) in paired_end_edges: 
            paired_score1+=paired_end_edges[(this_node,succ_node)]
        elif (succ_node,this_node) in paired_end_edges:
            paired_score1+=paired_end_edges[(succ_node,this_node)]
    if (path[-1],succ_node) in paired_end_edges:
        paired_score2+=paired_end_edges[(path[-1],succ_node)]
    elif (succ_node,path[-1]) in paired_end_edges:
        paired_score2+=paired_end_edges[(succ_node,path[-1])]

    return paired_score1,paired_score2

def get_paired_score_2(path,succ_node,PE_G):
    score=0
    for N in path:
        if N in PE_G:
            if succ_node in PE_G[N]:
                o=int(PE_G[N][succ_node]['label'])
                score+=o
    return score

def get_paired_score_plan(path,plan_nodes,PE_G):
    score=0
    for N1 in plan_nodes:
        for N in path:
            if N in PE_G:
                if N1 in PE_G[N]:
                    o=int(PE_G[N][N1]['label'])
                    score+=o
    return score

def get_paired_connection_score(path,succ_node,PE_node_dict):
    score1,score2=0,0
    if not succ_node in PE_node_dict:
        return 0,0
    else:
        for this_node in path[0:-1]:
            if this_node in PE_node_dict and PE_node_dict[this_node]==PE_node_dict[succ_node]:
                score1+=1
        if path[-1] in PE_node_dict and PE_node_dict[path[-1]]==PE_node_dict[succ_node]:
            score2+=1
        return score1,score2

# using paired-end information for paths searching
# rule 1: if the new node has connection to the path except its direct predecessor (denote as pre-path), add the new node to the path
# rule 2: if the new node has no connection to the pre-path but the direct predecessor, do not add the new node to the path
def DFS_paths_paired_end(G, start_node, end_node,paired_end_edges,PE_node_dict):
    # PE_node_dict: the index of group for each node, key: node, value: index, group: pair-end connected groups
    stack=[(start_node,[start_node])]
    while stack:
        (vertex, path)=stack.pop()
        succ_node_all = set(G.successors(vertex)) - set(path)
        if len(succ_node_all)>1: # seeing a bifurcation node, multiple successors, judge which one to append
            flag=0
            for succ_node in succ_node_all:
                score1,score2=get_paired_score(path,succ_node,paired_end_edges)
                connect_score1,connect_score2=get_paired_connection_score(path,succ_node,PE_node_dict)
                if succ_node==end_node: # reaching the end
                    if score1>0:
                        yield path+[succ_node]
                    elif connect_score1>0:
                        yield path+[succ_node]
                    elif score2>0:
                        #yield path[-1:]+[succ_node]
                        print "1. Path:",path
                        print "succ_node:",succ_node
                    elif connect_score2>0:
                        print "2. Path:",path
                        print "succ_node:",succ_node
                    else:
                        flag+=1
                        print "No connection succint node (1):",path,succ_node
                else:                    # not reaching the end
                    if score1>0:
                        stack.append((succ_node,path+[succ_node]))
                    elif connect_score1>0:
                        stack.append((succ_node,path+[succ_node]))
                    elif score2>0:
                        #pdb.set_trace()
                        #stack.append((succ_node,path[-1:]+[succ_node]))
                        print "1_middle. Path:",path
                        print "succ_node:",succ_node
                    elif connect_score2>0:
                        print "2_middle. Path:",path
                        print "succ_node:",succ_node
                    else:
                        flag+=1
                        print "No connection succint node (2_middle):",path,succ_node
                if flag==len(succ_node_all): # cannot extend to any of the next node
                    #yield path
                    print "Break middle node found!",vertex
                    #pdb.set_trace()

        else: # just one successor, append to the path
            for succ_node in succ_node_all:
                if succ_node==end_node:
                    yield path+[succ_node]
                else:
                    stack.append((succ_node,path+[succ_node]))

def get_max_score_node(score_dict):
    max_score=0
    max_score_node=0
    for key in score_dict:
        if score_dict[key]>max_score:
            max_score=score_dict[key]
            max_score_node=key
    return max_score_node

def classify_SISO_score_node(SISO_dict,plan_dict,PE_plan_dict,path_dict,path_plan_dict):
    plan_flag=sum(i>0 for i in plan_dict.values())
    PE_plan_flag=sum(i>0 for i in PE_plan_dict.values())
    path_flag=sum(i>0 for i in path_dict.values())
    path_plan_flag=sum(i>0 for i in path_plan_dict.values())
    score_vec=SISO_dict.values()
    diff=[score_vec[i+1]-score_vec[i] for i in range(len(score_vec)-1)]
    
    max_score_node=0
    if sum(diff)==0: # same SISO_score for multiple successors
        if plan_flag>0:
            max_score_node=get_max_score_node(plan_dict)
        elif PE_plan_flag>0:
            max_score_node=get_max_score_node(PE_plan_dict)
        elif path_flag>0:
            max_score_node=get_max_score_node(path_dict)
        elif path_plan_flag>0:
            max_score_node=get_max_score_node(path_plan_dict)
        else:
            #pdb.set_trace()
            max_score_node=get_max_score_node(SISO_dict)
    else: # choose the node with maximum SISO score
        max_score_node=get_max_score_node(SISO_dict)
    return max_score_node


def bifurcation_classifier(N,G,PE_G,path,SISO,read_db,fragment_len):
    max_score_node=0
    if len(read_db[N])<fragment_len:
        SISO_flag=0
        SISO_score_dict={}
        path_flag=0
        path_score_dict={}
        for succ_node in G.successors(N):
            plan_nodes=[succ_node]
            while len(plan_nodes)<4 and G.out_degree(plan_nodes[-1])>0:
                plan_nodes.extend(G.successors(plan_nodes[-1]))
            SISO_score=get_paired_score_plan(SISO,plan_nodes,PE_G)
            path_score=get_paired_score_plan(path,plan_nodes,PE_G)
            if SISO_score>0:
                SISO_flag+=1
                SISO_score_dict[succ_node]=SISO_score
            if path_score>0:
                path_flag+=1
                path_score_dict[succ_node]=path_score
        if SISO_flag==1:
            max_score_node=get_max_score_node(SISO_score_dict)
        elif path_flag==1:
            max_score_node=get_max_score_node(path_score_dict)
        elif SISO_flag>1:
            max_score_node=get_max_score_node(SISO_score_dict)
        elif path_flag>1:
            max_score_node=get_max_score_node(path_score_dict)
        #else:
        #    pdb.set_trace()
        return max_score_node
    else:
        # classify by coverage, etc.
        # calculate the coverage of the SISO nodes, compare with bifurcations
        node_read_num=0
        node_len=0
        for this_node in SISO:
            node_read_num+=len(G.node[this_node]['read_ids'])
            node_len+=len(read_db[this_node])
        if node_len>0:
            cov_SISO=float(node_read_num)/node_len
        else:
            cov_SISO=0

        cov_diff_list=[]
        min_diff=10000
        min_diff_node=0
        for succ_node in G.successors(N):
            cov_succ_node=len(G.node[succ_node]['read_ids'])/float(len(read_db[succ_node]))
            diff=abs(cov_succ_node-cov_SISO)
            if diff<min_diff:
                min_diff=diff
                min_diff_node=succ_node
        # calculate the pair-end connections
        con_PE_num=0
        con_num=0
        min_diff_con_node=0
        #for pre in PE_G.predecessors(N):
        for pre in PE_G[N]:
            if pre in SISO:
                con_PE_num+=int(PE_G[pre][N]['label'])
                con_num+=1
        if con_num>0:
            con_SISO=float(con_PE_num)/con_num

            min_diff=10000
            for succ_node in G.successors(N):
                if succ_node in PE_G[N]:
                    con_succ_node=int(PE_G[N][succ_node]['label'])
                    diff=abs(con_SISO-con_succ_node)
                    if diff<min_diff:
                        min_diff=diff
                        min_diff_con_node=succ_node
            if min_diff_con_node!=0:
                min_diff_node=min_diff_con_node
        return min_diff_node

def DFS_paths_paired_end_2(G, start_node, PE_G):
    # PE_node_dict: the index of group for each node, key: node, value: index, group: pair-end connected groups
    stack=[(start_node,[start_node],[start_node])]
    while stack:
        (vertex,path,SISO)=stack.pop()
        succ_node_all = set(G.successors(vertex)) - set(path)
        if len(succ_node_all)<G.out_degree(vertex):
            print "Finding a cycle:",path
            yield path
            continue
        if len(succ_node_all)>1: # seeing a bifurcation node, multiple successors, judge which one to append
            flag=0
            plan_flag=0
            path_plan_flag=0
            path_flag=0
            path_score_dict={}
            plan_score_dict={}
            path_plan_score_dict={}
            for succ_node in succ_node_all:
                plan_nodes=[succ_node]
                while G.out_degree(plan_nodes[-1])==1 and G.out_degree(G.successors(plan_nodes[-1])[0])==1 and (not G.successors(plan_nodes[-1])[0] in plan_nodes): 
                    plan_nodes.append(G.successors(plan_nodes[-1])[0])
                #plan_nodes.extend(G.successors(plan_nodes[-1]))
                SISO_plan_score=get_paired_score_plan(SISO,plan_nodes,PE_G)
                path_plan_score=get_paired_score_plan(path[:-1],plan_nodes,PE_G)
                path_score=get_paired_score_2(path,succ_node,PE_G)
                if SISO_plan_score>0:
                    plan_flag+=1
                if path_plan_score>0:
                    path_plan_flag+=1
                if path_score>0:
                    path_flag+=1
                plan_score_dict[succ_node]=SISO_plan_score
                path_plan_score_dict[succ_node]=path_plan_score
                path_score_dict[succ_node]=path_score
            #pdb.set_trace()
            if plan_flag==0:
                if path_flag==1:
                    max_score_suc=get_max_score_node(path_score_dict)
                elif path_plan_flag==1:
                    max_score_suc=get_max_score_node(path_plan_score_dict)
                else:
                    #pdb.set_trace()
                    max_score_suc=bifurcation_classifier(vertex,G,PE_G,path,SISO,read_db,fragment_len)

            for succ_node in succ_node_all:
                #if succ_node=='[YU2][43.8][274]5893|59':
                #    pdb.set_trace()
                SISO_score=get_paired_score_2(SISO,succ_node,PE_G)
                if G.out_degree(succ_node)==0: # reaching the end
                    if SISO_score>0:
                        yield path+[succ_node]
                    elif len(path)>1:
                        if succ_node in PE_G[path[-2]]:
                            yield path+[succ_node]
                    elif len(path)==1:
                        if succ_node in PE_G[path[-1]]:
                            yield path+[succ_node]
                    else:
                        flag+=1
                else:                    # not reaching the end
                    if SISO_score>0:
                        if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                            if not succ_node in SISO:
                                SISO.append(succ_node)
                        stack.append((succ_node,path+[succ_node],SISO))
                    else:
                        if plan_flag>0:
                            if plan_score_dict[succ_node]>0:
                                if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                    if not succ_node in SISO:
                                        SISO.append(succ_node)
                                stack.append((succ_node,path+[succ_node],SISO))
                        else:
                            if path_flag==1 and succ_node==max_score_suc:
                                if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                    if not succ_node in SISO:
                                        SISO.append(succ_node)
                                stack.append((succ_node,path+[succ_node],SISO))
                            elif path_plan_flag==1 and succ_node==max_score_suc:
                                if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                    if not succ_node in SISO:
                                        SISO.append(succ_node)
                                stack.append((succ_node,path+[succ_node],SISO))
                            elif max_score_suc!=0 and succ_node==max_score_suc:
                                if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                    if not succ_node in SISO:
                                        SISO.append(succ_node)
                                stack.append((succ_node,path+[succ_node],SISO))
                            else:
                                flag+=1
            if flag==len(succ_node_all):
                print "Break middle node found!",path,vertex

        else: # just one successor, append to the path
            for succ_node in succ_node_all:
                if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                    if not succ_node in SISO:
                        SISO.append(succ_node)
                if G.out_degree(succ_node)==0:
                    yield path+[succ_node]
                else:
                    stack.append((succ_node,path+[succ_node],SISO))

def DFS_paths_paired_end_3(G, start_node, PE_G, read_db, fragment_len=300):
    # PE_node_dict: the index of group for each node, key: node, value: index, group: pair-end connected groups
    #if start_node=='53805|4':
    #if start_node=='41868|2':
    if start_node=='42075|8':
        pdb.set_trace()
    stack=[(start_node,[start_node],[start_node])]
    backtrack_flag=0
    backtrack_dict={}
    while stack:
        (vertex,path,SISO)=stack.pop()
        if vertex==0:
            pdb.set_trace()
        succ_node_all = set(G.successors(vertex)) - set(path)

        flag=0 #indicate whether the node is appended
        SISO_flag,plan_flag,PE_plan_flag,path_flag,path_plan_flag=0,0,0,0,0
        SISO_score_dict,plan_score_dict,PE_plan_score_dict,path_score_dict,path_plan_score_dict={},{},{},{},{}
        for succ_node in succ_node_all:
            plan_nodes=[succ_node]
            while G.out_degree(plan_nodes[-1])==1 and G.out_degree(G.successors(plan_nodes[-1])[0])==1 and (not G.successors(plan_nodes[-1])[0] in plan_nodes): 
                plan_nodes.append(G.successors(plan_nodes[-1])[0])
            #if len(PE_G.predecessors(succ_node))==0 and len(plan_nodes)==1: #if succ_node has no predecessors in PE_G graph, consider its successors
            if len(plan_nodes)==1: #if cannot plan to succ_node's successors, plan_flag will be same as SISO_flag, consider its successors directly
                plan_nodes=G.successors(succ_node)
            ##
            plan_nodes=list(set(plan_nodes)-set(path))
            #plan_nodes.extend(G.successors(plan_nodes[-1]))
            SISO_score=get_paired_score_2(SISO,succ_node,PE_G)
            SISO_plan_score=get_paired_score_plan(SISO,plan_nodes,PE_G)
            PE_plan_nodes=list(set(PE_G.successors(succ_node))-set(path))
            PE_plan_score=get_paired_score_plan(SISO,PE_plan_nodes,PE_G)
            path_score=get_paired_score_2(path[:-1],succ_node,PE_G)
            path_plan_score=get_paired_score_plan(path[:-1],plan_nodes,PE_G)
            if SISO_score>0:
                SISO_flag+=1
            if SISO_plan_score>0:
                plan_flag+=1
            if PE_plan_score>0:
                PE_plan_flag+=1
            if path_score>0:
                path_flag+=1
            if path_plan_score>0:
                path_plan_flag+=1
            SISO_score_dict[succ_node]=SISO_score
            plan_score_dict[succ_node]=SISO_plan_score
            PE_plan_score_dict[succ_node]=PE_plan_score
            path_score_dict[succ_node]=path_score
            path_plan_score_dict[succ_node]=path_plan_score
        #pdb.set_trace()
        max_score_suc=0
        if SISO_flag>=1:
            #max_score_suc=get_max_score_node(SISO_score_dict)
            max_score_suc=classify_SISO_score_node(SISO_score_dict,plan_score_dict,PE_plan_score_dict,path_score_dict,path_plan_score_dict)
        elif plan_flag>=1:
            max_score_suc=get_max_score_node(plan_score_dict)
        elif PE_plan_flag>=1:
            max_score_suc=get_max_score_node(PE_plan_score_dict)
        elif path_flag>=1:
            max_score_suc=get_max_score_node(path_score_dict)
        elif path_plan_flag>=1:
            max_score_suc=get_max_score_node(path_plan_score_dict)
        else:  ## all the flags are zero
            #pdb.set_trace()
            max_score_suc=bifurcation_classifier(vertex,G,PE_G,path,SISO,read_db,fragment_len)
      
        if backtrack_flag==1:
            path=path[:-1]
        else:
            if len(succ_node_all)<G.out_degree(vertex):
                #pdb.set_trace()
                print "Finding a cycle:",path
                #yield path
                #continue
                #'''
                if G.out_degree(vertex)==1 and path.count(vertex)==1:
                    stack.append((G.successors(vertex)[0],path+[G.successors(vertex)[0]],SISO))
                    continue
                elif SISO_flag>=1 or plan_flag>=1:
                    for succ_node in succ_node_all:
                        if succ_node==max_score_suc:
                            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                SISO.append(succ_node)
                            stack.append((succ_node,path+[succ_node],SISO))
                            if G.out_degree(succ_node)==0:
                                yield path
                                break
                    continue
                else:
                    yield path
                    continue
                #'''

        if len(succ_node_all)>1: # seeing a bifurcation node, multiple successors, judge which one to append
            for succ_node in succ_node_all:
                if G.out_degree(succ_node)==0: # reaching the end, bug fixed
                    if SISO_score>0:
                        yield path+[succ_node]
                    elif len(path)>1 and succ_node in PE_G[path[-2]]:
                        yield path+[succ_node]
                    elif len(path)==1 and succ_node in PE_G[path[-1]]:
                        yield path+[succ_node]
                    else:
                        flag+=1
                        ## remove the flag from the ending node
                        if SISO_flag>0 and SISO_score_dict[succ_node]>0:
                            SISO_score_dict[succ_node]=0
                            SISO_flag-=1
                        if plan_flag>0 and plan_score_dict[succ_node]>0:
                            plan_score_dict[succ_node]=0
                            plan_flag-=1
                        if PE_plan_flag>0 and PE_plan_score_dict[succ_node]>0:
                            PE_plan_score_dict[succ_node]=0
                            PE_plan_flag-=1
                        if path_flag>0 and path_score_dict[succ_node]>0:
                            path_score_dict[succ_node]=0
                            path_flag-=1
                        if path_plan_flag>0 and path_plan_score_dict[succ_node]>0:
                            path_plan_score_dict[succ_node]=0
                            path_plan_flag-=1
                else:                    # not reaching the end
                    if SISO_flag==1:
                        if SISO_score_dict[succ_node]>0:
                            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                SISO.append(succ_node)
                            stack.append((succ_node,path+[succ_node],SISO))
                    elif SISO_flag>1:
                        if succ_node==max_score_suc:
                            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                SISO.append(succ_node)
                            stack.append((succ_node,path+[succ_node],SISO))
                    else:
                        if plan_flag==1:
                            if plan_score_dict[succ_node]>0:
                                if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                    SISO.append(succ_node)
                                stack.append((succ_node,path+[succ_node],SISO))
                        elif plan_flag>1:
                            if succ_node==max_score_suc:
                                if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                    SISO.append(succ_node)
                                stack.append((succ_node,path+[succ_node],SISO))
                        else:
                            if PE_plan_flag==1:
                                if PE_plan_score_dict[succ_node]>0:
                                    if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                        SISO.append(succ_node)
                                    stack.append((succ_node,path+[succ_node],SISO))
                            elif PE_plan_flag>1:
                                if succ_node==max_score_suc:
                                    if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                        SISO.append(succ_node)
                                    stack.append((succ_node,path+[succ_node],SISO))
                            else:
                                if path_flag==1:
                                    if path_score_dict[succ_node]>0:
                                        if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                            SISO.append(succ_node)
                                        stack.append((succ_node,path+[succ_node],SISO))
                                elif path_flag>1:
                                    if succ_node==max_score_suc:
                                        if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                            SISO.append(succ_node)
                                        stack.append((succ_node,path+[succ_node],SISO))
                                else:
                                    if path_plan_flag==1:
                                        if path_plan_score_dict[succ_node]>0:   ## bug fixed
                                            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                                SISO.append(succ_node)
                                            stack.append((succ_node,path+[succ_node],SISO))
                                    elif path_plan_flag>1:
                                        if succ_node==max_score_suc:
                                            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                                SISO.append(succ_node)
                                            stack.append((succ_node,path+[succ_node],SISO))
                                    else:
                                        if max_score_suc!=0 and succ_node==max_score_suc:
                                            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                                SISO.append(succ_node)
                                            stack.append((succ_node,path+[succ_node],SISO))
                                        else:
                                            #pdb.set_trace()
                                            flag+=1
            if flag==len(succ_node_all):
                #pdb.set_trace()
                print "Break middle node found!",path,vertex
                if len(path)>1:
                    ## backtrack
                    #if backtrack_flag==0 and G.successors(path[-2])>1:
                    if not path[-2] in backtrack_dict:
                        backtrack_dict[path[-2]]=[]
                    if backtrack_flag==0 and len(set(G.successors(path[-2]))-set(backtrack_dict[path[-2]])-set(path[-1]))>0:
                        stack.append((path[-2],path,SISO))
                        backtrack_flag=1
                        if not path[-1] in backtrack_dict[path[-2]]:
                            backtrack_dict[path[-2]].append(path[-1])
                    else: # cannot backtrack or backtrack failed, yield the path
                        backtrack_flag=0
                        yield path
                else:
                    yield path
            elif backtrack_flag>0: #successfully backtracked
                backtrack_flag=0

        elif len(succ_node_all)==1: # just one successor, append to the path
            if backtrack_flag==0:
                for succ_node in succ_node_all:
                    if not SISO_flag and not plan_flag and not PE_plan_flag and not path_flag and not path_plan_flag and len(path)>1 and len(read_db[vertex])<fragment_len:

                        ## backtrack
                        #if backtrack_flag==0 and G.successors(path[-2])>1:
                        if not path[-2] in backtrack_dict:
                            backtrack_dict[path[-2]]=[]
                        if backtrack_flag==0 and len(set(G.successors(path[-2]))-set(backtrack_dict[path[-2]])-set(path[-1]))>0:
                            stack.append((path[-2],path,SISO))
                            backtrack_flag=1
                            if not path[-1] in backtrack_dict[path[-2]]:
                                backtrack_dict[path[-2]].append(path[-1])
                        else: # backtrack failed, yield the path
                            backtrack_flag=0
                            yield path
                    elif SISO_flag or plan_flag:
                        if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                            SISO.append(succ_node)
                        if G.out_degree(succ_node)==0:
                            yield path+[succ_node]
                        else:
                            stack.append((succ_node,path+[succ_node],SISO))
                    else:
                        yield path
            else:
                if SISO_flag or plan_flag or PE_plan_flag or path_flag or path_plan_flag:
                    backtrack_flag=0
                    if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                        SISO.append(succ_node)
                    if G.out_degree(succ_node)==0:
                        yield path+[succ_node]
                    else:
                        stack.append((succ_node,path+[succ_node],SISO))
                else:
                    yield path
        else: # isolated node
            yield path

def DFS_paths_paired_end_4(G, start_node, PE_G, read_db, fragment_len=300):
    # PE_node_dict: the index of group for each node, key: node, value: index, group: pair-end connected groups
    stack=[(start_node,[start_node],[start_node])]
    while stack:
        (vertex,path,SISO)=stack.pop()
        succ_node_all = set(G.successors(vertex)) - set(path)

        flag=0 #indicate whether the node is appended
        SISO_flag,plan_flag,PE_plan_flag,path_flag,path_plan_flag=0,0,0,0,0
        SISO_score_dict,plan_score_dict,PE_plan_score_dict,path_score_dict,path_plan_score_dict={},{},{},{},{}
        for succ_node in succ_node_all:
            plan_nodes=[succ_node]
            while G.out_degree(plan_nodes[-1])==1 and G.out_degree(G.successors(plan_nodes[-1])[0])==1 and (not G.successors(plan_nodes[-1])[0] in plan_nodes): 
                plan_nodes.append(G.successors(plan_nodes[-1])[0])
            #if len(PE_G.predecessors(succ_node))==0 and len(plan_nodes)==1: #if succ_node has no predecessors in PE_G graph, consider its successors
            if len(plan_nodes)==1: #if cannot plan to succ_node's successors, plan_flag will be same as SISO_flag, consider its successors directly
                plan_nodes=G.successors(succ_node)
            ##
            plan_nodes=list(set(plan_nodes)-set(path))
            SISO_score=get_paired_score_2(SISO,succ_node,PE_G)
            SISO_plan_score=get_paired_score_plan(SISO,plan_nodes,PE_G)
            PE_plan_nodes=list(set(PE_G.successors(succ_node))-set(path))
            PE_plan_score=get_paired_score_plan(SISO,PE_plan_nodes,PE_G)
            path_score=get_paired_score_2(path[:-1],succ_node,PE_G)
            path_plan_score=get_paired_score_plan(path[:-1],plan_nodes,PE_G)
            if SISO_score>0:
                SISO_flag+=1
            if SISO_plan_score>0:
                plan_flag+=1
            if PE_plan_score>0:
                PE_plan_flag+=1
            if path_score>0:
                path_flag+=1
            if path_plan_score>0:
                path_plan_flag+=1
            SISO_score_dict[succ_node]=SISO_score
            plan_score_dict[succ_node]=SISO_plan_score
            PE_plan_score_dict[succ_node]=PE_plan_score
            path_score_dict[succ_node]=path_score
            path_plan_score_dict[succ_node]=path_plan_score
        #pdb.set_trace()
        max_score_suc=0
        if SISO_flag>=1:
            #max_score_suc=get_max_score_node(SISO_score_dict)
            max_score_suc=classify_SISO_score_node(SISO_score_dict,plan_score_dict,PE_plan_score_dict,path_score_dict,path_plan_score_dict)
        elif plan_flag>=1:
            max_score_suc=get_max_score_node(plan_score_dict)
        elif PE_plan_flag>=1:
            max_score_suc=get_max_score_node(PE_plan_score_dict)
        elif path_flag>=1:
            max_score_suc=get_max_score_node(path_score_dict)
        elif path_plan_flag>=1:
            max_score_suc=get_max_score_node(path_plan_score_dict)
        else:  ## all the flags are zero
            #pdb.set_trace()
            max_score_suc=bifurcation_classifier(vertex,G,PE_G,path,SISO,read_db,fragment_len)
      
        if len(succ_node_all)<G.out_degree(vertex):
            print "Finding a cycle:",path
            yield path
            continue

        if len(succ_node_all)>1: # seeing a bifurcation node, multiple successors, judge which one to append
            for succ_node in succ_node_all:
                if G.out_degree(succ_node)==0: # reaching the end, bug fixed
                    if SISO_score>0:
                        yield path+[succ_node]
                    elif len(path)>1 and succ_node in PE_G[path[-2]]:
                        yield path+[succ_node]
                    elif len(path)==1 and succ_node in PE_G[path[-1]]:
                        yield path+[succ_node]
                    else:
                        flag+=1
                        ## remove the flag from the ending node
                        if SISO_flag>0 and SISO_score_dict[succ_node]>0:
                            SISO_score_dict[succ_node]=0
                            SISO_flag-=1
                        if plan_flag>0 and plan_score_dict[succ_node]>0:
                            plan_score_dict[succ_node]=0
                            plan_flag-=1
                        if PE_plan_flag>0 and PE_plan_score_dict[succ_node]>0:
                            PE_plan_score_dict[succ_node]=0
                            PE_plan_flag-=1
                        if path_flag>0 and path_score_dict[succ_node]>0:
                            path_score_dict[succ_node]=0
                            path_flag-=1
                        if path_plan_flag>0 and path_plan_score_dict[succ_node]>0:
                            path_plan_score_dict[succ_node]=0
                            path_plan_flag-=1
                else:                    # not reaching the end
                    if SISO_flag==1:
                        if SISO_score_dict[succ_node]>0:
                            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                if not succ_node in SISO:
                                    SISO.append(succ_node)
                            stack.append((succ_node,path+[succ_node],SISO))
                    elif SISO_flag>1:
                        if succ_node==max_score_suc:
                            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                if not succ_node in SISO:
                                    SISO.append(succ_node)
                            stack.append((succ_node,path+[succ_node],SISO))
                    else:
                        if plan_flag==1:
                            if plan_score_dict[succ_node]>0:
                                if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                    if not succ_node in SISO:
                                        SISO.append(succ_node)
                                stack.append((succ_node,path+[succ_node],SISO))
                        elif plan_flag>1:
                            if succ_node==max_score_suc:
                                if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                    if not succ_node in SISO:
                                        SISO.append(succ_node)
                                stack.append((succ_node,path+[succ_node],SISO))
                        else:
                            if PE_plan_flag==1:
                                if PE_plan_score_dict[succ_node]>0:
                                    if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                        if not succ_node in SISO:
                                            SISO.append(succ_node)
                                    stack.append((succ_node,path+[succ_node],SISO))
                            elif PE_plan_flag>1:
                                if succ_node==max_score_suc:
                                    if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                        if not succ_node in SISO:
                                            SISO.append(succ_node)
                                    stack.append((succ_node,path+[succ_node],SISO))
                            else:
                                if path_flag==1:
                                    if path_score_dict[succ_node]>0:
                                        if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                            if not succ_node in SISO:
                                                SISO.append(succ_node)
                                        stack.append((succ_node,path+[succ_node],SISO))
                                elif path_flag>1:
                                    if succ_node==max_score_suc:
                                        if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                            if not succ_node in SISO:
                                                SISO.append(succ_node)
                                        stack.append((succ_node,path+[succ_node],SISO))
                                else:
                                    if path_plan_flag==1:
                                        if path_plan_score_dict[succ_node]>0:   ## bug fixed
                                            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                                if not succ_node in SISO:
                                                    SISO.append(succ_node)
                                            stack.append((succ_node,path+[succ_node],SISO))
                                    elif path_plan_flag>1:
                                        if succ_node==max_score_suc:
                                            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                                if not succ_node in SISO:
                                                    SISO.append(succ_node)
                                            stack.append((succ_node,path+[succ_node],SISO))
                                    else:
                                        if max_score_suc!=0 and succ_node==max_score_suc:
                                            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                                if not succ_node in SISO:
                                                    SISO.append(succ_node)
                                            stack.append((succ_node,path+[succ_node],SISO))
                                        else:
                                            #pdb.set_trace()
                                            flag+=1
            if flag==len(succ_node_all):
                #pdb.set_trace()
                print "Break middle node found!",path,vertex
                yield path

        elif len(succ_node_all)==1: # just one successor, append to the path
            #if SISO_flag or plan_flag or PE_plan_flag or path_flag or path_plan_flag:
            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                if not succ_node in SISO:
                    SISO.append(succ_node)
            if G.out_degree(succ_node)==0:
                yield path+[succ_node]
            else:
                stack.append((succ_node,path+[succ_node],SISO))
        else: # isolated node
            yield path

def DFS_paths_single_pair_end(G, start_node, PE_G, read_db, fragment_len=300):
    # PE_node_dict: the index of group for each node, key: node, value: index, group: pair-end connected groups
    stack=[(start_node,[start_node],[start_node])]
    while stack:
        (vertex,path,SISO)=stack.pop()
        succ_node_all = set(G.successors(vertex)) - set(path)

        flag=0 #indicate whether the node is appended
        SISO_flag,plan_flag,PE_plan_flag,path_flag,path_plan_flag=0,0,0,0,0
        SISO_score_dict,plan_score_dict,PE_plan_score_dict,path_score_dict,path_plan_score_dict={},{},{},{},{}
        for succ_node in succ_node_all:
            plan_nodes=[succ_node]
            while G.out_degree(plan_nodes[-1])==1 and G.out_degree(G.successors(plan_nodes[-1])[0])==1 and (not G.successors(plan_nodes[-1])[0] in plan_nodes): 
                plan_nodes.append(G.successors(plan_nodes[-1])[0])
            #if len(PE_G.predecessors(succ_node))==0 and len(plan_nodes)==1: #if succ_node has no predecessors in PE_G graph, consider its successors
            if len(plan_nodes)==1: #if cannot plan to succ_node's successors, plan_flag will be same as SISO_flag, consider its successors directly
                plan_nodes=G.successors(succ_node)
            ##
            plan_nodes=list(set(plan_nodes)-set(path))
            SISO_score=get_paired_score_2(SISO,succ_node,PE_G)
            SISO_plan_score=get_paired_score_plan(SISO,plan_nodes,PE_G)
            #PE_plan_nodes=list(set(PE_G.successors(succ_node))-set(path))
            PE_plan_nodes=list(set(PE_G[succ_node].keys())-set(path))
            PE_plan_score=get_paired_score_plan(SISO,PE_plan_nodes,PE_G)
            path_score=get_paired_score_2(path[:-1],succ_node,PE_G)
            path_plan_score=get_paired_score_plan(path[:-1],plan_nodes,PE_G)
            if SISO_score>0:
                SISO_flag+=1
            if SISO_plan_score>0:
                plan_flag+=1
            if PE_plan_score>0:
                PE_plan_flag+=1
            if path_score>0:
                path_flag+=1
            if path_plan_score>0:
                path_plan_flag+=1
            SISO_score_dict[succ_node]=SISO_score
            plan_score_dict[succ_node]=SISO_plan_score
            PE_plan_score_dict[succ_node]=PE_plan_score
            path_score_dict[succ_node]=path_score
            path_plan_score_dict[succ_node]=path_plan_score
        #pdb.set_trace()
        max_score_suc=0
        if SISO_flag>=1:
            #max_score_suc=get_max_score_node(SISO_score_dict)
            max_score_suc=classify_SISO_score_node(SISO_score_dict,plan_score_dict,PE_plan_score_dict,path_score_dict,path_plan_score_dict)
        elif plan_flag>=1:
            max_score_suc=get_max_score_node(plan_score_dict)
        elif PE_plan_flag>=1:
            max_score_suc=get_max_score_node(PE_plan_score_dict)
        elif path_flag>=1:
            max_score_suc=get_max_score_node(path_score_dict)
        elif path_plan_flag>=1:
            max_score_suc=get_max_score_node(path_plan_score_dict)
        else:  ## all the flags are zero
            #pdb.set_trace()
            max_score_suc=bifurcation_classifier(vertex,G,PE_G,path,SISO,read_db,fragment_len)
      
        if len(succ_node_all)<G.out_degree(vertex):
            print "Finding a cycle:",path
            yield path
            continue

        if len(succ_node_all)>1: # seeing a bifurcation node, multiple successors, judge which one to append
            for succ_node in succ_node_all:
                if G.out_degree(succ_node)==0: # reaching the end, bug fixed
                    if SISO_score>0:
                        yield path+[succ_node]
                    elif len(path)>1 and succ_node in PE_G[path[-2]]:
                        yield path+[succ_node]
                    elif len(path)==1 and succ_node in PE_G[path[-1]]:
                        yield path+[succ_node]
                    else:
                        flag+=1
                        ## remove the flag from the ending node
                        if SISO_flag>0 and SISO_score_dict[succ_node]>0:
                            SISO_score_dict[succ_node]=0
                            SISO_flag-=1
                        if plan_flag>0 and plan_score_dict[succ_node]>0:
                            plan_score_dict[succ_node]=0
                            plan_flag-=1
                        if PE_plan_flag>0 and PE_plan_score_dict[succ_node]>0:
                            PE_plan_score_dict[succ_node]=0
                            PE_plan_flag-=1
                        if path_flag>0 and path_score_dict[succ_node]>0:
                            path_score_dict[succ_node]=0
                            path_flag-=1
                        if path_plan_flag>0 and path_plan_score_dict[succ_node]>0:
                            path_plan_score_dict[succ_node]=0
                            path_plan_flag-=1
                else:                    # not reaching the end
                    if SISO_flag==1:
                        if SISO_score_dict[succ_node]>0:
                            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                if not succ_node in SISO:
                                    SISO.append(succ_node)
                            stack.append((succ_node,path+[succ_node],SISO))
                    elif SISO_flag>1:
                        if succ_node==max_score_suc:
                            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                if not succ_node in SISO:
                                    SISO.append(succ_node)
                            stack.append((succ_node,path+[succ_node],SISO))
                    else:
                        if plan_flag==1:
                            if plan_score_dict[succ_node]>0:
                                if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                    if not succ_node in SISO:
                                        SISO.append(succ_node)
                                stack.append((succ_node,path+[succ_node],SISO))
                        elif plan_flag>1:
                            if succ_node==max_score_suc:
                                if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                    if not succ_node in SISO:
                                        SISO.append(succ_node)
                                stack.append((succ_node,path+[succ_node],SISO))
                        else:
                            if PE_plan_flag==1:
                                if PE_plan_score_dict[succ_node]>0:
                                    if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                        if not succ_node in SISO:
                                            SISO.append(succ_node)
                                    stack.append((succ_node,path+[succ_node],SISO))
                            elif PE_plan_flag>1:
                                if succ_node==max_score_suc:
                                    if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                        if not succ_node in SISO:
                                            SISO.append(succ_node)
                                    stack.append((succ_node,path+[succ_node],SISO))
                            else:
                                if path_flag==1:
                                    if path_score_dict[succ_node]>0:
                                        if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                            if not succ_node in SISO:
                                                SISO.append(succ_node)
                                        stack.append((succ_node,path+[succ_node],SISO))
                                elif path_flag>1:
                                    if succ_node==max_score_suc:
                                        if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                            if not succ_node in SISO:
                                                SISO.append(succ_node)
                                        stack.append((succ_node,path+[succ_node],SISO))
                                else:
                                    if path_plan_flag==1:
                                        if path_plan_score_dict[succ_node]>0:   ## bug fixed
                                            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                                if not succ_node in SISO:
                                                    SISO.append(succ_node)
                                            stack.append((succ_node,path+[succ_node],SISO))
                                    elif path_plan_flag>1:
                                        if succ_node==max_score_suc:
                                            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                                if not succ_node in SISO:
                                                    SISO.append(succ_node)
                                            stack.append((succ_node,path+[succ_node],SISO))
                                    else:
                                        if max_score_suc!=0 and succ_node==max_score_suc:
                                            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                                if not succ_node in SISO:
                                                    SISO.append(succ_node)
                                            stack.append((succ_node,path+[succ_node],SISO))
                                        else:
                                            #pdb.set_trace()
                                            flag+=1
            if flag==len(succ_node_all):
                #pdb.set_trace()
                print "Break middle node found!",path,vertex
                for succ_node in succ_node_all:
                    if G.in_degree(succ_node)==1:
                        if G.out_degree(succ_node)==1:
                            if not succ_node in SISO:
                                SISO.append(succ_node)
                            stack.append((succ_node, path+[succ_node], SISO))
                        elif G.out_degree(succ_node)==0:
                            yield path+[succ_node]
                        else:
                            stack.append((succ_node, path+[succ_node], SISO))

        elif len(succ_node_all)==1: # just one successor, append to the path
            #if SISO_flag or plan_flag or PE_plan_flag or path_flag or path_plan_flag:
            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                if not succ_node in SISO:
                    SISO.append(succ_node)
            if G.out_degree(succ_node)==0:
                yield path+[succ_node]
            else:
                stack.append((succ_node,path+[succ_node],SISO))
        else: # isolated node
            yield path

def DFS_paths_single_pair_end_unassembled(G, start_node, PE_G, read_db, fragment_len=300):
    # PE_node_dict: the index of group for each node, key: node, value: index, group: pair-end connected groups
    stack=[(start_node,[start_node],[start_node])]
    while stack:
        (vertex,path,SISO)=stack.pop()
        succ_node_all = set(G.successors(vertex)) - set(path)

        flag=0 #indicate whether the node is appended
        SISO_flag,plan_flag,PE_plan_flag,path_flag,path_plan_flag=0,0,0,0,0
        SISO_score_dict,plan_score_dict,PE_plan_score_dict,path_score_dict,path_plan_score_dict={},{},{},{},{}
        for succ_node in succ_node_all:
            plan_nodes=[succ_node]
            while G.out_degree(plan_nodes[-1])==1 and G.out_degree(G.successors(plan_nodes[-1])[0])==1 and (not G.successors(plan_nodes[-1])[0] in plan_nodes): 
                plan_nodes.append(G.successors(plan_nodes[-1])[0])
            #if len(PE_G.predecessors(succ_node))==0 and len(plan_nodes)==1: #if succ_node has no predecessors in PE_G graph, consider its successors
            if len(plan_nodes)==1: #if cannot plan to succ_node's successors, plan_flag will be same as SISO_flag, consider its successors directly
                plan_nodes=G.successors(succ_node)
            ##
            plan_nodes=list(set(plan_nodes)-set(path))
            SISO_score=get_paired_score_2(SISO,succ_node,PE_G)
            SISO_plan_score=get_paired_score_plan(SISO,plan_nodes,PE_G)
            #PE_plan_nodes=list(set(PE_G.successors(succ_node))-set(path))
            PE_plan_nodes=list(set(PE_G[succ_node].keys())-set(path))
            PE_plan_score=get_paired_score_plan(SISO,PE_plan_nodes,PE_G)
            path_score=get_paired_score_2(path[:-1],succ_node,PE_G)
            path_plan_score=get_paired_score_plan(path[:-1],plan_nodes,PE_G)
            if SISO_score>0:
                SISO_flag+=1
            if SISO_plan_score>0:
                plan_flag+=1
            if PE_plan_score>0:
                PE_plan_flag+=1
            if path_score>0:
                path_flag+=1
            if path_plan_score>0:
                path_plan_flag+=1
            SISO_score_dict[succ_node]=SISO_score
            plan_score_dict[succ_node]=SISO_plan_score
            PE_plan_score_dict[succ_node]=PE_plan_score
            path_score_dict[succ_node]=path_score
            path_plan_score_dict[succ_node]=path_plan_score
        #pdb.set_trace()
        max_score_suc=0
        if SISO_flag>=1:
            #max_score_suc=get_max_score_node(SISO_score_dict)
            max_score_suc=classify_SISO_score_node(SISO_score_dict,plan_score_dict,PE_plan_score_dict,path_score_dict,path_plan_score_dict)
        elif plan_flag>=1:
            max_score_suc=get_max_score_node(plan_score_dict)
        elif PE_plan_flag>=1:
            max_score_suc=get_max_score_node(PE_plan_score_dict)
        elif path_flag>=1:
            max_score_suc=get_max_score_node(path_score_dict)
        elif path_plan_flag>=1:
            max_score_suc=get_max_score_node(path_plan_score_dict)
        else:  ## all the flags are zero
            #pdb.set_trace()
            max_score_suc=bifurcation_classifier(vertex,G,PE_G,path,SISO,read_db,fragment_len)
      
        if len(succ_node_all)<G.out_degree(vertex):
            print "Finding a cycle:",path
            yield path
            continue

        if len(succ_node_all)>1: # seeing a bifurcation node, multiple successors, judge which one to append
            for succ_node in succ_node_all:
                if G.out_degree(succ_node)==0: # reaching the end, bug fixed
                    if SISO_score>0:
                        yield path+[succ_node]
                    elif len(path)>1 and succ_node in PE_G[path[-2]]:
                        yield path+[succ_node]
                    elif len(path)==1 and succ_node in PE_G[path[-1]]:
                        yield path+[succ_node]
                    else:
                        flag+=1
                        ## remove the flag from the ending node
                        if SISO_flag>0 and SISO_score_dict[succ_node]>0:
                            SISO_score_dict[succ_node]=0
                            SISO_flag-=1
                        if plan_flag>0 and plan_score_dict[succ_node]>0:
                            plan_score_dict[succ_node]=0
                            plan_flag-=1
                        if PE_plan_flag>0 and PE_plan_score_dict[succ_node]>0:
                            PE_plan_score_dict[succ_node]=0
                            PE_plan_flag-=1
                        if path_flag>0 and path_score_dict[succ_node]>0:
                            path_score_dict[succ_node]=0
                            path_flag-=1
                        if path_plan_flag>0 and path_plan_score_dict[succ_node]>0:
                            path_plan_score_dict[succ_node]=0
                            path_plan_flag-=1
                else:                    # not reaching the end
                    if SISO_flag==1:
                        if SISO_score_dict[succ_node]>0:
                            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                if not succ_node in SISO:
                                    SISO.append(succ_node)
                            stack.append((succ_node,path+[succ_node],SISO))
                    elif SISO_flag>1:
                        if succ_node==max_score_suc:
                            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                if not succ_node in SISO:
                                    SISO.append(succ_node)
                            stack.append((succ_node,path+[succ_node],SISO))
                    else:
                        if plan_flag==1:
                            if plan_score_dict[succ_node]>0:
                                if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                    if not succ_node in SISO:
                                        SISO.append(succ_node)
                                stack.append((succ_node,path+[succ_node],SISO))
                        elif plan_flag>1:
                            if succ_node==max_score_suc:
                                if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                    if not succ_node in SISO:
                                        SISO.append(succ_node)
                                stack.append((succ_node,path+[succ_node],SISO))
                        else:
                            if PE_plan_flag==1:
                                if PE_plan_score_dict[succ_node]>0:
                                    if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                        if not succ_node in SISO:
                                            SISO.append(succ_node)
                                    stack.append((succ_node,path+[succ_node],SISO))
                            elif PE_plan_flag>1:
                                if succ_node==max_score_suc:
                                    if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                        if not succ_node in SISO:
                                            SISO.append(succ_node)
                                    stack.append((succ_node,path+[succ_node],SISO))
                            else:
                                if path_flag==1:
                                    if path_score_dict[succ_node]>0:
                                        if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                            if not succ_node in SISO:
                                                SISO.append(succ_node)
                                        stack.append((succ_node,path+[succ_node],SISO))
                                elif path_flag>1:
                                    if succ_node==max_score_suc:
                                        if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                            if not succ_node in SISO:
                                                SISO.append(succ_node)
                                        stack.append((succ_node,path+[succ_node],SISO))
                                else:
                                    if path_plan_flag==1:
                                        if path_plan_score_dict[succ_node]>0:   ## bug fixed
                                            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                                if not succ_node in SISO:
                                                    SISO.append(succ_node)
                                            stack.append((succ_node,path+[succ_node],SISO))
                                    elif path_plan_flag>1:
                                        if succ_node==max_score_suc:
                                            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                                if not succ_node in SISO:
                                                    SISO.append(succ_node)
                                            stack.append((succ_node,path+[succ_node],SISO))
                                    else:
                                        if max_score_suc!=0 and succ_node==max_score_suc:
                                            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                                                if not succ_node in SISO:
                                                    SISO.append(succ_node)
                                            stack.append((succ_node,path+[succ_node],SISO))
                                        else:
                                            #pdb.set_trace()
                                            flag+=1

            # For firstly unassembled nodes, revise this part to make sure extend or yield
            if flag==len(succ_node_all):
                #pdb.set_trace()
                print "Break middle node found!",path,vertex
                flag_in_degree = 0
                for succ_node in succ_node_all:
                    if G.in_degree(succ_node)==1:
                        flag_in_degree+=1
                        if G.out_degree(succ_node)==1:
                            if not succ_node in SISO:
                                SISO.append(succ_node)
                            stack.append((succ_node, path+[succ_node], SISO))
                        elif G.out_degree(succ_node)==0:
                            yield path+[succ_node]
                        else:
                            stack.append((succ_node, path+[succ_node], SISO))
                if (not flag_in_degree) and (not stack): # in case cannot extend, yield the current path
                    yield path


        elif len(succ_node_all)==1: # just one successor, append to the path
            #if SISO_flag or plan_flag or PE_plan_flag or path_flag or path_plan_flag:
            if G.in_degree(succ_node)==1 and G.out_degree(succ_node)==1:
                if not succ_node in SISO:
                    SISO.append(succ_node)
            if G.out_degree(succ_node)==0:
                yield path+[succ_node]
            else:
                stack.append((succ_node,path+[succ_node],SISO))
        else: # isolated node
            pdb.set_trace()
            yield path

def get_assemblie2(G,read_db):
    contigs={}
    if len(G.nodes())>1:
        starting_nodes=[n for n in G.nodes() if G.in_degree(n)==0]
        ending_nodes=[n for n in G.nodes() if G.out_degree(n)==0]

        paths=[]
        for start_node in starting_nodes:
            for end_node in ending_nodes:
                two_nodes_paths=[]
                for path in DFS_paths_interative:
                    two_nodes_paths.append(path)

                for path in two_nodes_paths:
                    print path
                    contig_key='contig_'+':'.join(path)
                    contigs[contig_key]=read_db[path[0]]
                    for idx in range(1,len(path)):
                        prev,current=path[idx-1],path[idx]
                        seq=read_db[current]
                        #pdb.set_trace()
                        overlap=int(G[prev][current]["label"])
                        contigs[contig_key]+=seq[overlap:]
                    #contigs.append(contig)
    else:
        contig_key='contig_'+G.nodes()[0]
        contigs[contig_key]=read_db[G.nodes()[0]]

    return contigs

# assemble with knonw paths
def get_assemblie3(G,paths,read_db):
    contigs={}
    if len(G.nodes())>1:
        for path in paths:
            #print path
            contig_key='contig_'+':'.join(path)
            contigs[contig_key]=read_db[path[0]]
            for idx in range(1,len(path)):
                prev,current=path[idx-1],path[idx]
                seq=read_db[current]
                overlap=int(G[prev][current][0]["label"])
                contigs[contig_key]+=seq[overlap:]
            #contigs.append(contig)
    else:
        contig_key='contig_'+G.nodes()[0]
        contigs[contig_key]=read_db[G.nodes()[0]]
    return contigs

def contains_sublist(lst,sublst):
    return '--'.join(sublst) in '--'.join(lst)

def group_contain_sublist(lsts,sublst):
    flag=0
    for lst in lsts:
        if contains_sublist(lst,sublst):
            flag=1
            break
    return flag
