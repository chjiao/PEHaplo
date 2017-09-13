import re,sys,pdb
import networkx as nx
import subprocess
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
import pygraphviz
from Bio.Seq import Seq


def reverse_complement(seq):
    seq=Seq(seq)
    return str(seq.reverse_complement())

def get_seq_from_fa(fa_file,des_file):
    read_map={}  # key: read_name, value: read_index
    read_name_list=[] # read names
    count=0
    with open(des_file,'r') as f:
        for line in f:
            # HCV_1-163200/1
            read_map[line[:-1]]=str(count)
            read_name_list.append(line[:-1])
            count+=1

    seq_dict={}  # key: read_index, value: the corresponding sequence
    with open(fa_file,'r') as f:
        for line in f:
            if line.startswith('>'):
                read_id=line[1:].split()[0]
                if read_id in read_map:
                    read_num=read_map[read_id]
                    seq_dict[read_num]=""
            else:
                if read_id in read_map:
                    read_num=read_map[read_id]
                    seq_dict[read_num]+=line[:-1]
    return read_name_list,read_map,seq_dict

def read_pair_file(pair_file, read_map):
    pair1_dict, pair2_dict = {},{}
    with open(pair_file,'r') as f:
        for line in f:
            line=line.strip()
            pair1, pair2 = line.split()
            idx1,idx2 = read_map[pair1],read_map[pair2]
            pair1_dict[idx1]=idx2
            pair2_dict[idx2]=idx1
    return pair1_dict, pair2_dict

def get_seq_from_fq(fq_file,des_file):
    read_map={}  # key: read_name, value: read_index
    read_name_list=[] # read names
    count=0
    with open(des_file,'r') as f:
        for line in f:
            # HCV_1-163200/1
            read_map[line[:-1]]=str(count)
            read_name_list.append(line[:-1])
            count+=1

    seq_dict={}  # key: read_index, value: the corresponding sequence and quality
    lineno=0
    with open(fq_file,'r') as f:
        for line in f:
            lineno+=1
            if lineno%4==1:
                title=line[1:].split()[0]
            elif lineno%4==2:
                seq=line.strip()
            elif lineno%4==3:
                title2=line.strip()
            elif lineno%4==0:
                quality=line.strip()
                if title in read_map:
                    read_num=read_map[title]
                    seq_dict[read_num]=(seq,quality)
    return read_name_list,read_map,seq_dict

def create_graph_with_fq(edge_file,des_list):
    # create the initial graph, node name: read_index
    # read_node_dict: store the corresponding node for each read
    # des_list: the list storing the read_name

    G = nx.DiGraph()

    with open(edge_file,'r') as f:
        for line in f:
            # 23287 + 40135 - 196
            read_1, std1, read_2, std2, overlap_len=line.strip().split()
            connect_type=std1+std2
            #read_name=des_list[int(read_1)]
            G.add_node(read_1)
            #read_node_dict[read_name]=read_1
            #read_name=des_list[int(read_2)]
            #G.add_node(read_2)
            #read_node_dict[read_name]=read_2
            if not read_2 in G[read_1] and (not read_1==read_2): # bug fixed, judge whether the edge already exist, remove the edge connecting to self
                G.add_edge(read_1,read_2,label=overlap_len,con_type=connect_type)
            elif read_2 in G[read_1]:
                print "Duplicate edge found!",line.strip()
                if int(overlap_len)>int(G[read_1][read_2]['label']):
                    G[read_1][read_2]['label']=overlap
                #pdb.set_trace()
    return G 

def BFS_tranverse_graph(G,read_db,des_list,start_node):
    # G: a directed, connected graph, garantee that the sequences can be divided into two groups

    vertex_strand='+' #+: get the sequence as in the original dataset; '-': get the sequence as the reverse complement in the original dataset
    visited,queue=set(),[(start_node,vertex_strand)]
    plus_seq_dict={}
    minus_seq_dict={}
    plus_seq_dict[start_node]=read_db[start_node]
    minus_seq_dict[start_node]=reverse_complement(read_db[start_node])

    while queue:
        vertex,vertex_strand=queue.pop(0)
        if vertex not in visited:
            visited.add(vertex)
            #queue.extend(set(G.successors(vertex))-visited)
            candidate_nodes=G.successors(vertex)
            candidate_nodes.extend(G.predecessors(vertex))
            for succ_node in candidate_nodes:
                if not succ_node in plus_seq_dict:
                    ## get the edge strands
                    if succ_node in G[vertex]: # successor
                        edge_strands=G[vertex][succ_node]['con_type']
                    elif vertex in G[succ_node]: # predecessor
                        edge_strands=G[succ_node][vertex]['con_type']

                    if vertex_strand=='+':
                        #edge_strands=G[vertex][succ_node]['con_type']
                        if edge_strands[0]==edge_strands[1]:
                            plus_seq_dict[succ_node]=read_db[succ_node]
                            minus_seq_dict[succ_node]=reverse_complement(read_db[succ_node])
                            queue.append((succ_node,'+'))
                        else:
                            plus_seq_dict[succ_node]=reverse_complement(read_db[succ_node]) # reverse complement
                            minus_seq_dict[succ_node]=read_db[succ_node]
                            queue.append((succ_node,'-'))
                    elif vertex_strand=='-':
                        if edge_strands[0]==edge_strands[1]:
                            minus_seq_dict[succ_node]=read_db[succ_node]
                            plus_seq_dict[succ_node]=reverse_complement(read_db[succ_node])
                            queue.append((succ_node,'-'))
                        else:
                            minus_seq_dict[succ_node]=reverse_complement(read_db[succ_node]) # reverse complement
                            plus_seq_dict[succ_node]=read_db[succ_node]
                            queue.append((succ_node,'+'))
        
    ## look at the pair
    unvisited_pair=set()
    for node in visited:
        node=int(node)
        if node%2==0:  # the pair.1
            if str(node+1) in G:
                if not str(node+1) in visited:
                    if plus_seq_dict[str(node)]==read_db[str(node)]:
                        unvisited_pair.add((str(node+1),'-'))
                        plus_seq_dict[str(node+1)]=reverse_complement(read_db[str(node+1)])
                        minus_seq_dict[str(node+1)]=read_db[str(node+1)]
                    elif plus_seq_dict[str(node)]==reverse_complement(read_db[str(node)]):
                        unvisited_pair.add((str(node+1),'+'))
                        plus_seq_dict[str(node+1)]=read_db[str(node+1)]
                        minus_seq_dict[str(node+1)]=reverse_complement(read_db[str(node+1)])
            else:
                #print "Not paired read found!"
                if plus_seq_dict[str(node)] == read_db[str(node)]:
                    plus_seq_dict[str(node+1)]=reverse_complement(read_db[str(node+1)])
                    minus_seq_dict[str(node+1)]=read_db[str(node+1)]
                elif plus_seq_dict[str(node)]==reverse_complement(read_db[str(node)]):
                    plus_seq_dict[str(node+1)]=read_db[str(node+1)]
                    minus_seq_dict[str(node+1)]=reverse_complement(read_db[str(node+1)])
        else:  # the pair.2
            if str(node-1) in G:
                if not str(node-1) in visited:
                    if plus_seq_dict[str(node)]==read_db[str(node)]:
                        unvisited_pair.add((str(node-1),'-'))
                        plus_seq_dict[str(node-1)]=reverse_complement(read_db[str(node-1)])
                        minus_seq_dict[str(node-1)]=read_db[str(node-1)]
                    elif plus_seq_dict[str(node)]==reverse_complement(read_db[str(node)]):
                        unvisited_pair.add((str(node-1),'+'))
                        plus_seq_dict[str(node-1)]=read_db[str(node-1)]
                        minus_seq_dict[str(node-1)]=reverse_complement(read_db[str(node-1)])
            else:
                #print "Not paired read found!"
                if plus_seq_dict[str(node)]==read_db[str(node)]:
                    plus_seq_dict[str(node-1)]=reverse_complement(read_db[str(node-1)])
                    minus_seq_dict[str(node-1)]=read_db[str(node-1)]
                elif plus_seq_dict[str(node)]==reverse_complement(read_db[str(node)]):
                    plus_seq_dict[str(node-1)]=read_db[str(node-1)]
                    minus_seq_dict[str(node-1)]=reverse_complement(read_db[str(node-1)])

    while(unvisited_pair):
        queue=[unvisited_pair.pop()]
        while queue:
            vertex,vertex_strand=queue.pop(0)
            if vertex not in visited:
                visited.add(vertex)
                candidate_nodes=G.successors(vertex)
                candidate_nodes.extend(G.predecessors(vertex))
                for succ_node in candidate_nodes:
                    if not succ_node in plus_seq_dict:
                        # get the edge strands
                        if succ_node in G[vertex]:
                            edge_strands=G[vertex][succ_node]['con_type']
                        elif vertex in G[succ_node]:
                            edge_strands=G[succ_node][vertex]['con_type']

                        if vertex_strand=='+':
                            if edge_strands[0]==edge_strands[1]:
                                plus_seq_dict[succ_node]=read_db[succ_node]
                                minus_seq_dict[succ_node]=reverse_complement(read_db[succ_node])
                                queue.append((succ_node,'+'))
                            else:
                                plus_seq_dict[succ_node]=reverse_complement(read_db[succ_node]) # reverse complement
                                minus_seq_dict[succ_node]=read_db[succ_node]
                                queue.append((succ_node,'-'))
                        elif vertex_strand=='-':
                            if edge_strands[0]==edge_strands[1]:
                                minus_seq_dict[succ_node]=read_db[succ_node]
                                plus_seq_dict[succ_node]=reverse_complement(read_db[succ_node])
                                queue.append((succ_node,'-'))
                            else:
                                minus_seq_dict[succ_node]=reverse_complement(read_db[succ_node]) # reverse complement
                                plus_seq_dict[succ_node]=read_db[succ_node]
                                queue.append((succ_node,'+'))

        for node in visited:
            node=int(node)
            if node%2==0:  # the pair.1
                if str(node+1) in G:
                    if not str(node+1) in visited:
                        if plus_seq_dict[str(node)]==read_db[str(node)]:
                            unvisited_pair.add((str(node+1),'-'))
                            plus_seq_dict[str(node+1)]=reverse_complement(read_db[str(node+1)])
                            minus_seq_dict[str(node+1)]=read_db[str(node+1)]
                        elif plus_seq_dict[str(node)]==reverse_complement(read_db[str(node)]):
                            unvisited_pair.add((str(node+1),'+'))
                            plus_seq_dict[str(node+1)]=read_db[str(node+1)]
                            minus_seq_dict[str(node+1)]=reverse_complement(read_db[str(node+1)])
                else:
                    #print "Not paired read found!"
                    if plus_seq_dict[str(node)]==read_db[str(node)]:
                        plus_seq_dict[str(node+1)]=reverse_complement(read_db[str(node+1)])
                        minus_seq_dict[str(node+1)]=read_db[str(node+1)]
                    elif plus_seq_dict[str(node)]==reverse_complement(read_db[str(node)]):
                        plus_seq_dict[str(node+1)]=read_db[str(node+1)]
                        minus_seq_dict[str(node+1)]=reverse_complement(read_db[str(node+1)])
            else:  # the pair.2
                if str(node-1) in G:
                    if not str(node-1) in visited:
                        if plus_seq_dict[str(node)]==read_db[str(node)]:
                            unvisited_pair.add((str(node-1),'-'))
                            plus_seq_dict[str(node-1)]=reverse_complement(read_db[str(node-1)])
                            minus_seq_dict[str(node-1)]=read_db[str(node-1)]
                        elif plus_seq_dict[str(node)]==reverse_complement(read_db[str(node)]):
                            unvisited_pair.add((str(node-1),'+'))
                            plus_seq_dict[str(node-1)]=read_db[str(node-1)]
                            minus_seq_dict[str(node-1)]=reverse_complement(read_db[str(node-1)])
                else:
                    #print "Not paired read found!"
                    if plus_seq_dict[str(node)]==read_db[str(node)]:
                        plus_seq_dict[str(node-1)]=reverse_complement(read_db[str(node-1)])
                        minus_seq_dict[str(node-1)]=read_db[str(node-1)]
                    elif plus_seq_dict[str(node)]==reverse_complement(read_db[str(node)]):
                        plus_seq_dict[str(node-1)]=read_db[str(node-1)]
                        minus_seq_dict[str(node-1)]=reverse_complement(read_db[str(node-1)])
    G.remove_nodes_from(visited)
    return plus_seq_dict,minus_seq_dict

def BFS_tranverse_graph_pair(G,read_db,des_list,start_node, pair1_dict, pair2_dict):
    # G: a directed, connected graph, garantee that the sequences can be divided into two groups
    # pair1_dict, key: pair.1, value: pair.2
    # pair2_dict, key: pair.2, value: pair.1

    vertex_strand='+' #+: get the sequence as in the original dataset; '-': get the sequence as the reverse complement in the original dataset
    visited,queue=set(),[(start_node,vertex_strand)]
    plus_seq_dict={}
    minus_seq_dict={}
    plus_seq_dict[start_node]=read_db[start_node]
    minus_seq_dict[start_node]=reverse_complement(read_db[start_node])

    while queue:
        vertex,vertex_strand=queue.pop(0)
        if vertex not in visited:
            visited.add(vertex)
            #queue.extend(set(G.successors(vertex))-visited)
            candidate_nodes=G.successors(vertex)
            candidate_nodes.extend(G.predecessors(vertex))
            for succ_node in candidate_nodes:
                if not succ_node in plus_seq_dict:
                    ## get the edge strands
                    if succ_node in G[vertex]: # successor
                        edge_strands=G[vertex][succ_node]['con_type']
                    elif vertex in G[succ_node]: # predecessor
                        edge_strands=G[succ_node][vertex]['con_type']

                    if vertex_strand=='+':
                        #edge_strands=G[vertex][succ_node]['con_type']
                        if edge_strands[0]==edge_strands[1]:
                            plus_seq_dict[succ_node]=read_db[succ_node]
                            minus_seq_dict[succ_node]=reverse_complement(read_db[succ_node])
                            queue.append((succ_node,'+'))
                        else:
                            plus_seq_dict[succ_node]=reverse_complement(read_db[succ_node]) # reverse complement
                            minus_seq_dict[succ_node]=read_db[succ_node]
                            queue.append((succ_node,'-'))
                    elif vertex_strand=='-':
                        if edge_strands[0]==edge_strands[1]:
                            minus_seq_dict[succ_node]=read_db[succ_node]
                            plus_seq_dict[succ_node]=reverse_complement(read_db[succ_node])
                            queue.append((succ_node,'-'))
                        else:
                            minus_seq_dict[succ_node]=reverse_complement(read_db[succ_node]) # reverse complement
                            plus_seq_dict[succ_node]=read_db[succ_node]
                            queue.append((succ_node,'+'))
        
    ## look at the pair
    unvisited_pair=set()
    for node in visited:
        if node in pair1_dict:  # the pair.1
            node_pair=pair1_dict[node]
            if node_pair in G:
                if not node_pair in visited:
                    if plus_seq_dict[node]==read_db[node]:
                        unvisited_pair.add((node_pair,'-'))
                        plus_seq_dict[node_pair]=reverse_complement(read_db[node_pair])
                        minus_seq_dict[node_pair]=read_db[node_pair]
                    elif plus_seq_dict[node]==reverse_complement(read_db[node]):
                        unvisited_pair.add((node_pair,'+'))
                        plus_seq_dict[node_pair]=read_db[node_pair]
                        minus_seq_dict[node_pair]=reverse_complement(read_db[node_pair])
            else:
                #print "Not paired read found!"
                if plus_seq_dict[node] == read_db[node]:
                    plus_seq_dict[node_pair]=reverse_complement(read_db[node_pair])
                    minus_seq_dict[node_pair]=read_db[node_pair]
                elif plus_seq_dict[node]==reverse_complement(read_db[node]):
                    plus_seq_dict[node_pair]=read_db[node_pair]
                    minus_seq_dict[node_pair]=reverse_complement(read_db[node_pair])
        elif node in pair2_dict and pair2_dict[node] in G:  # the pair.2
            node_pair=pair2_dict[node]
            if node_pair in G:
                if not node_pair in visited:
                    if plus_seq_dict[node]==read_db[node]:
                        unvisited_pair.add((node_pair,'-'))
                        plus_seq_dict[node_pair]=reverse_complement(read_db[node_pair])
                        minus_seq_dict[node_pair]=read_db[node_pair]
                    elif plus_seq_dict[node]==reverse_complement(read_db[node]):
                        unvisited_pair.add((node_pair,'+'))
                        plus_seq_dict[node_pair]=read_db[node_pair]
                        minus_seq_dict[node_pair]=reverse_complement(read_db[node_pair])
            else:
                #print "Not paired read found!"
                if plus_seq_dict[node]==read_db[node]:
                    plus_seq_dict[node_pair]=reverse_complement(read_db[node_pair])
                    minus_seq_dict[node_pair]=read_db[node_pair]
                elif plus_seq_dict[node]==reverse_complement(read_db[node]):
                    plus_seq_dict[node_pair]=read_db[node_pair]
                    minus_seq_dict[node_pair]=reverse_complement(read_db[node_pair])

    while(unvisited_pair):
        queue=[unvisited_pair.pop()]
        while queue:
            vertex,vertex_strand=queue.pop(0)
            if vertex not in visited:
                visited.add(vertex)
                candidate_nodes=G.successors(vertex)
                candidate_nodes.extend(G.predecessors(vertex))
                for succ_node in candidate_nodes:
                    if not succ_node in plus_seq_dict:
                        # get the edge strands
                        if succ_node in G[vertex]:
                            edge_strands=G[vertex][succ_node]['con_type']
                        elif vertex in G[succ_node]:
                            edge_strands=G[succ_node][vertex]['con_type']

                        if vertex_strand=='+':
                            if edge_strands[0]==edge_strands[1]:
                                plus_seq_dict[succ_node]=read_db[succ_node]
                                minus_seq_dict[succ_node]=reverse_complement(read_db[succ_node])
                                queue.append((succ_node,'+'))
                            else:
                                plus_seq_dict[succ_node]=reverse_complement(read_db[succ_node]) # reverse complement
                                minus_seq_dict[succ_node]=read_db[succ_node]
                                queue.append((succ_node,'-'))
                        elif vertex_strand=='-':
                            if edge_strands[0]==edge_strands[1]:
                                minus_seq_dict[succ_node]=read_db[succ_node]
                                plus_seq_dict[succ_node]=reverse_complement(read_db[succ_node])
                                queue.append((succ_node,'-'))
                            else:
                                minus_seq_dict[succ_node]=reverse_complement(read_db[succ_node]) # reverse complement
                                plus_seq_dict[succ_node]=read_db[succ_node]
                                queue.append((succ_node,'+'))

        for node in visited:
            if node in pair1_dict:  # the pair.1
                node_pair=pair1_dict[node]
                if node_pair in G:
                    if not node_pair in visited:
                        if plus_seq_dict[node]==read_db[node]:
                            unvisited_pair.add((node_pair,'-'))
                            plus_seq_dict[node_pair]=reverse_complement(read_db[node_pair])
                            minus_seq_dict[node_pair]=read_db[node_pair]
                        elif plus_seq_dict[node]==reverse_complement(read_db[node]):
                            unvisited_pair.add((node_pair,'+'))
                            plus_seq_dict[node_pair]=read_db[node_pair]
                            minus_seq_dict[node_pair]=reverse_complement(read_db[node_pair])
                else:
                    #print "Not paired read found!"
                    if plus_seq_dict[node]==read_db[node]:
                        plus_seq_dict[node_pair]=reverse_complement(read_db[node_pair])
                        minus_seq_dict[node_pair]=read_db[node_pair]
                    elif plus_seq_dict[node]==reverse_complement(read_db[node]):
                        plus_seq_dict[node_pair]=read_db[node_pair]
                        minus_seq_dict[node_pair]=reverse_complement(read_db[node_pair])
            elif node in pair2_dict:  # the pair.2
                node_pair=pair2_dict[node]
                if node_pair in G:
                    if not node_pair in visited:
                        if plus_seq_dict[node]==read_db[node]:
                            unvisited_pair.add((node_pair,'-'))
                            plus_seq_dict[node_pair]=reverse_complement(read_db[node_pair])
                            minus_seq_dict[node_pair]=read_db[node_pair]
                        elif plus_seq_dict[node]==reverse_complement(read_db[node]):
                            unvisited_pair.add((node_pair,'+'))
                            plus_seq_dict[node_pair]=read_db[node_pair]
                            minus_seq_dict[node_pair]=reverse_complement(read_db[node_pair])
                else:
                    #print "Not paired read found!"
                    if plus_seq_dict[node]==read_db[node]:
                        plus_seq_dict[node_pair]=reverse_complement(read_db[node_pair])
                        minus_seq_dict[node_pair]=read_db[node_pair]
                    elif plus_seq_dict[node]==reverse_complement(read_db[node]):
                        plus_seq_dict[node_pair]=read_db[node_pair]
                        minus_seq_dict[node_pair]=reverse_complement(read_db[node_pair])
    G.remove_nodes_from(visited)
    return plus_seq_dict,minus_seq_dict


def output_reads(plus_dict, minus_dict,des_list, f1, f2):
    plus_reads=plus_dict.keys()
    plus_reads=[int(x) for x in plus_reads]
    minus_reads=minus_dict.keys()
    minus_reads=[int(x) for x in minus_reads]
    plus_reads=sorted(plus_reads)
    minus_reads=sorted(minus_reads)
    for plus in plus_reads:
        if plus>len(des_list):
            pdb.set_trace()
        f1.write('>'+des_list[plus]+'\n'+str(plus_dict[str(plus)])+'\n')
    for minus in minus_reads:
        if minus>len(des_list):
            pdb.set_trace()
        f2.write('>'+des_list[minus]+'\n'+str(minus_dict[str(minus)])+'\n')

def compare_list(list1,list2):
    return sorted(list1)==sorted(list2)

def plot_graph(G, figname):
    G_plot=nx.drawing.nx_agraph.to_agraph(G)
    G_plot.draw(figname,prog='dot')

###########################################################################
des_file=sys.argv[1]
edge_file=sys.argv[2]
fa_file=sys.argv[3]
pair_file=sys.argv[4]
plus_name='Plus_strand_reads.fa'
minus_name='Minus_strand_reads.fa'
f1=open(plus_name,'w')
f2=open(minus_name,'w')

des_list,read_map,read_db=get_seq_from_fa(fa_file,des_file) # read dictionary
pair1_dict, pair2_dict=read_pair_file(pair_file, read_map)

read_node_dict={}
G=create_graph_with_fq(edge_file,des_list)
print "Graph construction finished!"
#plot_graph(G, 'strand_correction_overlap_graph.png')
#pdb.set_trace()
#subgraphs=nx.weakly_connected_components(G)

#f_out=open('Plus_strand_reads.txt','w')
#f_out2=open('Minus_strand_reads.txt','w')
idx=0
while(len(G)>0):
    idx+=1
    starting_nodes=[n for n in G.nodes() if G.in_degree(n)==0]
    if len(starting_nodes)==0:
        break
        #pdb.set_trace()
    plus_reads_dict,minus_reads_dict=BFS_tranverse_graph_pair(G,read_db,des_list,starting_nodes[0], pair1_dict, pair2_dict)
    if len(plus_reads_dict)>4:
        output_reads(plus_reads_dict, minus_reads_dict, des_list, f1, f2)
f1.close()
f2.close()

