import re,sys,pdb
import networkx as nx
import subprocess
import pygraphviz
from __future__ import division

def get_seq_from_fa(fa_file, des_file):
    """
    read_map, # key: read_name, value: read_index
    read_name_list, # read names
    seq_dict, # key: read_index, value, sequence
    """

    read_map={} 
    read_name_list=[] 
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
                read_name=line[1:].split()[0]
                if read_name in read_map:
                    read_idx=read_map[read_name]
                    seq_dict[read_idx]=""
            else:
                if read_name in read_map:
                    read_idx = read_map[read_name]
                    seq_dict[read_idx]+=line[:-1]
    return read_name_list,read_map, seq_dict

def create_graph(des_list, edge_file):

    """ 
    create the initial graph, node name: read_index
    des_list: the list storing the read_name
    read_node_dict: store the corresponding node for each read, key: read_index, value: corresponding node
    """

    G = nx.MultiDiGraph()
    read_node_dict = {}
    unconcordant_edge = 0
    with open(edge_file,'r') as f:
        for line in f:
            if "-" in line:
                unconcordant_edge+=1
                continue
            
            read_1, read_2, overlap_len=line.strip().split(' + ')
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

        if unconcordant_edge>0:
            print "Unconcordant overlap found:%d."%(unconcordant_edge)
    return G, read_node_dict 

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
            #pair1_dict[idx1]=idx2
            #pair2_dict[idx2]=idx1
    return pair_dict


def get_seq_from_fa_mix(fa_file, des_file):
    """
    read_map, # key: read_name, value: read_index
    read_name_list, # read names
    seq_dict, # key: read_index, value, sequence
    """

    read_map_single = {}
    read_map_pair = {}
    read_single = {}
    read_pair = {}
    read_name_list = []
    count=0
    with open(des_file,'r') as f:
        for line in f:
            if not line.find('/')>-1: # single-end reads
                read_map_single[line[:-1]]=str(count)
                read_single[str(count)] = line[:-1]
                read_name_list.append(line[:-1])
                count+=1
            else:
                break
        print line

        read_map_pair[line[:-1]] = str(count)
        read_pair[str(count)] = line[:-1]
        read_name_list.append(line[:-1])
        count += 1

        for line in f:
            read_map_pair[line[:-1]] = str(count)
            read_pair[str(count)] = line[:-1]
            read_name_list.append(line[:-1])
            count+=1
    
    not_found=0
    seq_dict={}  # key: read_index, value: the corresponding sequence
    with open(fa_file,'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                read_name=line[1:].split()[0]
                if read_name in read_map_single or read_name in read_map_pair:
                    if read_name in read_map_single:
                        read_idx=read_map_single[read_name]
                    else:
                        read_idx = read_map_pair[read_name]
                    seq_dict[read_idx]=""
                else:
                    not_found+=1

            else:
                if read_name in read_map_single or read_name in read_map_pair:
                    if read_name in read_map_single:
                        read_idx=read_map_single[read_name]
                    else:
                        read_idx = read_map_pair[read_name]
                    seq_dict[read_idx]+=line
    pdb.set_trace()
    print not_found
    print len(read_map_single), len(read_map_pair), len(seq_dict)
    return read_map_single, read_map_pair, read_single, read_pair, read_name_list, seq_dict

def create_graph_mix(des_list, edge_file, single_end_reads_num):
    """
    create graph from both singe-end reads and pair-end reads
    des_list: the list storing the read_name
    read_node_dict: store the corresponding node for each read, key: read_index, value: corresponding node
    """

    G = nx.MultiDiGraph()
    read_node_dict = {}
    unconcordant_edge = 0
    with open(edge_file,'r') as f:
        for line in f:
            if "-" in line:
                unconcordant_edge+=1
                continue
            
            read_1, read_2, overlap_len=line.strip().split(' + ')
            if not read_1 in G:
                G.add_node(read_1,read_ids=[read_1])
            read_node_dict[read_1] = read_1

            if not read_2 in G[read_1] and (not read_1==read_2): # bug fixed, judge whether the edge already exist, remove the edge connecting to self
                G.add_node(read_2,read_ids=[read_2])
                G.add_edge(read_1,read_2,label=overlap_len)
                read_node_dict[read_2] = read_2

            elif read_2 in G[read_1]:
                print "Duplicate edge found!",line.strip()
                if int(overlap_len)>int(G[read_1][read_2][0]['label']):
                    G[read_1][read_2][0]['label']=overlap_len

        if unconcordant_edge>0:
            print "Unconcordant overlap found:%d."%(unconcordant_edge)
    return G, read_node_dict 


def get_sub_reads(edge_file, read_db):
    sub_reads=set([])
    with open(edge_file,'r') as f:
        for line in f:
            if "-" in line:
                continue
            lmap=line.strip().split()
            # 9005 + 38127 + 230
            if lmap[0]==lmap[2]:
                continue
            #max_read_len=max(len(read_db[lmap[0]]),len(read_db[lmap[2]]))
            if int(lmap[4])==len(read_db[lmap[0]]):
                sub_reads.add(lmap[0])
            elif int(lmap[4])==len(read_db[lmap[2]]):
                sub_reads.add(lmap[2])
    return sub_reads

def compare_list(list1,list2):
    return sorted(list1)==sorted(list2)


