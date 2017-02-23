import re,sys,pdb,subprocess
import math
import networkx as nx
import pygraphviz
import numpy as np
import matplotlib
matplotlib.use('Agg')
from scipy import stats
from matplotlib import pyplot as plt
import matplotlib.patches as patches

# use the midpoint instead of the whole profile
# imply the PECC methods to identify misjoin contigs
# 2016.09.20
# save all the contig fragments
# 2016.10.05

def cal_contig_pair_end_connections(read_ids_1,read_ids_2):
    read1_dict={}
    read2_dict={}
    pairs1=[]
    pairs2=[]
    for read_id in read_ids_1:
        base=read_id.split('/')[0]
        if not base in read1_dict:
            read1_dict[base]=read_id
        else:
            read1_dict.pop(base)
    for read_id in read_ids_2:
        base=read_id.split('/')[0]
        if not base in read2_dict:
            read2_dict[base]=read_id
        else:
            read2_dict.pop(base)

    ol_1,ol_2=0,0
    for base in read1_dict:
        if base in read2_dict:
            read1=read1_dict[base]
            read2=read2_dict[base]
            if read1.endswith('1') and read2.endswith('2'):
                ol_1+=1
                pairs1.append([read1,read2])
            elif read1.endswith('2') and read2.endswith('1'):
                ol_2+=1
                pairs2.append([read2,read1])
    return ol_1,ol_2,pairs1,pairs2

def find_interval(A): #profile_array
    interval=[]
    A[A>1]=1
    A_diff=np.diff(A)
    des_loc=np.where(A_diff==-1)[0]
    asc_loc=np.where(A_diff==1)[0]
    #pdb.set_trace()
    for loc in des_loc:
        asc_pair=asc_loc[asc_loc>loc]
        if len(asc_pair)>0:
            interval.append([loc+1,asc_pair[0]])
    return interval

def cluster_intervals(con_invs,cutoff):
    cluster_invs=[]
    inv_cluster=[]
    for i in range(len(con_invs)):
        new_inv=con_invs[i]
        if len(inv_cluster)==0:
            inv_cluster=new_inv
        elif new_inv[0]-inv_cluster[1]<cutoff:
            inv_cluster[1]=new_inv[1]
        else:
            cluster_invs.append(inv_cluster)
            inv_cluster=con_invs[i]
    if len(inv_cluster)>0:
        cluster_invs.append(inv_cluster)
    return cluster_invs

def cal_discordant_interval(read_ids, reads_loc, single_reads_dict, inv, con_len, fragment_len, read_len):
    """
    read_ids: reads aligned to the contig
    reads_loc: strand and locations of reads aligned to the contig
    single_reads_dict: single-end reads aligned to the contig
    """
    pl_num,pr_num,sl_num,sr_num = 0,0,0,0
    pairs_len_l=[]
    pairs_len_r=[]
    left_margin=max(0,int(inv[0]-fragment_len/2))
    right_margin=min(con_len,int(inv[1]+fragment_len/2))
    for read_id in read_ids:
        base=read_id.split('/')[0]
        if base in single_reads_dict: # not paired reads
            #pdb.set_trace()
            if read_id.endswith('1') and int(reads_loc[read_id])>=left_margin and int(reads_loc[read_id])<inv[0]:
                sl_num+=1
            elif read_id.endswith('2') and int(reads_loc[read_id])>inv[1] and int(reads_loc[read_id])<right_margin:
                sr_num+=1
        else: # paired reads
            if read_id.endswith('1') and int(reads_loc[read_id])>=left_margin and int(reads_loc[read_id])<inv[0]:
                pl_num+=1
                pair=read_id.split('/')[0]+'/2'
                pair_len=int(reads_loc[pair])-int(reads_loc[read_id])+read_len+1
                pairs_len_l.append(pair_len)
            elif read_id.endswith('2') and int(reads_loc[read_id])>inv[1] and int(reads_loc[read_id])<right_margin:
                pr_num+=1
                pair=read_id.split('/')[0]+'/1'
                if not read_id in reads_loc or (not pair in reads_loc):
                    pdb.set_trace()
                pair_len=int(reads_loc[read_id])-int(reads_loc[pair])+read_len+1
                pairs_len_r.append(pair_len)
    return sl_num,sr_num,pl_num,pr_num,pairs_len_l,pairs_len_r

def cal_discordant_interval_2(read_ids, reads_loc, single_reads_dict, inv, con_len, fragment_len, read_len):
    """
    read_ids: reads aligned to the contig
    reads_loc: strand and locations of reads aligned to the contig
    single_reads_dict: single-end reads aligned to the contig
    Changes: for either .1 or .2 reads, if the alignment strand is plus, this read should be the left read and its corresponding pair will be right read
    """
    pl_num,pr_num,sl_num,sr_num = 0,0,0,0
    pairs_len_l=[]
    pairs_len_r=[]
    left_margin=max(0,int(inv[0]-fragment_len/2))
    right_margin=min(con_len,int(inv[1]+fragment_len/2))
    for read_id in read_ids:
        strand, loc = reads_loc[read_id]
        base=read_id.split('/')[0]
        if base in single_reads_dict: # not paired reads
            if strand=='0' and int(loc)>=left_margin and int(loc)<inv[0]: # left reads
                sl_num+=1
            elif strand=='16' and int(loc)>inv[1] and int(loc)<right_margin: # right reads
                sr_num+=1
        else: # paired reads
            if read_id.endswith('1'):
                pair = read_id.split('/')[0]+'/2'
                pair_loc = reads_loc[pair][1]
            elif read_id.endswith('2'):
                pair = read_id.split('/')[0]+'/1'
                pair_loc = reads_loc[pair][1]
            else:
                print "Read_id is not pair.1 or pair.2: ", read_id
            if strand=='0' and int(loc)>=left_margin and int(loc)<inv[0]:
                pl_num+=1
                pair_len=abs(int(pair_loc)-int(loc))+read_len+1
                pairs_len_l.append(pair_len)
            elif strand=='16' and int(loc)>inv[1] and int(loc)<right_margin:
                pr_num+=1
                pair_len=abs(int(pair_loc)-int(loc))+read_len+1
                pairs_len_r.append(pair_len)
    return sl_num,sr_num,pl_num,pr_num,pairs_len_l,pairs_len_r

def cal_discordant_interval_3(read_ids, reads_loc, single_reads_dict, left_pairs, right_pairs, inv, con_len, fragment_len, read_len):
    """
    read_ids: reads aligned to the contig
    reads_loc: strand and locations of reads aligned to the contig
    single_reads_dict: single-end reads aligned to the contig
    Changes: for either .1 or .2 reads, if the alignment strand is plus, this read should be the left read and its corresponding pair will be right read
    """
    pl_num,pr_num,sl_num,sr_num = 0,0,0,0
    pairs_len_l=[]
    pairs_len_r=[]
    left_margin=max(0,int(inv[0]-fragment_len/2))
    right_margin=min(con_len,int(inv[1]+fragment_len/2))
    for read_id in read_ids:
        strand, loc = reads_loc[read_id]
        base=read_id.split('/')[0]
        if base in single_reads_dict: # not paired reads
            if strand=='0' and int(loc)>=left_margin and int(loc)<inv[0]: # left reads
                sl_num+=1
            elif strand=='16' and int(loc)>inv[1] and int(loc)<right_margin: # right reads
                sr_num+=1
        else: # paired reads
            if read_id in left_pairs:
                pair = left_pairs[read_id]
            elif read_id in right_pairs:
                pair = right_pairs[read_id]
            pair_strand, pair_loc = reads_loc[pair]
            if strand=='0' and int(loc)>=left_margin and int(loc)<inv[0]:
                pl_num+=1
                pair_len=abs(int(pair_loc)-int(loc))+read_len+1
                pairs_len_l.append(pair_len)
            elif strand=='16' and int(loc)>inv[1] and int(loc)<right_margin:
                pr_num+=1
                pair_len=abs(int(pair_loc)-int(loc))+read_len+1
                pairs_len_r.append(pair_len)
    return sl_num,sr_num,pl_num,pr_num,pairs_len_l,pairs_len_r

def cal_SR(reads_ids,reads_loc,single_reads_dict,inv,con_len,fragment_len,read_len):
    pairs_len=[]
    for read_id in read_ids:
        strand, loc = reads_loc[read_id]
        base=read_id.split('/')[0]
        if read_id.endswith('1'):
            pair = base+'/2'
        elif read_id.endswith('2'):
            pair = base+'/1'
        else:
            pdb.set_trace()
        if not base in single_reads_dict: # paired reads
            pair_strand, pair_loc = reads_loc[pair]
            if strand=='16' and int(loc)>=inv[0] and int(loc)<=inv[1]:
                pair_len=abs(int(loc)-int(pair_loc))+read_len+1
                pairs_len.append(pair_len)
    return pairs_len

def cal_SL(reads_ids,reads_loc,single_reads_dict,inv,con_len,fragment_len,read_len):
    pairs_len=[]
    for read_id in read_ids:
        strand, loc = reads_loc[read_id]
        base=read_id.split('/')[0]
        if read_id.endswith('1'):
            pair = base+'/2'
        elif read_id.endswith('2'):
            pair = base+'/1'
        else:
            pdb.set_trace()
        if not base in single_reads_dict: # paired reads
            pair_strand, pair_loc = reads_loc[pair]
            if strand=='0' and int(loc)>=inv[0] and int(loc)<=inv[1]:
                pair_len=abs(int(loc)-int(pair_loc))+read_len+1
                pairs_len.append(pair_len)
    return pairs_len

def get_wide_clip(clip_locs):
    # get top 3 wide inv
    loc_len = [loc[1]-loc[0] for loc in clip_locs]
    loc_width_idx = [i[0] for i in sorted(enumerate(loc_len), reverse=True, key = lambda x:x[1])]
    if (len(loc_width_idx)>3):
        loc_width_idx = loc_width_idx[:3]
    loc_width_idx = sorted(loc_width_idx)
    wide_clip_locs = []
    for idx in loc_width_idx:
        if loc_len[idx]>10:
            wide_clip_locs.append(clip_locs[idx])
    return wide_clip_locs

############################################################################################
contigs_file = sys.argv[1]
sam_file = sys.argv[2]  # the reads that are aligned

# parameters 
read_len = int(sys.argv[3])
Fragment_len = int(sys.argv[4])
std = int(sys.argv[5])
mu = Fragment_len


#Fragment_len=500
#read_len=200
#mu=500
#std=200

lineno=0
list_idx=0
contig_idx={} # key: contig title value: contig index
contigs_dict={} #key: contig index, value: contig sequence
contig_len_list=[]

reads_all_loc=[]

contig_read_loc={} # contains dictionary for each contig, key: aligned reads locations, strand and pair-end, value: read ids
reads_loc=[] # contains dictionary for each contig, key: reads_ids; value: aligned locations

with open(contigs_file,'r') as f:
    for line in f:
        lineno+=1
        if lineno%2==1:
            title=line[1:-1]
            title = title.split()[0]
            #pdb.set_trace()
        elif lineno%2==0:
            seq=line.strip()
            contigs_dict[list_idx]=seq
            if not title in contig_idx:
                contig_idx[title]=list_idx
                contig_len_list.append(len(seq))
                contig_read_loc[list_idx]={}
                reads_loc.append({})
                reads_all_loc.append({})
                list_idx+=1

# read pairs
"""
left_pair_dict={}
right_pair_dict={}
with open(pair_file, 'r') as f:
    for line in f:
        pair1, pair2 = line.strip().split()
        if not pair1 in left_pair_dict:
            left_pair_dict[pair1] = pair2
        else:
            pdb.set_trace()
        if not pair2 in right_pair_dict:
            right_pair_dict[pair2] = pair1
        else:
            pdb.set_trace()
"""

# save all the aligned results
with open(sam_file,'r') as f:
    for line in f:
        if not line.startswith('@'):
            lmap=line.split()
            strand=lmap[1]
            loc=lmap[3]
            #pair_end=lmap[0].split('/')[-1]
            if lmap[2] in contig_idx:
                reads_all_loc[contig_idx[lmap[2]]][lmap[0]]=(strand,loc)

#pdb.set_trace()
contig_read_ids={} # save the read ids for each contig
# look at the reads that both ends are aligned
unaligned = 0
with open(sam_file,'r') as f:
    for line in f:
        if not line.startswith('@'):
            lmap=line.split()
            strand=lmap[1]
            loc=lmap[3]
            pair_end=lmap[0].split('/')[-1]
            if lmap[2] in contig_idx: # the aligned contig
                con_index=contig_idx[lmap[2]]
                if not (pair_end,strand,loc) in contig_read_loc[contig_idx[lmap[2]]]:
                    ## look at the pair
                    if lmap[0].endswith('1'):
                        base=lmap[0].split('/')[0]
                        read_pair=base+'/2'
                        if read_pair in reads_all_loc[contig_idx[lmap[2]]]:
                            (pair_strand,pair_loc)=reads_all_loc[con_index][read_pair]
                            if not ('2',pair_strand,pair_loc) in contig_read_loc[con_index]: # both ends are aligned and not duplicated, save them simultaniously
                                reads_loc[con_index][lmap[0]] = (strand, loc)
                                reads_loc[con_index][read_pair] = (pair_strand, pair_loc)
                                if not con_index in contig_read_ids:
                                    contig_read_ids[con_index]=[lmap[0]]
                                    contig_read_ids[con_index].append(read_pair)
                                else:
                                    contig_read_ids[con_index].append(lmap[0])
                                    contig_read_ids[con_index].append(read_pair)
                                contig_read_loc[con_index][(pair_end,strand,loc)]=lmap[0]
                                contig_read_loc[con_index][('2',pair_strand,pair_loc)]=read_pair
                    else:
                        base=lmap[0].split('/')[0]
                        read_pair=base+'/1'
                        if read_pair in reads_all_loc[contig_idx[lmap[2]]]:
                            (pair_strand, pair_loc)=reads_all_loc[con_index][read_pair]
                            if not ('1',pair_strand,pair_loc) in contig_read_loc[con_index]: # both ends are aligned and not duplicated, save them simultaniously
                                reads_loc[con_index][lmap[0]] = (strand, loc)
                                reads_loc[con_index][read_pair] = (pair_strand, pair_loc)
                                if not con_index in contig_read_ids:
                                    contig_read_ids[con_index]=[lmap[0]]
                                    contig_read_ids[con_index].append(read_pair)
                                else:
                                    contig_read_ids[con_index].append(lmap[0])
                                    contig_read_ids[con_index].append(read_pair)
                                contig_read_loc[con_index][(pair_end,strand,loc)]=lmap[0]
                                contig_read_loc[con_index][('1',pair_strand,pair_loc)]=read_pair
            else:
                unaligned+=1

print "Unaligned reads number: %d." % unaligned


#pdb.set_trace()
# look at the reads that are only one end is aligned
with open(sam_file,'r') as f:
    for line in f:
        if not line.startswith('@'):
            lmap=line.split()
            strand=lmap[1]
            loc=lmap[3]
            pair_end=lmap[0].split('/')[-1]
            if lmap[2] in contig_idx: # the aligned contig
                con_index=contig_idx[lmap[2]]
                if not (pair_end,strand,loc) in contig_read_loc[contig_idx[lmap[2]]]:
                    ## look at the pair
                    if lmap[0].endswith('1'):
                        base=lmap[0].split('/')[0]
                        read_pair=base+'/2'
                        if not read_pair in reads_all_loc[contig_idx[lmap[2]]]: # the pair is not aligned
                            reads_loc[con_index][lmap[0]]= (strand, loc)
                            if not con_index in contig_read_ids:
                                contig_read_ids[con_index]=[lmap[0]]
                            else:
                                contig_read_ids[con_index].append(lmap[0])
                            contig_read_loc[con_index][(pair_end,strand,loc)]=lmap[0]
                    else:
                        base=lmap[0].split('/')[0]
                        read_pair=base+'/1'
                        if not read_pair in reads_all_loc[contig_idx[lmap[2]]]: # the pair is not aligned
                            reads_loc[con_index][lmap[0]]= (strand, loc)
                            if not con_index in contig_read_ids:
                                contig_read_ids[con_index]=[lmap[0]]
                            else:
                                contig_read_ids[con_index].append(lmap[0])
                            contig_read_loc[con_index][(pair_end,strand,loc)]=lmap[0]

## clip the contigs and output the results
f_out=open('Contigs_clipped.fa','w')
con_len_range=max(contig_len_list)

for i in range(len(contig_read_ids)):
    read_ids=contig_read_ids[i]
    con_len=contig_len_list[i]
    profile_array=np.zeros(con_len)
    con_cov=len(read_ids)*read_len/float(con_len)
    con_cov=round(con_cov,2)
    con_name='Contig_'+str(i+1)+'_'+str(con_len)+'_'+str(con_cov)

    reads_dict={}
    contig_pairs=[]
    count_pair=0
    pair_distance=0
    for read_id in read_ids: # look at the reads aligned to one contig
        base=read_id.split('/')[0]
        if not base in reads_dict:
            reads_dict[base]=read_id
        else:
            if read_id.endswith('1'):
                strand1, align1 = reads_loc[i][read_id][0],int(reads_loc[i][read_id][1])
                strand2, align2 = reads_loc[i][reads_dict[base]][0],int(reads_loc[i][reads_dict[base]][1])
                count_pair+=1
                pair_distance += abs(align2-align1)+read_len
                contig_pairs.append((read_id,reads_dict[base]))
            elif read_id.endswith('2'):
                strand1, align1 = reads_loc[i][reads_dict[base]][0], int(reads_loc[i][reads_dict[base]][1])
                strand2, align2 = reads_loc[i][read_id][0], int(reads_loc[i][read_id][1])
                count_pair+=1
                pair_distance += abs(align2-align1)+read_len
                contig_pairs.append((reads_dict[base],read_id))
            reads_dict.pop(base)
    
    ## calculate the interval cutoff
    inv_cut=-math.log(0.01)*con_len/float(count_pair)
    print i,inv_cut
    ## step1: find the intervals
    for pair in contig_pairs:
        strand1, align1=reads_loc[i][pair[0]][0], int(reads_loc[i][pair[0]][1])
        strand2, align2=reads_loc[i][pair[1]][1], int(reads_loc[i][pair[1]][1])
        mid_loc=int((align1+align2)/2)
        profile_array[mid_loc]+=1

    # plot the profile
    #plt.plot(profile_array, 'k')
    #plt.xlabel('Contig (bp)', fontsize = 12)
    #plt.ylabel('read pair midpoints coverage', fontsize = 12)
    #plt.savefig(con_name+'.pdf',format = 'pdf')
    
    con_intervals=find_interval(profile_array)
    tmp_intervals=[]
    for inv in con_intervals:
        if inv[1]-inv[0]+1>=max(10,inv_cut/2.0): # the interval is set to be at least 10 bp
            tmp_intervals.append(inv)
    inv_clusters=cluster_intervals(tmp_intervals,10)
    tmp_intervals2=[]
    for inv in inv_clusters:
        if inv[1]-inv[0]+1>=max(20, inv_cut): # the clustered interval is set to be at least 20 bp
            tmp_intervals2.append(inv)
    print i,tmp_intervals2
   
    ## step2: validate the misjoin regions
    tmp_intervals3=[]
    l_discordant,r_discordant=0,0
    for inv in tmp_intervals2:
        inv_temp=[]
        inv_temp[:]=inv
        sl,sr,pl,pr,pairs_len_l,pairs_len_r=cal_discordant_interval_2(read_ids,reads_loc[i],reads_dict,inv,con_len,Fragment_len,read_len)
        #pdb.set_trace()
        l_discordant, r_discordant = 0,0
        if pl>0:
            l_discordant=float(sl)/pl
        if pr>0:
            r_discordant=float(sr)/pr
        if l_discordant>0.5: # hypothesis testing
            if len(pairs_len_l)==0:
                pdb.set_trace()
            D,pValue=stats.kstest(pairs_len_l,'norm',(mu,std))
            if pValue<0.01:
                inv_temp[0]=max(0, inv[0]-Fragment_len/2)
        if pl==0 and sl>0:
            inv_temp[0] = max(0, inv[0]-Fragment_len/2)
        if r_discordant>0.5:
            if len(pairs_len_r)==0:
                pdb.set_trace()
            D,pValue=stats.kstest(pairs_len_r,'norm',(mu,std))
            if pValue<0.01:
                inv_temp[1]= min(con_len, inv[1]+Fragment_len/2)
        if pr==0 and sr>0:
            inv_temp[1] = min(con_len, inv[1] + Fragment_len/2)
        tmp_intervals3.append(inv_temp)
    print i,tmp_intervals3
    
    #pdb.set_trace()

    ## step3: identify the clipping boundary
    clip_loc=[]
    for inv in tmp_intervals3:
        ds=Fragment_len/2
        left_margin=min(inv[1],inv[0]+ds)
        while (ds>10):
            SR_left=cal_SR(read_ids,reads_loc[i],reads_dict,[inv[0],left_margin],con_len,Fragment_len,read_len)
            #pdb.set_trace()
            if len(SR_left)==0:
                break
            D,pValue=stats.kstest(SR_left,'norm',(mu,std))
            #pdb.set_trace()
            if pValue>=0.01:
                break
            else:
                ds=ds/2
                left_margin=min(inv[1],inv[0]+ds)

        ds=Fragment_len/2 
        right_margin=max(inv[0],inv[1]-ds)
        while (ds>10):
            SL_right=cal_SL(read_ids,reads_loc[i],reads_dict,[right_margin,inv[1]],con_len,Fragment_len,read_len)
            if len(SL_right)==0:
                break
            D,pValue=stats.kstest(SL_right,'norm',(mu,std))
            if pValue>=0.01:
                break
            else:
                ds=ds/2
                right_margin=max(inv[0],inv[1]-ds)
        
        if left_margin<right_margin: # cut the contig
            clip_loc.append([left_margin,right_margin])
        else:
            cut_loc=(left_margin+right_margin)/2
            clip_loc.append([cut_loc,cut_loc])

    # clip the contig
    #pdb.set_trace()
    if len(clip_loc)>0:
        seq=contigs_dict[i]
        if len(clip_loc)==1:
            seq1=seq[0:clip_loc[0][0]]
            seq2=seq[clip_loc[0][1]:]
            con_name1='Contig_'+str(i+1)+'-1'+'_'+str(len(seq1))
            con_name2='Contig_'+str(i+1)+'-2'+'_'+str(len(seq2))
            if len(seq1)>read_len:
                f_out.write('>'+con_name1+'\n'+seq1+'\n')
            if len(seq2)>read_len:
                f_out.write('>'+con_name2+'\n'+seq2+'\n')
        else:
            #pdb.set_trace()
            clip_loc = get_wide_clip(clip_loc)
            print i, clip_loc
            if not clip_loc:
                f_out.write('>'+con_name+'\n'+seq+'\n')
                continue
            clip_parts = []
            for k in range(len(clip_loc)):
                if k ==0:
                    clip_parts.append([0,clip_loc[k][0]])
                else:
                    if clip_loc[k-1][1] < clip_loc[k][0]:
                        clip_parts.append([clip_loc[k-1][1],clip_loc[k][0]])
            clip_parts.append([clip_loc[-1][1],len(seq)])
            print i, clip_parts
                
            clip_idx=1
            for clip in clip_parts:
                clip_seq = seq[clip[0]:clip[1]]
                clip_con_name = 'Contig_'+str(i+1)+'_'+str(clip_idx)+'_'+str(len(clip_seq))
                f_out.write('>'+clip_con_name+'\n'+clip_seq+'\n')
                clip_idx+=1

    else: # output the original sequence
        seq=contigs_dict[i]
        f_out.write('>'+con_name+'\n'+seq+'\n')
f_out.close()
