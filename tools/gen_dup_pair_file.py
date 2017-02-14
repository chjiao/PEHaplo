import re,sys,pdb

'''
For single-end reads in the rmdup file, find the sequence of the other end if it is duplicated, link this read to that read in the non-dup file
Algorithm:
    for  each single-end reads, get the sequence of its pair, if the sequence is in the rmdup file, keep this pair
'''

rmdup_file=sys.argv[1]  # the duplicated number of these reads is at least N
raw_rmdup_file=sys.argv[2]
pair1_file=sys.argv[3]
pair2_file=sys.argv[4]
f_out=open('pair_end_connections.txt','w')

rmdup_single_end_reads={}
single_seq={}
count=0
with open(rmdup_file,'r') as f:
    for line in f:
        line=line.strip()
        if line.startswith('>'):
            count+=1
            title=line.split()[0]
            base=title.split('/')[0]
        else:
            seq=line
            if not base in rmdup_single_end_reads:
                rmdup_single_end_reads[base]=title
                single_seq[base]=seq
            else:
                rmdup_single_end_reads.pop(base)
                single_seq.pop(base)
print "Single-end reads number: %d." % len(rmdup_single_end_reads)
print "Total reads number: %d." % count

raw_rmdup_seq={}
with open(raw_rmdup_file,'r') as f:
    for line in f:
        line=line.strip()
        if line.startswith('>'):
            title=line.split()[0]
        else:
            #assert line in raw_rmdup_seq, "Duplicated reads found!"
            if line in raw_rmdup_seq:
                pdb.set_trace()
            raw_rmdup_seq[line]=title

unconcordant1 = 0
unconcordant2 = 0
pair_count=0
with open(pair1_file,'r') as f:
    for line in f:
        line=line.strip()
        if line.startswith('>'):
            title=line.split()[0]
            base=title.split('/')[0]
        else:
            seq=line
            if base in rmdup_single_end_reads and title!=rmdup_single_end_reads[base]:
                if seq in raw_rmdup_seq: # duplicated sequence
                    pair_count+=1
                    f_out.write(raw_rmdup_seq[seq][1:]+"\t"+rmdup_single_end_reads[base][1:]+'\n')
                    if raw_rmdup_seq[seq].endswith('2'):
                        unconcordant1+=1
                    #    pdb.set_trace()
                    #    print raw_rmdup_seq[seq]

with open(pair2_file,'r') as f:
    for line in f:
        line=line.strip()
        if line.startswith('>'):
            title=line.split()[0]
            base=title.split('/')[0]
        else:
            seq=line
            if base in rmdup_single_end_reads and title!=rmdup_single_end_reads[base]:  # meet with a pair
                if seq in raw_rmdup_seq: # duplicated sequence
                    pair_count+=1
                    f_out.write(rmdup_single_end_reads[base][1:]+'\t'+raw_rmdup_seq[seq][1:]+'\n')
                    if raw_rmdup_seq[seq].endswith('1'):
                        unconcordant2+=1
                    #    pdb.set_trace()
                    #    print raw_rmdup_seq[seq]
print "%d single-end reads have been found their pairs's sequences." % pair_count
print "Unconcordant duplicated pair on .1 file: %d; on .2 file: %d." % (unconcordant1, unconcordant2)
f_out.close()


