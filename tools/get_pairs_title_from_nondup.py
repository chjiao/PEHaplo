import re,sys,pdb
from __future__ import division

# the output file of sga contains single-end and pair-end reads
# only keep the pair-end reads 
nondup_file=sys.argv[1]
out_file=sys.argv[2]
f_out=open(out_file,'w')

title_base_dict={}
reads_count=0
title_dict={}
with open(nondup_file,'r') as f:
    for line in f:
        if line.startswith('>'):
            title=line.split()[0]
            base=title.split('/')[0]
        else:
            seq=line.strip()
            reads_count+=1
            if not base in title_base_dict:
                title_base_dict[base]=1
                title_dict[base]=title
            else:
                title_base_dict[base]+=1
                f_out.write(title_dict[base][1:]+'\t'+title[1:]+'\n')
                assert title_dict[base].endswith('1') and title.endswith('2')

unpaired_list=[x for x in title_base_dict.keys() if title_base_dict[x]==1]
paired_list=[x for x in title_base_dict.keys() if title_base_dict[x]==2]
print "Number of total reads:\t",2*len(paired_list)+len(unpaired_list)
print "Number of paired reads:\t",2*len(paired_list)
print "Number of unpaired reads:\t",len(unpaired_list)

f_out.close()

