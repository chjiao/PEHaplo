import sys
from __future__ import division


fa_file=sys.argv[1]
num = int(sys.argv[2])
f_out=open(sys.argv[3],'w')
total_num=0
kept_num=0
with open(fa_file,'r') as f:
    for line in f:
        line=line.strip()
        if line.startswith('>'):
            total_num+=1
            title=line
            dup_num=title.split('NumDuplicates=')[1]
            dup_num=int(dup_num)
        else:
            seq=line
            if dup_num>=num:
                kept_num+=1
                f_out.write(title+'\n'+seq+'\n')

print "Total reads number:\t",total_num
print "Kept reads number:\t",kept_num
f_out.close()
