import re,sys,pdb
from __future__ import division

fq_file=sys.argv[1]
fa_file=sys.argv[2]

lineno=0
fa_out=open(fa_file,'w')

with open(fq_file,'r') as f:
    for line in f:
        lineno+=1
        if lineno%4==1:
            if line.startswith('@'):
                line=re.sub('@','>',line)
                fa_out.write(line)
            else:
                print 'sequence title error!'
                pdb.set_trace()
        elif lineno%4==2:
            fa_out.write(line)

fa_out.close()
