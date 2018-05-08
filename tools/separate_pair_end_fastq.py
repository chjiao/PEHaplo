import re,sys,pdb
from __future__ import division

fq_file=sys.argv[1]
f1=open('separate1.fq','w')
f2=open('separate2.fq','w')

title_dict={}
title_list=[]
lineno=0
with open(fq_file,'r') as f:
    for line in f:
        lineno+=1
        line=line.strip()
        if lineno%4==1:
            title=line
        elif lineno%4==2:
            seq=line
        elif lineno%4==3:
            qtitle=line
        elif lineno%4==0:
            quality=line
            if title.endswith('1'):
                f1.write(title+'\n'+seq+'\n'+qtitle+'\n'+quality+'\n')
            else:
                f2.write(title+'\n'+seq+'\n'+qtitle+'\n'+quality+'\n')

f1.close()
f2.close()
