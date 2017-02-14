import re,sys,pdb

fa_file=sys.argv[1]

f1=open('fa1.fa','w')
f2=open('fa2.fa','w')
title_dict={}
title_list=[]

lineno=0
with open(fa_file,'r') as f:
    for line in f:
        lineno+=1
        if lineno%4==1:
            title1=line.strip()
        elif lineno%4==2:
            seq1=line.strip()
        elif lineno%4==3:
            title2=line.strip()
        elif lineno%4==0:
            seq2=line.strip()
            f1.write(title1+'\n'+seq1+'\n')
            f2.write(title2+'\n'+seq2+'\n')
f1.close()
f2.close()
