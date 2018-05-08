import re,sys,pdb
from __future__ import division

# join two paired-end fq files to one fastq file
fq1=sys.argv[1]
fq2=sys.argv[2]
join_fq=sys.argv[3]

f_out=open(join_fq,'w')

lineno=0
f1=open(fq1,'r')
f2=open(fq2,'r')

line1=f1.readline()
line2=f2.readline()
while line1:
    lineno+=1
    if lineno%4==1:
        title1=line1.strip()
        title2=line2.strip()
    elif lineno%4==2:
        seq1=line1.strip()
        seq2=line2.strip()
    elif lineno%4==3:
        q_title1=line1.strip()
        q_title2=line2.strip()
    elif lineno%4==0:
        quality1=line1.strip()
        quality2=line2.strip()

        f_out.write(title1+'\n'+seq1+'\n'+q_title1+'\n'+quality1+'\n')
        f_out.write(title2+'\n'+seq2+'\n'+q_title2+'\n'+quality2+'\n')
    line1=f1.readline()
    line2=f2.readline()
f1.close()
f2.close()
f_out.close()

