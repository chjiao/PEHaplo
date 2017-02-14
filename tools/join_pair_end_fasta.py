import re,sys,pdb

# join two paired-end fq files to one fastq file
fa1=sys.argv[1]
fa2=sys.argv[2]
join_fa=sys.argv[3]

f_out=open(join_fa,'w')
with open(fa1, 'r') as f1, open(fa2, 'r') as f2:
    for line1, line2 in zip(f1, f2):
        line1 = line1.strip()
        line2 = line2.strip()
        if line1.startswith('>') and line2.startswith('>'):
            title1 = line1
            title2 = line2
        else:
            seq1 = line1
            seq2 = line2
            f_out.write(title1+'\n'+seq1+'\n'+title2+'\n'+seq2+'\n')
f_out.close()

