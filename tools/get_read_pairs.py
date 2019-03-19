import re,sys,pdb

fa_file = sys.argv[1]
f0 = open("single_end.fa", 'w')
f1 = open("pair1.fa", 'w')
f2 = open("pair2.fa", 'w')

fa_dict = {}
seq_dict = {}
cnt = 0
title, seq = "", ""
with open(fa_file, 'r') as f:
    for line in f:
        if line.startswith('>'):
            if seq: # process the previous sequence
                base = title.split('/')[0]
                if not base in fa_dict:
                    fa_dict[base] = title
                    seq_dict[base] = seq
                else:
                    cnt+=1
                    if title.endswith('1'):
                        f1.write('>'+title+'\n'+seq+'\n')
                        f2.write('>'+fa_dict[base]+'\n'+seq_dict[base]+'\n')
                    else:
                        f1.write('>'+fa_dict[base]+'\n'+seq_dict[base]+'\n')
                        f2.write('>'+title+'\n'+seq+'\n')
                    fa_dict.pop(base)
                    seq_dict.pop(base)
            seq = ""
            title = line[1:-1]
        else:
            seq += line[:-1]
    
    if seq: # the last sequence
        base = title.split('/')[0]
        if not base in fa_dict:
            fa_dict[base] = title
            seq_dict[base] = seq
        else:
            cnt+=1
            if title.endswith('1'):
                f1.write('>'+title+'\n'+seq+'\n')
                f2.write('>'+fa_dict[base]+'\n'+seq_dict[base]+'\n')
            else:
                f1.write('>'+fa_dict[base]+'\n'+seq_dict[base]+'\n')
                f2.write('>'+title+'\n'+seq+'\n')
            fa_dict.pop(base)
            seq_dict.pop(base)


seq, base = "", ""
with open(fa_file, 'r') as f:
    for line in f:
        if line.startswith('>'):
            if seq and (base in fa_dict):
                f0.write('>'+title+'\n'+seq+'\n')
            seq = ""
            title = line[1:-1]
            base = title.split('/')[0]
        else:
            seq += line[:-1]
    if seq and (base in fa_dict):
        f0.write('>'+title+'\n'+seq+'\n')

f0.close()
f1.close()
f2.close()

print "The number of read pairs is: %d, single-end reads is: %d" % (cnt, len(fa_dict))

