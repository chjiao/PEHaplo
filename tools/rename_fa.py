import re,sys,pdb

fa_file=sys.argv[1]
f_out=open(sys.argv[2],'w')
pair_id = int(sys.argv[3])
assert pair_id==1 or pair_id==2

lineno=0
with open(fa_file,'r') as f:
    for line in f:
        lineno+=1
        if lineno%2==1:
            line=line.strip()
            new_title=line.split()[0]+'/'+str(pair_id)
            f_out.write(new_title+'\n')
        elif lineno%2==0:
            f_out.write(line)
f_out.close()

