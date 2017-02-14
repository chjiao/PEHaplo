import re,sys,os,pdb,subprocess
import argparse

# -------------------------------
__author__ = "Jiao Chen"
__date__ = "10 Feb, 2017"
__license__ = "GPL"

usage = """
Program: PEHaplo - De novo haplotype reconstruction in virus quasispecies using paired-end reads
Version: 0.1
Contact: Jiao Chen - chjiao3456@gmail.com
"""

def main():
    parser = argparse.ArgumentParser(description = usage)
    #basic = parser.add_argument('-f', dest='input_file', type=str, help='input fastq or fasta file')
    basic = parser.add_argument('-f1', dest='input_f1', type=str, required=True, help='input .1 part of paired-end fastq or fasta file')
    basic = parser.add_argument('-f2', dest='input_f2', type=str, required=True, help='input .2 part of paired-end fastq or fasta file')
    basic = parser.add_argument('-l', '--overlap_len', dest='overlap_len', type=int, required=True, help='overlap threshold between reads for reads orientation adjustment')
    basic = parser.add_argument('-n', '--dup_n', dest='dup_n', type=int, help='the reads kept should be duplicated at least n times')

    if len(sys.argv[1:])==0:
        print usage
        parser.exit()
    args = parser.parse_args()
    
    base_path = os.path.dirname(os.path.abspath(__file__))

    if args.input_f1 and args.input_f2:
        # ---------------- Preprocess ------------------

        ## error correction
        #"""
        subprocess.check_call("karect -correct -inputfile=%s -inputfile=%s -celltype=haploid -matchtype=hamming -aggressive=5.0 -numstages=2 -errorrate=0.25 -errorratesec=0.25" % (args.input_f1, args.input_f2), shell=True)
        subprocess.check_call('python "%s"/tools/join_pair_end_fasta.py karect_%s karect_%s karect_whole.fa' % (base_path, args.input_f1, args.input_f2), shell=True)

        ## remove duplicates and substrings
        subprocess.check_call("sga preprocess karect_whole.fa >karect_whole_preprocessed.fa", shell=True)
        subprocess.check_call("sga index karect_whole_preprocessed.fa", shell=True)
        subprocess.check_call("sga rmdup karect_whole_preprocessed.fa", shell=True)
        #"""

        ## keep reads duplicated n times
        if args.dup_n:
            subprocess.check_call('python "%s"/tools/get_dup_num_reads.py karect_whole_preprocessed.rmdup.fa %s kept_num.fa' %(base_path, args.dup_n), shell=True)
        else:
            subprocess.check_call("cp karect_whole_preprocessed.rmdup.fa kept_num.fa", shell=True)
        
        ## readjoiner for overlap graph
        subprocess.check_call("gt readjoiner prefilter -q -des -readset samp -db kept_num.fa", shell=True)
        subprocess.check_call("gt readjoiner overlap -readset samp -l %s" % args.overlap_len, shell=True)
        subprocess.check_call("gt readjoiner spmtest -readset samp.0 -test showlist >samp_edges.txt", shell=True)
        subprocess.check_call("rm samp.0.spm samp.rlt samp.ssp samp.0.cnt samp.esq samp.sds", shell=True)

        ## reads orientation adjustment
        subprocess.check_call('python "%s"/tools/get_pairs_title_from_nondup.py kept_num.fa kept_pairs.txt' % base_path, shell=True)
        subprocess.check_call('python "%s"/tools/Seperate_strand_single_pair.py samp.des samp_edges.txt kept_num.fa kept_pairs.txt' % base_path, shell=True)

        # ---------------- Assembly ------------------
        subprocess.check_call('python "%s"/tools/gen_dup_pair_file.py kept_num.fa karect_whole_preprocessed.rmdup.fa karect_%s karect_%s' %(base_path, args.input_f1, args.input_f2), shell=True)
        subprocess.check_call('python "%s"/apsp_overlap_clique_flu.py Plus_strand_reads.fa pair_end_connections.txt %s' % (base_path, args.overlap_len), shell=True)
                

if __name__=='__main__':
    main()
