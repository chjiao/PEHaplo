import re,sys,os,pdb,subprocess
import argparse
from __future__ import division

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
    basic = parser.add_argument('-l1', '--overlap_stage1', dest='overlap_len1', type=int, help='overlap cutoff to remove potentially wrong overlaps after merging linked cliques, default same as -l')
    basic = parser.add_argument('-r', '--read_len', dest='read_len', type=int, required=True, help='reads length')
    basic = parser.add_argument('-F', '--fragment_len', dest='fragment_len', type=int, help='paired-end reads insert size, default as 2.5*read_len')
    basic = parser.add_argument('-std', dest='fragment_std', type=int, help='standard deviation of paired-end reads insert size, default as 100') 
    basic = parser.add_argument('-n', '--dup_n', dest='dup_n', type=int, help='the reads kept should be duplicated at least n times, default as keep all the duplicates removed reads')
    basic = parser.add_argument('-correct', dest='contig_correct', type=str, help='whether apply alignment based contigs correction(yes/no), default = no')
    basic = parser.add_argument('-t', '--threads', dest='threads', type=int, help='threads for karect, sga, bowtie2')
    basic = parser.add_argument('-m', '--memory', dest='memory', type=str, help='maximal memory for readjoiner (xGB), default: 2GB')

    if len(sys.argv[1:])==0:
        parser.print_help()        #print usage
        parser.exit()

    args = parser.parse_args()
    
    base_path = os.path.dirname(os.path.abspath(__file__))

    if args.input_f1 and args.input_f2 and args.overlap_len and args.read_len:
        # ---------------- Preprocess ------------------
        if not args.fragment_len:
            args.fragment_len = int(2.5*args.read_len)
        if not args.fragment_std:
            args.fragment_std = 100

        ## error correction
        #"""
        print "Error correction begins -------------------------------------------------"
        if args.threads:
            subprocess.check_call("karect -correct -inputfile=%s -inputfile=%s -threads=%s -celltype=haploid -matchtype=hamming -aggressive=5.0 -numstages=2 -errorrate=0.25 -errorratesec=0.25" % (args.input_f1, args.input_f2, args.threads), shell=True)
        else:
            subprocess.check_call("karect -correct -inputfile=%s -inputfile=%s -celltype=haploid -matchtype=hamming -aggressive=5.0 -numstages=2 -errorrate=0.25 -errorratesec=0.25" % (args.input_f1, args.input_f2), shell=True)

        fa1_file = os.path.basename(args.input_f1)
        fa2_file = os.path.basename(args.input_f2)
        subprocess.check_call('python "%s"/tools/join_pair_end_fasta.py karect_%s karect_%s karect_whole.fa' % (base_path, fa1_file, fa2_file), shell=True)

        ## remove duplicates and substrings
        print "Removing duplicates and substrings -------------------------------------"
        subprocess.check_call("sga preprocess karect_whole.fa >karect_whole_preprocessed.fa", shell=True)
        if args.threads:
            subprocess.check_call("sga index -t %s karect_whole_preprocessed.fa" % args.threads, shell=True)
            subprocess.check_call("sga rmdup -t %s karect_whole_preprocessed.fa" % args.threads, shell=True)
        else:
            subprocess.check_call("sga index karect_whole_preprocessed.fa", shell=True)
            subprocess.check_call("sga rmdup karect_whole_preprocessed.fa", shell=True)

        ## keep reads duplicated n times
        if args.dup_n:
            subprocess.check_call('python "%s"/tools/get_dup_num_reads.py karect_whole_preprocessed.rmdup.fa %s kept_num.fa' %(base_path, args.dup_n), shell=True)
        else:
            subprocess.check_call("cp karect_whole_preprocessed.rmdup.fa kept_num.fa", shell=True)
        
        ## readjoiner for overlap graph
        print "Reads orientation adjustment ------------------------------------------"
        subprocess.check_call("gt readjoiner prefilter -q -des -readset samp -db kept_num.fa", shell=True)
        if args.memory:
            subprocess.check_call("gt readjoiner overlap -readset samp -memlimit %s -l %s" % (args.memory, args.overlap_len), shell=True)
        else: 
            subprocess.check_call("gt readjoiner overlap -readset samp -memlimit 2GB -l %s" % args.overlap_len, shell=True)
        subprocess.check_call("gt readjoiner spmtest -readset samp.0 -test showlist >samp_edges.txt", shell=True)
        #subprocess.check_call("rm samp.0.spm samp.rlt samp.ssp samp.0.cnt samp.esq samp.sds", shell=True)

        ## reads orientation adjustment
        subprocess.check_call('python "%s"/tools/get_pairs_title_from_nondup.py kept_num.fa kept_pairs.txt' % base_path, shell=True)
        subprocess.check_call('python "%s"/tools/Seperate_strand_single_pair.py samp.des samp_edges.txt kept_num.fa kept_pairs.txt' % base_path, shell=True)


        # ---------------- Assembly ------------------
        print "Begin assembly --------------------------------------------------------"
        subprocess.check_call('python "%s"/tools/gen_dup_pair_file.py kept_num.fa karect_whole_preprocessed.rmdup.fa karect_%s karect_%s' %(base_path, fa1_file, fa2_file), shell=True)
        
        if args.overlap_len1:
            subprocess.check_call('python "%s"/apsp_overlap_clique.py Plus_strand_reads.fa pair_end_connections.txt %s %s %s %s' % (base_path, args.overlap_len, args.read_len, args.fragment_len, args.overlap_len1), shell=True)
        else:
            subprocess.check_call('python "%s"/apsp_overlap_clique.py Plus_strand_reads.fa pair_end_connections.txt %s %s %s %s' % (base_path, args.overlap_len, args.read_len, args.fragment_len, args.overlap_len), shell=True)

        # ---------------- remove duplicated contigs ------------------
        if args.threads:
            subprocess.check_call('sga index -t %s Contigs.fa' % args.threads, shell=True)
            subprocess.check_call('sga rmdup -t %s -e 0.005 Contigs.fa' % args.threads, shell=True)
        else:
            subprocess.check_call('sga index Contigs.fa', shell=True)
            subprocess.check_call('sga rmdup -e 0.005 Contigs.fa', shell=True)
        subprocess.check_call('mv Contigs.rmdup.fa Contigs.fa', shell=True)
        #"""

        # ---------------- error correction -------------------------
        if args.contig_correct == 'yes':
            print "Begin contig correction ------------------------------------------- "
            index_path = 'index'
            if not os.path.exists(index_path):
                os.makedirs(index_path)
            subprocess.check_call('bowtie2-build -f Contigs.fa index/contigs_index', shell=True)

            subprocess.check_call('bowtie2 -x index/contigs_index -f -k 1 --score-min L,0,-0.15 -t -p 4 -S contigs_alignment.sam karect_whole.fa', shell=True)
            subprocess.check_call('samtools view -F 4 -S contigs_alignment.sam > contigs_mapped.sam', shell=True)
            
            subprocess.check_call('python "%s"/tools/identify_misjoin_contigs.py Contigs.fa contigs_mapped.sam %s %s %s' %(base_path, args.read_len, args.fragment_len, args.fragment_std), shell=True)
            
            if args.threads:
                subprocess.check_call('sga index -t %s Contigs_clipped.fa' % args.threads, shell=True)
                subprocess.check_call('sga rmdup -t %s -e 0.005 Contigs_clipped.fa' % args.threads, shell=True)
            else:
                subprocess.check_call('sga index Contigs_clipped.fa', shell=True)
                subprocess.check_call('sga rmdup -e 0.005 Contigs_clipped.fa', shell=True)
            subprocess.check_call('mv Contigs_clipped.rmdup.fa Contigs_clipped.fa', shell=True)
        
        # ---------------- clean up -------------------------
        remove_files = "kept_pairs.txt kept_num.fa input_file.txt Plus_strand_reads.fa Minus_strand_reads.fa pair_end_connections.txt sequences.txt overlap_whole.txt Paths.txt Contigs.bwt Contigs_clipped.rbwt Contigs_clipped.rmdup.rbwt Contigs_clipped.rsai Contigs.rbwt Contigs.rmdup.rbwt Contigs.rsai Contigs_clipped.bwt Contigs_clipped.rmdup.bwt Contigs_clipped.rmdup.rsai Contigs_clipped.sai Contigs.rmdup.bwt Contigs.rmdup.rsai  Contigs.sai  Contigs_clipped.rmdup.dups.fa Contigs_clipped.rmdup.sai Contigs.rmdup.dups.fa Contigs.rmdup.sai samp.0.cnt samp.0.spm samp.des samp_edges.txt samp.esq samp.rlt samp.sds samp.ssp karect_whole.fa karect_whole_preprocessed.bwt karect_whole_preprocessed.rmdup.bwt karect_whole_preprocessed.rmdup.rbwt karect_whole_preprocessed.rsai karect_whole_preprocessed.fa karect_whole_preprocessed.rmdup.dups.fa karect_whole_preprocessed.rmdup.rsai karect_whole_preprocessed.sai karect_whole_preprocessed.rbwt karect_whole_preprocessed.rmdup.fa karect_whole_preprocessed.rmdup.sai res_graph_a.txt res_graph_b.txt contigs_alignment.sam contigs_mapped.sam PEG_nodes_sequences_before_removing.fa"
        remove_files = remove_files.split()
        remove_files.append('karect_'+fa1_file)
        remove_files.append('karect_'+fa2_file)
        remove_files.append('temp_res_'+fa1_file)
        remove_files.append('temp_res_'+fa2_file)
        for rev_file in remove_files:
            try:
                os.remove(rev_file)
            except OSError:
                pass
        try:
            subprocess.check_call('rm -r index', shell=True)
        except:
            pass

if __name__=='__main__':
    main()
