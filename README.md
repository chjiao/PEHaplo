# PEHaplo
PEHaplo is a de novo assembly tool for recovering virus haplotypes from virus quasispecies sequencing data. It utilizes overlap graph and paired-end information to recover virus haplotypes.

# Dependencies
Python module: networkx
Karect
Readjoiner
Apsp
SGA

# Usage
preprocess.py -f1 INPUT_F1 -f2 INPUT_F2 -l OVERLAP_LEN [-n DUP_N]

Arguments:
  -f1 INPUT_F1          input .1 part of paired-end fasta file
  -f2 INPUT_F2          input .2 part of paired-end fasta file
  -l OVERLAP_LEN, --overlap_len OVERLAP_LEN
                        overlap threshold between reads for overlap graph construction
  -n DUP_N, --dup_n DUP_N
                        keep the reads that are duplicated at least n times
