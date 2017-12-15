# PEHaplo

PEHaplo is a _de novo_ assembly tool for recovering virus haplotypes from virus quasispecies sequencing data. It utilizes overlap graph and paired-end information to recover virus haplotypes.

# Dependencies
PEHaplo is developed based on Python 2.7

Python module: [networkx](https://networkx.github.io)

[Karect](https://github.com/aminallam/karect)

[Readjoiner](http://www.zbh.uni-hamburg.de/forschung/gi/software/readjoiner.html)

[Apsp](https://github.com/chjiao/Apsp)

[SGA](https://github.com/jts/sga)

For contigs correction based on alignment:

[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

# Usage
```python 
python pehaplo.py [-h] -f1 INPUT_F1 -f2 INPUT_F2 -l OVERLAP_LEN -r READ_LEN [-F FRAGMENT_LEN] [-std FRAGMENT_STD] [-n DUP_N] [-correct CONTIG_CORRECT] [-t THREADS]
```

Arguments:

  -f1 INPUT_F1          input .1 part of paired-end fasta file
  
  -f2 INPUT_F2          input .2 part of paired-end fasta file
  
  -l OVERLAP_LEN, --overlap_len OVERLAP_LEN
                        overlap threshold between reads for overlap graph construction
                        
  -n DUP_N, --dup_n DUP_N
                        keep the reads that are duplicated at least n times
                        
  -r READ_LEN, --read_len READ_LEN
                        reads length
                        
  -F FRAGMENT_LEN, --fragment_len FRAGMENT_LEN
                        paired-end reads insert size, default as 2.5*read_len
                        
  -std FRAGMENT_STD     standard deviation of paired-end reads insert size,
                        default as 100
                        
  -n DUP_N, --dup_n DUP_N
                        the reads kept should be duplicated at least n times,
                        default as keep all the duplicates removed reads
                        
  -correct CONTIG_CORRECT
                        whether apply alignment based contigs
                        correction(yes/no), default: no
                        
  -t THREADS, --threads THREADS
                        threads for karect, sga, bowtie2
