# PAMmaker: pipelines creating allelic presence-absence matrices from results of gene screen

This repository consists of Python and R scripts that produce an allelic presence-absence matrix (PAM) from the outcome of a gene screen on bacterial samples (strains/isolates). Assuming _m_ alleles are identified in _n_ samples, the allelic PAM _**A** = (a<sub>ij</sub>)_ is an n-by-m binary matrix, where _a<sub>ij</sub> = 1_ when the _j_-th allele is present in the _i_-th sample, and _a<sub>ij</sub> = 0_ otherwise. This tool assumes the gene screen process produces files in the output format of [SRST2](github.com/katholt/srst2) and [screen\_genes\_in\_assemblies](github.com/wanyuac/screen_genes_in_assemblies).

## Installation
```
git clone https://github.com/wanyuac/PAMmaker.git
```

### Dependencies

* [R](https://www.r-project.org) (>=3.0)
* Python versions 2 and 3 compatible
* Linux bash

### Sub-directories

* reliability\_assessment: scripts used for evaluating reliability of allele calls
* unicity\_assessment: scripts used for collecting evidence for identifying co-occurrence of alleles of the same gene in each sample

## Usage
### Gene screen
Assuming a user uses SRST2 for screening antimicrobial resistance genes (ARGs) in bacterial samples, then an example SRST2 command line is (see SRST2 manual for details):

```
python srst2/scripts/slurm_srst2.py --script srst2/scripts/srst2.py --output demo --input_pe Reads/*_[1,2].fastq.gz --walltime '0-8:0:0' --threads 4 --memory 8192 --rundir ARGs --other_args "--gene_db srst2/data/ARGannot_r2.fasta --save_scores --report_all_consensus" > srst2_AMR.log
```




