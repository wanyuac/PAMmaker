# PAMmaker: creating allelic presence-absence matrices from results of gene detection

**Table of Contents**

- [Installation](#Installation)
    - [Dependencies](#Dependencies)
    - [Subdirectories of code](#Subdirectories)

- [A step-by-step guide for creating a PAM from SRST2 outputs](#guide_srst2)
    - [Running SRST2 for targeted gene detection](#srst2)
    - [Uncertainty assessment of allele calls](#uncertainty)

<br/>

This repository consists of Python and R scripts that produce an allelic presence-absence matrix (PAM) from [SRST2](https://github.com/katholt/srst2)-compatible outputs. Specifically, assuming _m_ alleles are identified in _n_ samples (strains or isolates), the allelic PAM _**A** = (a<sub>ij</sub>)_ is an n-by-m binary matrix, where _a<sub>ij</sub> = 1_ when the _j_-th allele is present in the _i_-th sample, and _a<sub>ij</sub> = 0_ otherwise.

PAMmaker is a helper tool of R package [GeneMates](https://github.com/wanyuac/GeneMates) for creating mandatory inputs of main functions `findPhysLink` and `lmm`. Moreover, it offers two pipelines for evaluating quality of SRST2's allele calls:  

* Reliability assessment, which uses allele-call scores to determine whether an allele call is reliable;  
* Unicity assessment, which summarises minor-allele frequencies (MAFs) of each allele call across samples to show whether it has homologous alleles in some samples.  

**Note**

PAMmaker is not applicable to results of [geneDetector](https://github.com/wanyuac/geneDetector) (gene detection for genome assemblies), because the alignment of reference alleles (query sequences) against contigs (subject sequences) does not generate dubious allele calls or estimates of MAFs.

**Citation**

Wan, Y., Wick, R.R., Zobel, J., Ingle, D.J., Inouye, M., Holt, K.E. GeneMates: an R package for detecting horizontal gene co-transfer between bacteria using gene-gene  associations controlled for population structure. *BMC Genomics* 21, 658 (2020). https://doi.org/10.1186/s12864-020-07019-6.

<br/>

## 1. Installation<a name = "Installation"/>
```
git clone https://github.com/wanyuac/PAMmaker.git
```

### Dependencies<a name = "Dependencies"/>

PAMmaker requires three code interpreters and a Linux-compatible operating system:

- [GNU bash](https://www.gnu.org/software/bash/)

* [R](https://www.r-project.org) (>=3.0)
* [Python](https://www.python.org/) (versions 2 and 3 compatible)

### Subdirectories of code<a name = "Subdirectories"/>
There are two subdirectories offering code for evaluating two characteristics of allele calls from SRST2, respectively:  

* `reliability`: scripts used for evaluating reliability of allele calls  
* `unicity`: scripts used for collecting evidence for identifying co-occurrence of alleles of the same gene in each sample  

<br/>

## 2. A step-by-step guide for creating a PAM from SRST2 outputs<a name = "guide_srst2"/>

This section demonstrates a typical procedure for processing SRST2-formatted results. In this demonstration, we create an allelic PAM and a genetic PAM from detected antimicrobial resistance genes (ARGs). Particularly, we use the SRST2-compatible ARG-ANNOT database version 2 ([ARGannot_r2.fasta](https://github.com/katholt/srst2/blob/master/data/ARGannot_r2.fasta)) as a reference database (hence the gene detection process to be launched is a targeted analysis).

### 2.1. Targeted gene detection with SRST2<a name = "srst2"/>
Assuming that an [SLURM Workload Manager](https://slurm.schedmd.com/documentation.html) has been installed on a Linux cluster, users can run SRST2 using the following command line to detect ARGs in multiple samples for which Illumina reads are accessible (see [SRST2 manual](https://github.com/katholt/srst2)):

```
python srst2/scripts/slurm_srst2.py --script srst2/scripts/srst2.py --output demo --input_pe Reads/*_[1,2].fastq.gz --walltime '0-8:0:0' --threads 4 --memory 8192 --rundir ARGs --other_args "--gene_db srst2/data/ARGannot_r2.fasta --save_scores --report_all_consensus" > srst2_AMR.log
```

Eventually, as demonstrated below, the user compiles allele profiles from individual samples into a single table, in which row names denote samples and column names are gene names. Each entry of the table is an allele call to a gene in the reference database. The command line for compiling individual results is as follows:

```
python srst2/scripts/srst2.py --prev_output *_demo__genes__ARGannot_r2__results.txt --output Demo
```

In this demonstration, we only analyse two samples and assume their compiled allele profile is 

| Sample  | FloR\_Phe    | SulI\_Sul    |
|---------|-------------|-------------|
| strain1 | FloR\_1212\*? | SulI\_1616\*? |
| strain2 | FloR\_1212\*  | SulI\_1616?  |

Note that this table shows three dubious allele calls (that is, entries with a question mark). Now our question is, whether allele calls FloR\_1212*?, SulI\_1616\*? and SulI\_1616? indicate alleles that are actually present in our samples or merely partial hits to reference sequences?

### 2.2. Uncertainty assessment of allele calls<a name = "uncertainty"/>

PAMmaker extracts summary statistics from score files produced by SRST2 to determine whether a dubious allele call (such as "FloR\_1212*?" and "SulI\_1616?") can be treated. Readers may see a [blog post](https://www.microbialsystems.cn/en/post/srst2/) for a detailed explanation of SRST2's summary statistics that are used for determining allele calls.

Supposing we can determine that SulI\_1616\*? and SulI\_1616? are reliable while FloR\_1212\*? is unreliable, then the gene profile becomes

| Sample  | FloR\_Phe    | SulI\_Sul    |
|---------|-------------|-------------|
| strain1 | - | SulI\_1616\* |
| strain2 | FloR\_1212\*  | SulI\_1616  |

We can use scripts under the sub-directory reliability\_assessment of PAMmaker to perform this uncertainty assessment. The procedure is comprised of two steps.

```
python ~/PAMmaker/reliability_assessment/collate_topAllele_scores.py --allele_calls ./Gene_calls/*_aEPEC__genes__ARGannot_r2__results.txt --allele_scores ./Scores/*_aEPEC__*.ARGannot_r2.scores --prefix mergedScores
```

```
Rscript ~/PAMmaker/reliability_assessment/assessAlleleCallUncertainty.R --profiles aEPEC__compiledResults.txt --scores mergedScores__gene.scores --output aEPEC_srst2__reliableCalls
```

