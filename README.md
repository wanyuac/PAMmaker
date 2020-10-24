# PAMmaker: pipelines creating allelic presence-absence matrices from SRST2 results

**Table of Contents**

1. [Installation](#Installation)
    - [Dependencies](#Dependencies)
    - [Subdirectories of code](#Subdirectories)
2. [A step-by-step guide for creating PAMs](#guide)
    - [Running SRST2 for targeted gene detection](#srst2)
    - [Uncertainty assessment of allele calls](#uncertainty)

<br/>

This repository consists of Python and R scripts that produce an allelic presence-absence matrix (PAM) from the outcome of SRST2-based gene detection for bacterial samples (strains or isolates). Assuming _m_ alleles are identified in _n_ samples, the allelic PAM _**A** = (a<sub>ij</sub>)_ is an n-by-m binary matrix, where _a<sub>ij</sub> = 1_ when the _j_-th allele is present in the _i_-th sample, and _a<sub>ij</sub> = 0_ otherwise.

PAMmaker assumes that the gene detection process produces files in the output format of [SRST2](https://github.com/katholt/srst2). PAMmaker is a helper tool of R package [GeneMates](https://github.com/wanyuac/GeneMates) and it implements two pipelines for allele calls:  

* Reliability assessment, which uses allele-call scores to determine whether an allele call is reliable or not;  
* Unicity assessment, which summarises minor-allele frequencies (MAFs) of each allele call across samples to show whether it has homologous alleles in some samples.  

Notice PAMmaker is not applicable to results of [geneDetector](https://github.com/wanyuac/geneDetector) (gene detection for genome assemblies), because every alignment of an assembly against a reference genome does not generate dubious allele calls or estimates of MAFs.

<br/>

## 1. Installation<a name = "Installation"/>
```
git clone https://github.com/wanyuac/PAMmaker.git
```

### Dependencies<a name = "Dependencies"/>

* [R](https://www.r-project.org) (>=3.0)
* Python (versions 2 and 3 compatible)
* Linux bash

### Subdirectories of code<a name = "Subdirectories"/>
There are two subdirectories under PAMmaker:  

* reliability\_assessment: scripts used for evaluating reliability of allele calls  
* unicity\_assessment: scripts used for collecting evidence for identifying co-occurrence of alleles of the same gene in each sample  

<br/>

## 2. A step-by-step guide for creating PAMs<a name = "guide"/>

This section demonstrates a typical procedure for processing SRST2-formatted results. In this demonstration, we create an allelic PAM and a genetic PAM from detected antimicrobial resistance genes (ARGs). Particularly, we use the SRST2-compatible ARG-ANNOT database version 2 ([ARGannot_r2.fasta](https://github.com/katholt/srst2/blob/master/data/ARGannot_r2.fasta)) as a reference database (hence the gene detection process to be launched is a targeted analysis).

### 2.1. Running SRST2 for targeted gene detection<a name = "srst2"/>
Assuming that an [SLURM Workload Manager](https://slurm.schedmd.com/documentation.html) has been installed on a Linux cluster, users can run SRST2 using the following command line to detect ARGs in multiple samples for which Illumina reads are accessible (see [SRST2 manual](https://github.com/katholt/srst2)):

```
python srst2/scripts/slurm_srst2.py --script srst2/scripts/srst2.py --output demo --input_pe Reads/*_[1,2].fastq.gz --walltime '0-8:0:0' --threads 4 --memory 8192 --rundir ARGs --other_args "--gene_db srst2/data/ARGannot_r2.fasta --save_scores --report_all_consensus" > srst2_AMR.log
```

Eventually, the user compiles allele profiles from individual samples into a single table, in which row names denote samples and column names are gene names. Each entry of the table is an allele call to a gene in the reference database. The command line for compiling individual results is as follows:

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

