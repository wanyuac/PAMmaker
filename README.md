# PAMmaker: creating allelic presence-absence matrices from results of gene detection

**Table of Contents**

- [Installation](#Installation)
    - [Dependencies](#Dependencies)
    - [Subdirectories of code](#Subdirectories)

- [A step-by-step guide for creating a PAM from SRST2 outputs](#guide_srst2)
    - [Targeted gene detection with SRST2](#srst2)
    - [Reliability assessment of allele calls](#uncertainty)
    - [Creating the allelic PAM](#make_PAM)

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

- [bash](https://www.gnu.org/software/bash/)

* [Rscript](https://www.r-project.org) (>=3.0)
* [Python](https://www.python.org/) (versions 2 and 3 compatible)

### Subdirectories of code<a name = "Subdirectories"/>

Subdirectory `utility` stores scripts used by pipeline `mk_allele_matrix.sh`. In addition, there are two subdirectories offering code for evaluating two characteristics of allele calls from SRST2, respectively:  

* `reliability`: scripts used for evaluating reliability of allele calls  
* `unicity`: scripts used for collecting evidence for identifying co-occurrence of alleles of the same gene in each sample  

<br/>

## 2. A step-by-step guide for creating a PAM from SRST2 outputs<a name = "guide_srst2"/>

This section demonstrates a typical procedure for processing SRST2-formatted results. In this demonstration, we create an allelic PAM and a genetic PAM from detected antimicrobial resistance genes (ARGs). Particularly, we use the SRST2-compatible ARG-ANNOT database version 2 ([ARGannot_r2.fasta](https://github.com/katholt/srst2/blob/master/data/ARGannot_r2.fasta)) as a reference database (hence the gene detection process to be launched is a targeted analysis).  



### 2.1. Targeted gene detection with SRST2<a name = "srst2"/>
The way for running SRST2 on a number of read sets varies between computer systems (see [SRST2 manual](https://github.com/katholt/srst2)). Assuming an [SLURM Workload Manager](https://slurm.schedmd.com/documentation.html) has been installed on a Linux cluster, users can use the following command line to detect ARGs in multiple genomes whose Illumina reads are accessible:

```
python srst2/scripts/slurm_srst2.py --script srst2/scripts/srst2.py --output demo --input_pe Reads/*_[1,2].fastq.gz --walltime '0-8:0:0' --threads 4 --memory 8192 --rundir ARGs --other_args "--gene_db srst2/data/ARGannot_r2.fasta --save_scores --report_all_consensus" > srst2_AMR.log
```

SRST2 creates an allele profile for each genome. Then using the command line below, users compile allele profiles into a single table, in which row names denote genomes, column names denote allele clusters (often known as genes) defined in the reference database, and each entry denotes an allele call.

```
python srst2/scripts/srst2.py --prev_output *_demo__genes__ARGannot_r2__results.txt --output Demo
```

Supposing we only analyse genomes of two samples, the output table or gene table may look like:

**Table 1**. Gene table of compiled allele calls from SRST2

| Sample  | floR\_Phe    | sul1\_Sul   |
|---------|-------------|-------------|
| strain1 | floR\_1212\*? | sul1\_1616\*? |
| strain2 | floR\_1212\*  | sul1\_1616?  |

This table shows three dubious allele calls (entries with a question mark). Now users may want to know whether allele calls `floR_1212*?`, `sul1_1616*?` and `sul1_1616?` indicate alleles that are actually present in these samples or merely partial hits to reference sequences. Section 2.2 addresses this question.  



### 2.2. Reliability assessment of allele calls<a name = "uncertainty"/>

To determine whether a dubious allele call can be accepted for further analysis, PAMmaker extracts summary statistics from score files (one per sample) produced by SRST2. Users may read a [blog post](https://www.microbialsystems.cn/en/post/srst2/) for details of SRST2's summary statistics that are used by PAMmaker for this assessment.

Supposing for Table 1 we can determine that `sul1_1616*?` and `sul1_1616?` are reliable while `floR_1212*?` remains suspicious, then the table of allele calls becomes

**Table 2**. Gene table of allele calls passed reliability filter

| Sample  | FloR\_Phe    | SulI\_Sul    |
|---------|-------------|-------------|
| strain1 | - | sul1\_1616\* |
| strain2 | floR\_1212\*  | sul1\_1616 |

Scripts in sub-directory `reliability` perform this reliability assessment in two steps:

- Retrieving and appending scores to allele calls in the compiled gene table (namely, top allele calls per gene/cluster per genome). Output file: `mergedScores__gene.scores`.

    ```python
    python ~/PAMmaker/reliability/collate_topAllele_scores.py --allele_calls ./Gene_calls/*_aEPEC__genes__ARGannot_r2__results.txt --allele_scores ./Scores/*_aEPEC__*.ARGannot_r2.scores --prefix mergedScores
    ```

    - Output file: `*.scores`, top-allele scores of clusters or genes.

- Filtering entries in the gene table by their quality scores. The output is a modified version of the input table. Of note, script `filter_allele_table.R` was known as `assessAlleleCallUncertainty.R`.

    ```R
    Rscript ~/PAMmaker/reliability/filter_allele_table.R --profiles aEPEC__compiledResults.txt --scores mergedScores__gene.scores --output aEPEC_srst2__reliableCalls
    ```

    - Two TSV-format output files: (1) `*.scores`, scores of top alleles passed the quality filter implemented in this script; (2) `*.txt`, the filtered gene table.

    

### 2.3. Creating the PAM<a name = "make_PAM" />

Shell script `mk_allele_matrix.sh` implements a pipeline for creating an allelic PAM from a gene table.

```bash
sh ./mk_allele_matrix.sh

"This pipeline modifies a combined gene table that is made by SRST2 to make an allele table in accordance with sequence clustering.
Dependencies: CD-HIT, Python 3 and R. Use this pipeline on a gene table after reliability assessment of its allele calls. In other
words, the input gene table should not have any question marks (laid by SRST2 to denote uncertain allele calls).

Arguments:
    -p: the program you want to run for sequence clustering, including its path when it is not in the current working directory
    -d: (optional) the directory where other scripts of this pipeline are located. Default: the current script directory.
    -o: (optional) the output directory name not followed by a forward slash. Default: clusters
    -a: arguments for CD-HIT-EST
    -f: input FASTA files for clustering (consensus allele sequences from SRST2)
    -g: a gene table from SRST2
    -c: an optional tab-delimited file of allele-hit scores. It is produced by assessAlleleCallUncertainty.R in the directory reliability_assessment.
    -s: stages to run. -s=all or -s=1,3,4 (maximum: 5) etc. Comma-delimited, no whitespace is allowed.
    -e: whether allele names contain unique identifiers. 0: no; 1 (default): yes. For example, set -e=1 when allele names look like dfrA12_1.
    -r: whether replace allele names with the representative allele name when their consensus sequences belong to the same cluster.

Usage:
    bash mk_allele_matrix.sh -p='apps/seq/cd-hit-est' -a='-c 1 -d 0 -s 1 -aL 1 -aS 1 -A 1 -uL 0 -uS 0 -p 1 -g 1' -f=*.fasta -g='profiles_res_genes.txt' -s=all -e=0
    bash mk_allele_matrix.sh -p='apps/seq/cd-hit' -d='../PAMmaker' -o='clrst' -a='-c 1 -d 0 -s 1 -aL 1 -aS 1 -A 1 -uL 0 -uS 0 -p 1 -g 1' -f=*.fasta -g='profiles_res_genes.txt' -s='1,3,4' -e=1

Warning:
    Sequence headers in original input FASTA files will be changed by this program, where whitespaces are replaced with '|'. So you may want to compress and backup your raw data before running this pipeline."
```

Example command line:

```bash
sh /vlsci/SG0006/shared/wan/scripts/cdhit_pip/mk_allele_matrix.sh -p="$HOME/software/cd-hit/bin/cd-hit-est" -a='-c 1 -d 0 -S 0 -AL 0 -AS 0 -U 0 -p 1 -g 1' -f='srst2_out/seq/*.fasta' -g='gene_table.txt' -s='all' -c='srst2__reliableCalls.scores'
```



This pipeline is comprised of four steps:

- Appending a sample name to every sequence name in each FASTA file. This step is carried out with the `sed` command of Linux.
- Clustering consensus sequences from SRST2 using `cd-hit-est`. The pipeline only supports clustering by perfect sequence identity at present.
- Tabulating the `*.clstr` file from `cd-hit-est` (script `./utility/tabulate_cdhit.py`), appending suffices to allele names as extended identifiers when multiple alleles are identified under an original allele name in the input gene table (script `./utility/clustering_allele_variants.py`, and creating a database of representative sequences named by new allele names (script `./utility/mk_allele_db.py`).
- Creating the allele matrix (script `./utility/convert_matrix.R`) and when the score file is provided, modifying allele names in the score file and produce summary statistics of `maxMAF` for each allele call in the score file.



Outputs:

- `allele_paMatrix.txt`, the allelic PAM;
- `modified_gene_table.txt`, gene table with allele names updated for cluster information (namely, the table comprised of extended allele names);
- `cluster_table.txt`, tabulated cluster information from the result of `cd-hit-est`;
- `allele_db.fna`, a FASTA file of alleles in the allelic PAM;
- `allele_name_replacement.txt`, a table linking original allele names to extended allele names based on sequence clustering.