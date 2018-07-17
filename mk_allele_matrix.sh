#!/bin/bash
# Author: Yu Wan <wanyuac@gmail.com>
# Copyright 2017 Yu Wan
# Licensed under the Apache License, Version 2.0
# Creation: 5 Sep 2016, latest edition: 2 Feb 2018

display_usage() {
echo "
This pipeline modifies a combined gene table that is made by SRST2 to make an allele table in accordance with sequence clustering.
Dependencies: CD-HIT, Python 3 and R. Use this pipeline on a gene table after reliability assessment of its allele calls. In other
words, the input gene table should not have any question marks (laid by SRST2 to denote uncertain allele calls).
Arguments
    -p: the program you want to run for sequence clustering, including its path when it is not in the current working directory
	-d: (optional) the directory where other scripts of this pipeline are located. Default: the current script directory.
	-o: (optional) the output directory name not followed by a forward slash. Default: clusters
    -a: arguments for CD-HIT
    -f: input FASTA files for clustering
	-g: a gene table from SRST2
	-c: an optional tab-delimited file of allele-hit scores. It is produced by assessAlleleCallUncertainty.R in the directory reliability_assessment.
    -s: stages to run. -s=all or -s="1,3,4" (maximum: 5) etc. Comma-delimited, no whitespace is allowed.
    -e: whether allele names contain unique identifiers. 0: no; 1 (default): yes. For example, set -e=1 when allele names look like dfrA12_1.
Usage:
    bash mk_allele_matrix.sh -p='apps/seq/cd-hit-est' -a='-c 1 -d 0 -s 1 -aL 1 -aS 1 -A 1 -uL 0 -uS 0 -p 1 -g 1' -f=*.fasta -g='profiles_res_genes.txt' -s=all -e=0
    bash mk_allele_matrix.sh -p='apps/seq/cd-hit' -d='../PAMmaker' -o='clrst' -a='-c 1 -d 0 -s 1 -aL 1 -aS 1 -A 1 -uL 0 -uS 0 -p 1 -g 1' -f=*.fasta -g='profiles_res_genes.txt' -s='1,3,4' -e=1
Warning:
    Sequence headers in original input FASTA files will be changed by this program, where whitespaces are replaced with '|'.
    So you may want to compress and backup your raw data before running this pipeline.
"
}

# Display usage information if a null argument, "-h" or "--help" is encountered.
if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	display_usage
	exit
fi

# Defaults
code_path=$(dirname "$0")  # get the directory name of this script
out_path="clusters"  # default name of the output directory under $PWD
id_ext=1

# Read arguments from the command line ##########
for i in "$@"; do  # loop through every argument
    case $i in
        -p=*)
        program="${i#*=}"
        ;;
		-d=*)
		code_path="${i#*=}"
		;;
		-o=*)
		out_path="${i#*=}"
		;;
        -a=*)
        arguments="${i#*=}"  # get the string after the equal sign in $i
        ;;
        -f=*)
        fasta="${i#*=}"  # get a list of FASTA files
        ;;
        -g=*)
        allele_table="${i#*=}"  # an allele table from SRST2 after reliability assessment
        ;;
		-c=*)
		score_file="${i#*=}"  # an optional tab-delimited file of allele-hit scores
		;;
        -s=*)
        stages="${i#*=}"  # get stages (comma-delimited, no whitespaces) to run
        ;;
        -e=*)
        id_ext="${i#*=}"  # "0" or "1"
		;;
    esac
done

# Parse stages ##########
if [[ "$stages" -eq "all" ]]; then
    stages=(1, 2, 3, 4, 5)
else
    IFS="," read -r -a stages <<< $stages  # split 
fi

# make an output directory if it does not exist
if [ ! -d ${out_path} ]; then
	mkdir $out_path
fi

# 1. Append a sample name to every sequence name in each FASTA file ##########
if [[ "${stages[@]}" =~ "1" ]]; then  # if the array contains 1
	echo "Appending sample names to sequences."
    for fna in ${fasta}; do
        sed -i 's/ /\|/g' ${fna}  # replace whitespaces with '|' in-place because CD-HIT does not read information beyond the first space
    done
fi

# 2. Cluster consensus sequences ##########
if [[ "${stages[@]}" =~ "2" ]]; then
	echo "Clustering nucleotide sequences with CD-HIT-EST."
    cat $fasta > fasta.tmp
    ${program} -i fasta.tmp -o ${out_path}/clusters.fna $arguments > ${out_path}/$(basename "${program}").log  # run either cd-hit or cd-hit-est
    rm fasta.tmp
fi

# 3. Append suffices to allele names ##########
if [[ "${stages[@]}" =~ "3" ]]; then
    if [ ! -f "${out_path}/cluster_table.txt" ]; then
        echo "Tabulating CD-HIT-EST outputs."
        # As it is better to keep the cluster_table.txt for additional analyses, I do not make a longer pipe here.
        python3 ${code_path}/tabulate_cdhit.py ${out_path}/clusters.fna.clstr > ${out_path}/cluster_table.txt  # run the next command as long as this command successed
    else
        echo "Skip tabulating CD-HIT-EST outputs because it has been done already. Delete the file if you want to rerun this step."
	fi
    
	echo "Assigning allele identifiers into the genotype table."
    if [[ ${id_ext} = 1 ]]; then  # Allele names contain the unique identifiers.
        cat ${out_path}/cluster_table.txt | python3 ${code_path}/clustering_allele_variants.py $allele_table > ${out_path}/modified_allele_matrix.txt
    else
        echo "No unique allele identifiers are present."
        cat ${out_path}/cluster_table.txt | python3 ${code_path}/clustering_allele_variants.py $allele_table 1 > ${out_path}/modified_allele_matrix.txt
    fi
	
	AM_present=`[ -f "${out_path}/modified_allele_matrix.txt" ]`
	
	if [ ! -f allele_db.fna ] && ${AM_present}; then
		echo "Creating a database of representative sequences named by allele names."
		python3 ${code_path}/mk_allele_db.py -a ${out_path}/modified_allele_matrix.txt -c ${fasta} -o ${out_path}/allele_db.fna
	fi
fi

if [[ "${stages[@]}" =~ "4" ]] && ${AM_present}; then
	echo "Converting the genotype table into an allelic presence/absence matrix."
    Rscript ${code_path}/convert_matrix.R ${out_path}/modified_allele_matrix.txt ${out_path}
fi

# Modify allele names in the score file and produce summary statistics of maxMAF for each allele call
# Please ensure the score file is accessible.
# This step does not require every allele in the score file to be present in the allelic PAM. Only those in the PAM
# will be retained in the new score file.
aPAM="${out_path}/allele_paMatrix.txt"  # output of convert_matrix.R
if [[ "${stages[@]}" =~ "5" ]] && [ ! -z "$score_file" ] && [ -f "${aPAM}" ]; then
	echo "Summarising maximal minor-allele-frequencies (MAF) for each allele call."
	Rscript ${code_path}/unicity_assessment/maxMAFstats.R "${out_path}/allele_paMatrix.txt" "${score_file}" "${out_path}/allele"
fi
