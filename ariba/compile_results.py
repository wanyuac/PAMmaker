#!/usr/bin/env python

"""
This script compiles ARIBA's outputs for individual samples. Specifically, it performs two tasks:
    1. Create a table of identified alleles for all samples with columns: sample, cluster, allele, exact match/variant
       of the reference allele, percent of nucleotide identity, coverage of the assembled allele to its reference.
    2. Pool FASTA files of ARIBA's output allele sequences into one file and append sample names to sequence IDs (in
       the same format as that for SRST2's output consensus sequences).

Notes:
    1. Sample names from filenames of input reports and FASTA files must match, or the script will terminate.
    2. Dependencies: packages pandas and Bio (BioPython).

Example command:
    python compile_results.py --in_fastas gene/*_genes.fna --in_reports report/*_report.tsv --out_fasta alleles.fna \
        --out_table alleles.tsv --ext_fastas '_genes.fna' --ext_reports '_report.tsv'

Copyright (C) 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>
Publication: 11 Nov 2020; the latest modification: 29 Dec 2020
Previous name: pool_seqs.py
"""

import os
import sys
import pandas
from argparse import ArgumentParser
from collections import namedtuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_arguments():
    parser = ArgumentParser(description = "Pool ARIBA's output allele sequences into one FASTA file and append sample names to sequence IDs",\
             epilog = "This is a helper script of R package GeneMates.")
    parser.add_argument("--in_fastas", nargs = "+", dest = "in_fastas", type = str, required = True, help = "Input FASTA files")
    parser.add_argument("--in_reports", nargs = "+", dest = "in_reports", type = str, required = True, help = "Input report files")
    parser.add_argument("--out_fasta", dest = "out_fasta", type = str, required = False, default = "alleles.fna", help = "Output FASTA file of pooled allele sequences")
    parser.add_argument("--out_table", dest = "out_table", type = str, required = False, default = "alleles.tsv", help = "Output table about identified alleles")
    parser.add_argument("--ext_fastas", dest = "ext_fastas", type = str, required = False, default = "_genes.fna", help = "Filename extension of input FASTA files to be removed for sample names")
    parser.add_argument("--ext_reports", dest = "ext_reports", type = str, required = False, default = "_report.tsv", help = "Filename extension of input report files to be removed for sample names")
    return parser.parse_args()

def main():
    args = parse_arguments()
    sample_count = 0
    samples = pair_input_files(args.in_fastas, args.in_reports, args.ext_fastas, args.ext_reports)  # Returns a dictionary {sample : namedtuple(fasta, report)}
    fasta_out = open(args.out_fasta, "w")
    table_out = open(args.out_table, "w")
    table_out.write("\t".join(["Sample", "Reference", "Cluster", "Ref_len", "Ref_assembled", "Contig", "Identity", "Coverage"]) + "\n")  # Print the table header
    for s, fs in samples.items():  # Go through result files of each sample
        report = import_report(fs.report)  # A report comprised of selected columns and non-redundant rows
        fasta = import_fasta(fs.fasta)  # A simpler dictionary (id : seq) than that created by SeqIO.to_dict (seq.id : SeqRecord)
        for _, row in report.iterrows():  # Iterate through rows of the data frame
            ctg = row["ctg"]  # ID used for searching against keys of dictionary 'fasta'
            cov_perc = round(row["ref_base_assembled"] / row["ref_len"], 4) * 100  # Coverage (%) of the reference allele
            is_var = 0 if (row["pc_ident"] == 100 and cov_perc == 100) else 1  # Whether this hit is an exact match to its reference or not (variant)
            if ctg in fasta.keys():
                fasta_out.write(">%s|%i|%s %s\n%s\n" % (row["ref_name"], is_var, s, ctg, fasta[ctg]))
                table_out.write("\t".join([s] + list(map(str, list(row) + [cov_perc]))) + "\n")  # The map function applies str() to every element of a list.
            else:
                """
                It is normal that an allele from the record file is absent from the FASTA file due to insufficient nucleotide
                identity or coverage following cut-offs set for ARIBA runs.
                """
                print("Notice: allele ID %s is not found in %s. Skipped this allele." % (ctg, s + args.ext_fastas), file = sys.stderr)
        sample_count += 1
    fasta_out.close()
    table_out.close()
    print(str(sample_count) + " samples have been successfully processed.", file = sys.stdout)
    return

def pair_input_files(input_fastas, input_reports, ext_fastas, ext_reports):
    """
    Extract sample names from filenames of input FASTA files, assuming these sample names are also present in filenames of input reports.
    """
    samples = dict()  # Outcome of this function
    Sample = namedtuple("Sample", ["fasta", "report"])
    samples_fastas = get_sample_names(input_fastas, ext_fastas)
    samples_reports = get_sample_names(input_reports, ext_reports)
    for s, f in samples_fastas.items():
        if s in samples_reports.keys():
            samples[s] = Sample(fasta = f, report = samples_reports[s])
        else:
            print("Warning: skipped sample " + s + " due to absence of its report file.", file = sys.stderr)
    return samples

def get_sample_names(input_files, ext):
    """ A subordinate function of pair_input_files """
    fs = dict()
    for f in input_files:
        f_base = os.path.basename(f)
        sample = f_base.replace(ext, "")  # No change applies if ext is not found in the filename
        fs[sample] = f
    return fs

def import_report(report_file):
    """
    Reads and processes ARIBA's report file for a given sample
    Data type of column pc_ident: numpy.float64
    """
    report = pandas.read_csv(report_file, sep = "\t", index_col = None)  # Returns a data frame
    report = report[["ref_name", "cluster", "ref_len", "ref_base_assembled", "ctg", "pc_ident"]]  # Select columns using a list of column names
    report = report.drop_duplicates()  # Remove duplicated rows as ARIBA prints each mutation on a separate row
    return report

def import_fasta(fasta_file):
    """ Reads ARIBA's output gene file for a given sample and parse sequence headers """
    fasta = dict()
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        seqid = seq_record.id
        seqid_fields = seqid.split(".")  # First element: cluster name; 1st - 5th fields make up the contig name.
        allele_id = ".".join(seqid_fields[0 : 5])  # This ID can handle scenarios where multiple alleles of the same cluster are found, although these scenarios are rare.
        fasta[allele_id] = str(seq_record.seq)  # Usual scenario
    return fasta

if __name__ == "__main__":
    main()
