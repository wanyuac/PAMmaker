#!/usr/bin/env python

"""
Making a database of representative sequences of alleles. It takes as two main inputs (1) -c: a list of FASTA files
consisting of consensus sequences per strain. A strain name must be attached to the header through a "|" character:
>[allele name and relevant information]|[strain name]. For example, >256__SulII_Sul__SulII__1219.variant|234-12.
(2) -a: a tab-delimited table of allele profiles, where extended allele IDs are incorporated.

Notice since this script works on representative sequences identified using CD-HIT-EST, there is no redundancy in
both input and output database.

Dependency: Python versions 2 and 3 compatible, BioPython
Usage
    python mk_allele_db.py -a modified_allele_matrix.txt -c *_all_consensus.fasta -o allele_db.fna
    python mk_allele_db.py -a modified_allele_matrix.txt -r allele_name_replacement.txt -c *_all_consensus.fasta -o allele_db.fna

Copyright 2017 Yu Wan <wanyuac@gmail.com>
Licensed under the Apache License, Version 2.0
Created on 28 Sep 2017. Lastest edition: 7 Aug 2018.
"""

from __future__ import print_function
from argparse import ArgumentParser
from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from os import path
import sys


def parse_arguments():
    parser = ArgumentParser(description = "Making a database of representative alleles produced by CD-HIT-EST")
    parser.add_argument("-a", type = str, required = True, help = "Allele profiles or presence/absence matrix, with extended allele IDs attached to allele names")
    parser.add_argument("-p", action = "store_true", help = "Flag it when a presence/absence matrix is provided instead of allele profiles")
    parser.add_argument("-r", type = str, required = False, default = "", help = "A tab-delimited table produced by clustering_allele_variants.py to replace allele names")
    parser.add_argument("-c", nargs = "+", type = str, required = True, help = "A list of FASTA files for consensus sequences")
    parser.add_argument("-o", type = str, required = False, default = "allele_db.fna", help = "Output file (name and path)")
    parser.add_argument("-k", action = "store_true", help = "Set it to keep the name of the strain from which the allele sequence is drawn")
    parser.add_argument("-n", action = "store_true", help = "Flag it when allele IDs are not assigned based on perfect sequence identity")
    return parser.parse_args()


def import_allele_table(tab_file):
    """
    Read and parse allelic profiles (with CD-HIT-EST cluster IDs attached to allele names)
    allele_tab = {allele1 : [strain11, ..., strain1N], ..., alleleM : [strainM1, ..., strainMN]}
    """
    allele_tab = {}
    strains = []
    with open(tab_file, "rU") as f:
        lines = f.read().splitlines()  # read all lines once
    header = lines[0]
    content = lines[1 : ]
    genes = header.split("\t")[1 : ]  # drop the first entry "Sample", leaving allele cluster (gene) names
    gene_num = len(genes)
    for line in content:  # one line per strain, no duplicates
        entries = line.split("\t")
        strain = entries[0]
        strains.append(strain)
        entries = entries[1 : ]  # remove the first entry
        for i in range(0, gene_num):
            allele = entries[i]  # allele name
            if allele != "-":
                if allele in list(allele_tab.keys()):
                    allele_tab[allele].append(strain)
                else:
                    allele_tab[allele] = [strain]
    return [allele_tab, strains]  # a mixed list comprised of a dictionary and a list


def import_PAM(pam_file):
    """
    This function is added to read allele information from a presence/absence matrix (PAM) instead of from allele profiles.
    Here, I assume that in every PAM, one represents the presence whereas zero represents the absence of an allele in a strain.
    """
    allele_tab = {}
    strains = []
    with open(pam_file, "rU") as f:
        lines = f.read().splitlines()  # read all lines once
    header = lines[0]
    content = lines[1 : ]
    alleles = header.split("\t")[1 : ]  # drop the first entry "Sample", leaving allele cluster (gene) names
    allele_num = len(alleles)
    for line in content:
        entries = line.split("\t")
        strain = entries[0]
        strains.append(strain)
        entries = entries[1 : ]  # remove the first entry
        for i in range(0, allele_num):
            pa = entries[i]  # presence/absence status of the (i+1)-th allele
            if pa == "1":
                allele = alleles[i]
                if allele in list(allele_tab.keys()):  # when other strains have the same allele
                    allele_tab[allele].append(strain)
                else:
                    allele_tab[allele] = [strain]
    return [allele_tab, strains]


def get_allele_name(id):
    """
    This is a subordinate function of import_consensus_seqs.
    For instance, it returns "RmtB_1580" from "309__RmtB_AGly__RmtB__1580.consensus|ERR560522".
    """
    id = id.split(".")[0]
    _, _, allele, seq_id = id.split("__")
    return allele + "_" + seq_id  # e.g. RmtB__1580 -> RmtB_1580
    

def import_consensus_seqs(fastas, strains_accept = None):
    """
    Return value: {allele : {strain : sequence}}. It can be a huge dictionary.
    strains_accept: a list of strains to be accepted for the output
    """
    cons_seqs = defaultdict(dict)
    for fasta in fastas:
        records = list(SeqIO.parse(fasta, "fasta"))       
        
        # Process the first record to determine if the current strain should be processed
        record = records[0]
        strain = record.id.split("|")[1]
        if strains_accept != None:
            if strain not in strains_accept:    
                continue  # stop the current "for" iteration and read the next FASTA file
        allele = get_allele_name(record.id)
        cons_seqs[allele][strain] = record.seq  # A strain must have a single allele per gene. So there is no append method required.
        
        # Process other records of the current strain when the "continue" command is not called
        if len(records) > 1:  # when there are other sequences in the current FASTA file
            for record in records[1 : ]:
                strain = record.id.split("|")[1]
                allele = get_allele_name(record.id)
                cons_seqs[allele][strain] = record.seq
    return cons_seqs


def import_mapping_table(tab_file):
    """
    Import the mapping table of allele names when it is specified
    output data structure: {new_allele_name : {strain : previous_allele_name, ...}}
    """
    with open(tab_file, "rU") as f:
        lines = f.read().splitlines()
    lines = lines[1:]  # omit the header line
    mapping = defaultdict(dict)
    for line in lines:
        strain, prev_name, _, new_name = line.split("\t")
        mapping[new_name][strain] = prev_name
    
    return mapping


def retrieve_prev_allele(allele, strain, mapping):
    """
    Find out the allele name in a FASTA file. The new allele name does not contain
    any extended allele ID.
    """
    if allele in list(mapping.keys()):
        if strain in list(mapping[allele].keys()):
            prev_name = mapping[allele][strain]
        else:
            prev_name = allele
    else:
        prev_name = allele
    
    return prev_name


def main():
    args = parse_arguments()
    
    # Sanity check
    replace_names = args.r != ""
    if replace_names:
        files = [args.a, args.r] + list(args.c)
    else:
        files = [args.a] + list(args.c)
    for f in files:  # args.c is a list when len(args.c) > 1, otherwise, it is a string. So use list() for safty.
        if not path.isfile(f):
            sys.exit("Error: file " + f + " is not accessible.")
    
    # Get the mapping table for allele name substitution
    if replace_names:
        mapping = import_mapping_table(args.r)
    
    # The following two import functions produce the same data structure.
    if args.p:
        allele_tab, strains = import_PAM(args.a)  # read allele assignment from an allelic PAM
    else:
        allele_tab, strains = import_allele_table(args.a)  # import from allele profiles otherwise
    seqs = import_consensus_seqs(args.c, strains)  # establish an internal database of consensus sequences
    
    # Construct the allele database
    out = open(args.o, "w")
    for allele in list(allele_tab.keys()):  # for each allele in the allele table, with the extended allele ID
        if args.n:  # "not clustered based on perfect matches"
            """
            Sometimes a user may not assign extended allele IDs based on 100% sequence identity. For example, consensus
            sequences may be clustered under a nucleotide identity of 95% and allows a few difference in the sequence length.
            In this scenario, the user need to keep all consensus sequences even they have been assigned to the same cluster
            (i.e. have the same extended allele ID).
            """
            strains = allele_tab[allele]
            for strain in strains:
                allele_name_root = allele.split(".")[0]
                if replace_names:
                    prev_allele = retrieve_prev_allele(allele_name_root, strain, mapping)
                    seq = seqs[prev_allele][strain]
                else:
                    seq = seqs[allele_name_root][strain]
                print(">" + allele + " " + strain, file = out)
                print(seq, file = out)
        else:  # Sequences are clustered based on 100% nucleotide identity and length.
            strain = allele_tab[allele][0]  # take the first strain to retrive the sequence as all other strains have the same sequence
            """
            The next command uses the split method to drop the extended allele ID from the complete allele ID as there is not
            such an ID in the FASTA file.
            """
            allele_name_root = allele.split(".")[0]
            if replace_names:
                prev_allele = retrieve_prev_allele(allele_name_root, strain, mapping)
                seq = seqs[prev_allele][strain]
            else:
                seq = seqs[allele_name_root][strain]
        
            # Print the sequence header first
            if args.k:
                print(">" + allele + " " + strain, file = out)
            else:
                print(">" + allele, file = out)  # allele: the new name
            
            # Then print the sequence itself to complete the current record
            print(seq, file = out)
    out.close()


if __name__ == "__main__":  # if this script is run as the main script. Not "main", or you will not get any output.
    main()
    