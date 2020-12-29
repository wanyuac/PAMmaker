#!/usr/bin/env python

"""
Converts output TSV of utility/tabulate_cdhit.py to an allelic presence-absence matrix (PAM) and saves
it in TSV format.

Command demonstration:
    python clusters2pam.py -i alleles.clstr.tsv -om alleles_pam.tsv -ot alleles_clstr_updated.tsv

Dependency: module pandas, Python v3

Copyright (C) 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Publication: 11 Nov 2020; the latest modification: 29 Dec 2020
"""

import os
import sys
import pandas
from argparse import ArgumentParser

def parse_arguments():
    parser = ArgumentParser(description = "Creating an allelic presence-absence matrix from a table of sequence clusters",\
             epilog = "This is a helper script of R package GeneMates.")
    parser.add_argument("-i", dest = "i", type = str, required = True, help = "Input tab-delimited table of clusters from CD-HIT-EST")
    parser.add_argument("-om", dest = "om", type = str, required = False, default = "allelic_PAM.tsv", help = "Output presence-absence matrix in TSV format")
    parser.add_argument("-ot", dest = "ot", type = str, required = False, default = "clusters.tsv", help = "Output table about sequence clusters")
    return parser.parse_args()

def main():
    args = parse_arguments()

    # Import the input TSV file as a data frame
    input_tsv = args.i
    if os.path.exists(input_tsv):
        clusters = parse_seqid(pandas.read_csv(input_tsv, sep = "\t", index_col = None))  # Import the table as a data frame of six columns and parse column 'seqid'
    else:
        print("Error: the input file is not accessible.", file = sys.stderr)
        sys.exit(1)
    
    # Assign allele IDs and create an allelic PAM from data frame "clusters"
    pam = create_allelic_pam(clusters)
    pam.to_csv(args.om, header = True, index = False, sep = "\t", float_format = None)
    clusters.to_csv(args.ot, header = True, index = False, sep = "\t", float_format = None)
    return

def parse_seqid(df):
    """ 
    Parses column 'seqid' into two new columns and returns a data frame of eight columns.
    """
    seqids = df["seqid"].tolist()  # Convert a column into a list
    genes = list()
    samples = list()

    # Create two lists from one column
    for item in seqids:
        g, s = item.split("|")
        genes.append(g)
        samples.append(s)

    # Append lists to the data frame as columns    
    df["gene"] = genes
    df["sample"] = samples
    return df

def create_allelic_pam(df):
    """
    Assign allele IDs and create an allelic PAM based on clustering results.
    """
    samples = get_unique_ids(df, "sample")
    cluster_ids = get_unique_ids(df, "cluster")  # 0, 1, 2, ...
    genes_visited = dict()  # A dictionary of genes in processed clusters. {gene : number of clusters}
    pam = pandas.DataFrame(samples, columns = ["sample"])  # Initalise the output PAM

    # Create a list about presence-absence of each allele and combine it to the output PAM as a column
    for c in cluster_ids:  # Type of elements: numpy.int64
        df_c = df[df["cluster"] == c]  # Select rows of the current cluster
        genes_c = df_c["gene"].tolist()  # All gene names should be the same when alleles are clustered under 100% nucleotide identity
        gene = genes_c[0]
        if gene in genes_visited.keys():
            genes_visited[gene] += 1
            allele = gene + "." + str(genes_visited[gene])  # Adding a suffix for making an allele name. Example result: sul1.1, sul1.2.
        else:
            genes_visited[gene] = 0  # Record a new gene encountered
            allele = gene  # The first allele of the gene will not have an extended index appended.
        pa_vec = list()  # A binary vector about presence-absence of the current allele across samples
        samples_c = df_c["sample"].tolist()  # Samples in which the current allele is detected
        for s in samples:
            pa = 1 if s in samples_c else 0
            pa_vec.append(pa)
        pam[allele] = pa_vec
    return pam

def get_unique_ids(df, col_name):
    """
    A subordinate function of create_allelic_pam for getting a list of unique and sorted values from a
    given column of input data frame df.
    """
    ids = list(df[col_name].unique())
    ids.sort(reverse = False)
    return ids

if __name__ == "__main__":
    main()
