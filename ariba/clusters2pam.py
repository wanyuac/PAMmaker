#!/usr/bin/env python

"""
Converts output TSV of utility/tabulate_cdhit.py to an allelic presence-absence matrix (PAM) and saves
it in TSV format. This script assumes consensus sequences of ARIBA hits are clustered under a minimum
nucleotide identity that is sufficently high (e.g., 100%) so reference sequences are the same within
the same cluster.

Command demonstration:
    python clusters2pam.py -i alleles.clstr.tsv -om alleles_pam.tsv -ot alleles_clstr_updated.tsv

Dependency: module pandas, Python v3

Copyright (C) 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Publication: 11 Nov 2020; the latest modification: 30 Dec 2020
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
    parser.add_argument("-ot", dest = "ot", type = str, required = False, default = "clusters.tsv", help = "Output table about sequence clusters and allele labels")
    
    return parser.parse_args()


def main():
    args = parse_arguments()

    # Import the input TSV file and make a nine-column data frame
    input_tsv = args.i
    if os.path.exists(input_tsv):
        clusters = parse_seqid(pandas.read_csv(input_tsv, sep = "\t", index_col = None))  # Import the table as a data frame of six columns and parse column 'seqid'
    else:
        print("Error: the input file is not accessible.", file = sys.stderr)
        sys.exit(1)
    
    # Assign allele IDs and create an allelic PAM from data frame "clusters"
    pam, tab = create_allelic_pam(clusters)
    pam.to_csv(args.om, header = True, index = False, sep = "\t", float_format = None)
    tab.to_csv(args.ot, header = True, index = False, sep = "\t", float_format = None)
    
    return


def parse_seqid(df):
    """ 
    Parses column 'seqid' into two new columns and returns a data frame of eight columns.
    """
    seqids = df["seqid"].tolist()  # Convert a column into a list
    ref_alleles = list()  # Names of reference alleles
    var_code = list()  # A binary list
    samples = list()  # Sample names
    
    # Create three lists from one column
    for item in seqids:
        a, v, s = item.split("|")  # Reference allele, variant code, sample name
        ref_alleles.append(a)
        var_code.append(int(v))  # Otherwise, each var_code is a character.
        samples.append(s)
    
    # Append these new lists to the data frame as additional columns
    df["sample"] = samples
    df["ref_allele"] = ref_alleles
    df["variant"] = var_code
    
    return df


def create_allelic_pam(df):
    """
    Assign allele IDs and create an allelic PAM based on clustering results.
    Algorithm for assigning an allele ID for each cluster:
        1. Use the name of the reference allele having exactly matched hits; (So the algorithm deals with scenarios
           where multiple reference alleles are present in the same cluster)
        2. Else, add an extended index (1, 2, 3, ...) to the allele name.
    """
    cluster_indices = get_unique_elements(df, "cluster")  # 0, 1, 2, ...
    samples = get_unique_elements(df, "sample")  # All samples from ARIBA results
    variants = dict()  # A dictionary of variants encountered. {referance allele : count of unique variants}
    pam = pandas.DataFrame(samples, columns = ["sample"])  # Initalise the output PAM
    tab = pandas.DataFrame([], columns = list(df.columns))  # Create an empty data frame of specific columns
    
    # Populate data frames 'tab' and 'pam'
    for c in cluster_indices:  # Type of elements: numpy.int64; Iterate through every cluster index.
        df_c = df[df["cluster"] == c]  # Extract rows of the current cluster
        df_c.reset_index(drop = True)
        
        # Determine the reference allele for the current cluster
        ref_allele, var_code = determine_ref_allele(df_c[["ref_allele", "variant"]], c)  # Reference-allele name and variant code for the current cluster
        if var_code == 0:  # Some or all hits are exact matches to the chosen reference allele
            a = ref_allele
        else:  # Approximate match: append an extended allele identifier to the reference-allele name
            if ref_allele in variants.keys():
                variants[ref_allele] += 1
                a = ref_allele + ".var" + str(variants[ref_allele])  # Adding a suffix for making an allele name. Example result: sul1.1, sul1.2.
            else:
                variants[ref_allele] = 1  # First variant of the reference allele is encountered
                a = ref_allele + ".var1"
        df_c = df_c.assign(allele = [a for i in range(0, df_c.shape[0])])  # To append a column to the data frame. The shape method returns row and column counts. Do not use pandas.Series, which does not work correctly here.
        tab = pandas.concat([tab, df_c], ignore_index = True, sort = False)  # Ignoring indexes on the concatenation axis: https://pandas.pydata.org/pandas-docs/stable/user_guide/merging.html
        
        # Add a column to the PAM
        pa_vec = list()  # Create a binary list about presence-absence of the current allele
        samples_c = get_unique_elements(df_c, "sample")  # Samples in the current sequence cluster
        for s in samples:
            pa = 1 if s in samples_c else 0
            pa_vec.append(pa)
        pam[a] = pa_vec
    
    return pam, tab


def get_unique_elements(df, col_name):
    """
    A subordinate function of create_allelic_pam for getting a list of unique and sorted values from a
    given column of input data frame df.
    """
    vals = list(df[col_name].unique())
    vals.sort(reverse = False)
    
    return vals


def determine_ref_allele(df, cluster_index):
    """
    This is a subordinate function of create_allelic_pam and it returns the first element in a list of
    allele names extracted from a column of data frame df_c["ref_allele"]. All names of reference alleles
    should be the same when hits are clustered under 100% nucleotide identity and a warning is raised when
    this is not the case.
    """
    # Check whether the current cluster comprises hits to more than one reference alleles
    if len(set(df["ref_allele"].tolist())) > 1:  # Normally, this is a rare scenario.
        print("Warning: reference alleles in cluster %d differ." % cluster_index)
    
    # Choose the first reference allele having an exact hit
    exact_matches = df.loc[df["variant"] == 0]  # 0 (exact match) or 1 (variant)
    if len(exact_matches) > 0:  # There is at least one exact hit
        allele = exact_matches["ref_allele"].values[0]  # Get the value of the first cell in this column
        var_code = 0
    else:  # No exact match is present
        allele = df["ref_allele"].values[0]
        var_code = 1
                                                   
    return allele, var_code


if __name__ == "__main__":
    main()
