#!/usr/bin/env python

"""
Converts vectors of items into a presence-absence matrix, assuming the first column records IDs.
This process can be considered as reversing the output of function mkFilterTSV in R package
GeneMates (https://github.com/wanyuac/GeneMates).

Input example (TSV or CSV, no header line should be contained; the file must contain at least two columns):
ID1,a,b,c,d
ID2,a,b,c
ID3,a,c,d
ID4,

Corresponding output:
ID,a,b,c,d
ID1,1,1,1,1
ID2,1,1,1,0
ID3,1,0,1,1
ID4,0,0,0,0

Command line: vec2csv.py --input [input file] --output [output file] --input_format [csv or tsv] --output_format [csv or tsv]
Example: vec2csv.py --input vectors.csv --output pam.tsv --input_format csv --output_format tsv

Dependancy: module pandas (https://pandas.pydata.org/), which can be easily installed using command line 'pip install pandas'.

Copyright (C) 2021 Yu Wan <wanyu@microbialsystems.cn>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 24 Apr 2021; the latest update: 24 Apr 2021.
"""

from argparse import ArgumentParser
import os
import sys
import pandas

def parse_arguments():
    parser = ArgumentParser(description= "Create a presence-absence matrix from a list of vectors")
    parser.add_argument("-i", "--input", dest = "input", type = str, required = True, help = "Path and name of the input file")
    parser.add_argument("-o", "--output", dest = "output", type = str, required = False, default = "pam.csv", help = "Path and name of the output file")
    parser.add_argument("-if", "--input_format", dest = "input_format", type = str, required = False, default = "csv", help = "Input format, csv/tsv")
    parser.add_argument("-of", "--output_format", dest = "output_format", type = str, required = False, default = "csv", help = "Output format, csv/tsv")
    return parser.parse_args()

def main():
    prog_args = parse_arguments()
    pam = make_pam(import_vectors(f = prog_args.input, sep = set_delimiter(prog_args.input_format)))
    pam.to_csv(prog_args.output, index = False, sep = set_delimiter(prog_args.output_format))
    return

def import_vectors(f, sep):
    """ Returns a dictionary {id : [a list of items]} """
    if not os.path.exists(f):
        print("Input error: input file " + f + " is not accessible.")
        sys.exit(1)
    else:
        d = dict()
        with open(f, "r") as input_file:
            vectors = input_file.read().splitlines()
        for v in vectors:
            fields = v.split(sep)
            d[fields[0]] = fields[1 : ]  # Assuming that the first item is the ID.
    return d

def set_delimiter(fmt):
    if fmt == "tsv":
        sep = "\t"
    else:
        sep = ","
    return(sep)

def get_item_names(d):
    """
    Pools items from the dictionary d into a set to eliminate redundancy.
    This is a subordinate function of make_pam.
    """
    colnames = set()
    for _, fields in d.items():
        colnames.update(fields)  # Add elements of a list to a set.
    return colnames

def make_pam(vectors):
    item_list = list(get_item_names(vectors))  # Gets a set of item names as column names of the output matrix
    colnames = ["ID"] + item_list
    pam = pandas.DataFrame(columns= colnames)
    n = len(item_list)
    for i, vals in vectors.items():
        if vals == "":
            pam = pam.append(pandas.Series([i] + [0] * n, index = colnames), ignore_index = True)  # [0] * n: create a list of n zeros.
        else:
            row = list()
            for j in item_list:
                if j in vals:
                    row.append(1)  # Item j is present under the current ID.
                else:
                    row.append(0)
            pam = pam.append(pandas.Series([i] + row, index = colnames), ignore_index = True)
    return pam

if __name__ == "__main__":
    main()
