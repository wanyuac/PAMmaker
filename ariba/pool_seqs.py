#!/usr/bin/env python

"""
Pool FASTA files of ARIBA's output allele sequences into one file and append sample names to sequence IDs (in
the same format as that for SRST2's output consensus sequences).

Command demonstration:
  python -i *_genes.fna -o alleles.fna -e '_genes.fna'

Copyright (C) 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Publication: 11 Nov 2020; the latest modification: 11 Nov 2020
"""

import os
from argparse import ArgumentParser

def parse_arguments():
    parser = ArgumentParser(description = "Pool ARIBA's output allele sequences into one FASTA file and append sample names to sequence IDs",\
             epilog = "This is a helper script of R package GeneMates.")
    parser.add_argument("-i", nargs = "+", dest = "i", type = str, required = True, help = "Input FASTA files")
    parser.add_argument("-o", dest = "o", type = str, required = False, default = "alleles.fna", help = "Output FASTA file of pooled allele sequences")
    parser.add_argument("-e", dest = "e", type = str, required = False, default = "_genes.fna", help = "Filename extension to be removed for sample names")
    return parser.parse_args()

def main():
    args = parse_arguments()
    fasta_out = open(args.o, "w")
    sample_count = 0
    for fasta_in in args.i:
        f_in = open(fasta_in, "r")
        fasta_in = os.path.basename(fasta_in)
        sample = fasta_in.replace(args.e, "")  # No change applies if args.e is not found in the filename
        line = f_in.readline()
        while line:
            if line.startswith(">"):  # A header is encountered
                fields = line.split(".")  # ARIBA uses full stops as delimiters in the header line
                new_id = fields[0] + "|" + sample  # Example value: ">cluster_1|sample_1"
                seq_descr = ".".join(fields[1 : ])
                fasta_out.write(new_id + " " + seq_descr)  # Note that the newline character of this line is not stripped off by the readline method.
            else:
                fasta_out.write(line)
            line = f_in.readline()  # Till the end of the file, a None value is returned.
        f_in.close()
        sample_count += 1
    fasta_out.close()
    print("%d samples have been processed." % sample_count)
    return

if __name__ == "__main__":
    main()
