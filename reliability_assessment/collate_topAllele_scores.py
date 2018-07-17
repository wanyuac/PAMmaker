"""
This script merges score files (*.scores) produced by SRST2. It only filters out score information of alleles that
are not called (namely, non-top-score alleles) by SRST2 and only keeps scores of allele calls (that is, top alleles)
in every sample. This is a utility that the current version of SRST2 does not provide.

As for the procedure of this script, first, it reads files of allele calls in order to construct a table of allele
calls in every sample. Then, it reads every score file and only add those corresponding to the allele calls into
a new tab-delimited file. Finally, two files of top-allele scores will be stored under the current working directory:
one for MLST allele calls and the other for non-MLST allele calls. You can concatenate them using bash commands or
other software for data analysis.

This script only works for SRST2-formatted gene databases (such as ARGannot.r1.fasta) and MLST-allele databases. You
need to modify this script if another kind of header formats are used in your database.

Typical input files and expected formats of their names:
    1. MLST allele calls: [sample name]_[prefix]__mlst__[MLST database name]__results.txt
    2. non-MLST allele calls: [sample name]_[prefix]__genes__[gene database name]__results.txt
    3. MLST allele-score files: [sample name]_[prefix]__[sample name].[MLST database name].scores
    4. non-MLST allele-score files: [sample name_[prefix]__[sample name].[gene database name].scores
    
Requirements for options: at least one pair of options, namely, mlst_calls and mlst_scores, or allele_calls and
allele_scores, must be configured.

Example command line:
    python collate_topAllele_scores.py --mlst_calls *_Kp__mlst__mlstDB__results.txt --mlst_scores *_Kp__*.mlstDB.scores
    --allele_calls *_Kp__genes__geneDB__results.txt --allele_scores *_Kp__*.geneDB.scores --mlst_delimiter '_' --prefix Kp

Other notes:
    Please check that whether the way for parsing file names in the function merge_allele_scores works for your
    file names because prefixes in the file names would violate the format assumption of this function:
        sample = (file_name.split(".")[0]).split("__")[-1]
    
    It is recommended that not to leave the prefix argument empty (namely, "") becasue this script erases contents
    of files if these files have the same name as those of the output files. You may override previous outputs if
    there is no prefix configured.

Python version: 2.7.10

Author: Yu Wan (wanyuac@gmail.com)
Development history: 8 April 2016, 17 July 2016, 16 August 2016, 27 June 2017, 28 Nov 2017 and 5 Dec 2017
"""

import os
import re
import sys
import copy
from collections import defaultdict
from argparse import ArgumentParser

def parse_arguments():
    # a function that creates an ArgumentParser object, which reads arguments for this script
    # use the command python merge_scores.py -h to show the list and description of arguments
    
    parser = ArgumentParser(description = "Merge scores of allele calls")
    parser.add_argument("--mlst_calls", type = str, nargs="+", required = False, default = "", help = "A tab-delimited output of SRST2 consisting of MLST calls")
    parser.add_argument("--mlst_scores", type = str, nargs="+", required = False, default = "", help = "*.scores files produced by SRST2 for MLST genes")
    parser.add_argument("--allele_calls", type = str, nargs="+", required = False, default = "", help = "A tab-delimited output of SRST2 consisting of allele calls for non-MLST genes")
    parser.add_argument("--allele_scores", type = str, nargs="+", required = False, default = "", help = "*.scores files produced by SRST2 for non-MLST genes")
    parser.add_argument("--mlst_delimiter", type = str, required = False, default = "-", help = "The delimiter separating gene names and allele numbers in your MLST database (default: '-')")
    parser.add_argument("--prefix", type = str, required = False, default = "", help = "The prefix of every output file")
    parser.add_argument("--no_unique_allele_id", action = "store_true", help = "Turn this option on when allele calls do not have the unique allele identifiers attached.")  # only applies to non-MLST allele calls
    
    return parser.parse_args()

def main():
    args = parse_arguments()  # read arguments
    with_unique_allele_id = not args.no_unique_allele_id
    
    # Validate arguments: at least one pair of arguments have to be completely configured
    if (args.mlst_calls == "" or args.mlst_scores == "") and (args.allele_calls == "" or args.allele_scores == ""):
        sys.exit("Error: both argument pairs for MLST results and genotyping results are incomplete.")  # exit with an error message sent to the sys.stderr object
    
    # Collate scores of MLST allele calls if the pair of arguments for processing MLST results is complete
    # This stage produces a file called [prefix]__mlst.scores.
    if (args.mlst_calls != "") and (args.mlst_scores != ""):
        merge_allele_scores(allele_table = read_allele_calls(files = args.mlst_calls, allele_type = "mlst",\
                                                             mlst_delimiter = args.mlst_delimiter),\
                            score_files = args.mlst_scores, allele_type = "mlst", output_prefix = args.prefix)
    
    # Collate scores of non-MLST allele calls if the pair of arguments for processing these results is complete
    # This stage produces a file called [prefix]__gene.scores.
    if (args.allele_calls != "") and (args.allele_scores != ""):
        merge_allele_scores(allele_table = read_allele_calls(files = args.allele_calls, allele_type = "gene",\
                                                             with_unique_allele_id = with_unique_allele_id),\
                            score_files = args.allele_scores, allele_type = "gene", output_prefix = args.prefix,\
                            with_unique_allele_id = with_unique_allele_id)
    
    return

def search(query, subject):
    # return the index of an element in a list (the subject) that matches a query value
    try:
        i = subject.index(query)  # the "subject" must be a list
        return i
    except ValueError:
        return -1  # not found

def read_allele_calls(files, allele_type, mlst_delimiter = "-", with_unique_allele_id = True):
    """
    Read allele calls of every sample and store this information in a table
    allele_type: either "mlst" or "gene"
    The argument mlst_delimiter will not be used when the allele_type is "gene".
    """
    
    allele_table = defaultdict(dict)  # [sample][gene] = allele
    
    for file_name in files:
        content = open(file_name, "rU").read().splitlines()  # read the whole file into a list of rows
        n_rows = len(content)  # number of rows in the target file
        
        if n_rows == 2:  # a proper file should always contain a header row and a value row
            for i in [0, 1]:
                # i = 0: the header line consisting of column names; i = 1: the row listing allele calls.
                content[i] = content[i].split("\t")  # Every element of content[i] becomes a list.
                
            sample = content[1][0]  # obtain the sample name from the first column in the second row
            
            if allele_type == "mlst":
                # fields of the content for MLST allele calls: Sample\tST\tGene1\t...\tGeneN\tmismatches\tuncertainty\tdepth\tmaxMAF
                gene_field_boundary = search("mismatches", content[0])  # Among column names, find out the index of the first element following gene names
                if gene_field_boundary > 0:  # if gene names are present
                    genes = content[0][2 : gene_field_boundary]  # obtain gene names, skipping the first two columns (Sample and ST)
                    allele_numbers = content[1][2 : gene_field_boundary]  # obtain digits representing every allele
            else:
                # fields of the content for non-MLST allele calls: Sample\tGene1\t...\tGeneN
                # Notice in some exetreme cases, where a tiny reference database is used and some samples do not have any allele calls,
                # the target file will only contain two rows "Sample\nsample". This does not apply to MLST results. So the following flow control
                # aims to deal with this case.
                if len(content[1]) == 1:
                    print("Reading allele calls: the file " + file_name + " is skipped because it does not contain any allele hits.")
                    continue  # skip the current sample; it will not be added into the dictionary of allele_table
                else:
                    genes = content[0][1 : ]  # Here, every gene ID contains a product name: [gene name]_[product]. Do not put [2 : ] here because the gene field starts at the second column.
                    allele_numbers = content[1][1 : ]  # Different to MLST results, these allele numbers are actually strings rather than digits, such as "AmpH_634*" or "-", in this scenario.
                    if with_unique_allele_id:  # when every allele call has a unique allele identifier, such as the "_1298" in "SHV-28__1298"
                        for j in range(0, len(allele_numbers)):
                            """
                            To insert an underscore character ("_") between the allele name and the internal sequence ID to fit the format of allele
                            names in the gene database and the score file. For example, "SHV-28_1298" becomes "SHV-28__1298"
                            after this process. However, sometimes an allele name contains underscores as well, such as "espF2_1-3_284",
                            a simple substitution of '_' does not work. Hence I need to find out an alternative way to achieve this.
                            """
                            allele = allele_numbers[j]
                            if allele.count("_") > 1:  # more than one underscores are present
                                pos = find_all(allele, '_')[-1]
                                allele_numbers[j] = "__".join([allele[0 : pos], allele[(pos + 1) : ]])  # skip the last "_"
                            else:
                                allele_numbers[j] = re.sub("_", "__", allele)
            
            # go through every gene of the current sample and store its corresponding allele name
            # Here, a gene name is always accompanied with an allele name. As a result, they share the same index.
            for j in range(0, len(genes)):
                if allele_numbers[j] != "-":
                    """
                    skip the gene without an allele call
                    Uncertainty marks ('*', '?', or '*?') are also included in resultant allele names.
                    An MLST allele name = gene name + delimiter + allele number, such as 'gapA_1'.
                    """
                    if allele_type == "mlst":
                        allele_table[sample][genes[j]] = genes[j] + mlst_delimiter + allele_numbers[j]
                    else:
                        allele_table[sample][genes[j]] = allele_numbers[j]  # It is unnecessary to attach a gene name to an allele number for non-MLST genes.
                        
                        #allele_table[sample][genes[j]] = genes[j] + "__" + allele_numbers[j]
                        """
                        Although it is logically correct, do not use the command above for ARGannot.r1.fasta and earlier versions of the database because there is a problem
                        in the definition line of the sequence "55__CmlA_PheCmlA5__CmlA5__1538", which causes a discrepancy between gene names in the allele profile and the
                        score file. You will lose alleles if you use this command, where the string "55__CmlA_Phe__CmlA5__1538" is not in "55__CmlA_PheCmlA5__CmlA5__1538".
                        """
        else:
            print("Reading allele calls: the file " + file_name + " is skipped because it does not contain any allelic information.")
    
    return allele_table


def merge_allele_scores(allele_table, score_files, allele_type, output_prefix, with_unique_allele_id = True):
    # filter scores in accordance with the allele table and save results to a text file
    # allele_type = "mlst" or "gene"
    
    with open(score_files[0], "rU") as f:
        header = f.readline().rstrip("\n").split("\t")  # read the header line of the first score file, which is actually the same in all score files
    
    # initialise the result file
    output_file_name = output_prefix + "__" + allele_type + ".scores"
    open(output_file_name, "w").close()  # create an empty file or erase an existing content
        
    output_file = open(output_file_name, "a")  # open this file again but for appending scores
    print >> output_file, "\t".join(["Sample"] + header)  # write the header into the output file
    
    valid_samples = allele_table.keys()  # all samples that have at least a single allele call
    
    for file_name in score_files:
        """
        Assumed format of the file name: [sample name]_[prefix]__[sample name].[MLST/gene database name].scores
        Note that this assumption may be violated by the structure of prefixes in your input file names.
        Please check or edit this command to fit your format:
            sample = (file_name.split(".")[0]).split("__")[-1]
        """
        sample = (os.path.basename(file_name).split(".")[0]).split("__")[-1]  # obtain the sample name from each file name
        
        if sample in valid_samples:
            # obtain genes and alleles of the current sample
            print("Transferring scores of allele calls in the sample " + sample + " into the result file.")
            genes = allele_table[sample].keys()  # a list of gene names of the current sample
            alleles = allele_table[sample].values()  # all allele names of this sample, including signs for uncertainty and variants
            
            # create a new list of allele names without uncertainty and variant marks to match those from the score file
            # Because score files do not contain these marks.
            alleles_unsigned = copy.deepcopy(alleles)
            for i in range(0, len(alleles_unsigned)):
                alleles_unsigned[i] = re.sub("[*?]", "", alleles_unsigned[i])  # removes "*" and "?" characters from the string "allele"
            
            content = open(file_name, "rU").read().splitlines()[1 : ]  # omit the header line in every score file
            
            # filtering every score files: go through every line of each file and only keep those called by SRST2
            alleles_found = []  # alleles in the current sample that have been found in the corresponding score file
            for line in content:
                fields = line.split("\t")
                """
                fields[0]: the allele name in the score file. No uncertainty sign is attached to this name.
                The following codes work if the current allele is recorded in the table of allele calls. Obviously, un-called alleles (denoted by '-') will not present
                in the merged score file.
                These codes also have to deal with allele names with uncertainty marks in the allele table.
                """
                if allele_type == "mlst":
                    j = search(fields[0], alleles_unsigned)
                    if j != -1:  # equivalent to "if field[0] in alleles_unsigned"
                        a = alleles[j]
                        fields[0] = a  # obtain the allele name with an uncertainty sign as the new allele name
                        alleles_found.append(a)
                        print >> output_file, "\t".join([sample] + fields)  # then append this line with the sample name into the output file
                else:
                    for j in range(0, len(alleles_unsigned)):  # go through every allele call of the current sample
                        a = alleles_unsigned[j]
                        b = fields[0]
                        pos = b.find(a)  # Whether the allele name is present in the first element, and where is it?
                        if pos != -1:  # if the former is a substring of the latter, then there may be an exact match
                            alleles_found.append(a)
                            if with_unique_allele_id:
                                if a == b[pos : ]:  # This is an exact match.
                                    """
                                    For example, the allele AadA24__1609 is found in both 229__AadA_AGly__AadA24__1609 (exact match) and
                                    229__AadA_AGly__AadA24__16091 (partial match), causing an error in the output.
                                    """
                                    fields[0] = alleles[j]  # obtain the allele name with an uncertainty sign as the new allele name
                                    print >> output_file, "\t".join([sample] + fields)
                            else:
                                """
                                In some extreme cases where no allele names are the same, SRST2 do not print the unique allele identifiers.
                                For example, the allele "intI1__1" is only shown as "intI1" rather than "intI1_1" in the output.
                                This is particularly possible when the reference database is small. Otherwise, using this option causes the
                                problem of partial matches. So use this option only when you can confirm that there is no unique allele
                                identifier in your results.
                                """
                                fields[0] = alleles[j]
                                print >> output_file, "\t".join([sample] + fields)
                 
            # Check whether there are any alleles in the gene table unmatched to those in the score file or not
            if len(alleles_found) < len(alleles_unsigned):
                missing_alleles = list(set(alleles_unsigned) - set(alleles_found))
                print("Warning: %i alleles is/are not found in the score file of the sample %s: %s." % (len(missing_alleles),\
                                                                                                       sample,\
                                                                                                       ",".join(missing_alleles)))
        else:
            print("Processing scores: the sample " + sample + " is skipped because it does not have any allele calls.")
                
    output_file.close()
    
    return


def find_all(s, ch):
    """ Find out all indices of the character ch in a string s """
    return [i for i, ltr in enumerate(s) if ltr == ch]


if __name__ == "__main__":
    main()
    