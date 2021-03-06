# This script produces summary statistics for the maxMAF information per allele in the output of SRST2.
# It first makes a matrix of maxMAF (maximum minor-allele frequencies) based on allele occurrence in samples,
# then it summarise maxMAF and produces a data frame for results.
#
# The allelic presence/absence matrix (PAM) and the score file are tab-delimited files produced by mk_allele_matrix.sh
# and assessAlleleCallUncertainty.R (under the directory PAMmaker/reliability_assessment), respectively. There are
# extended allele identifiers attached to allele names in the PAM, such as StrA_1501.173.
#
# Since the score file generated by mk_allele_matrix.sh also includes scores of
# unreliable allele calls, this script expects that some allele names from the
# score file to be missing in the PAM of reliable allele calls. Similarly, it does
# not require all strain names in the PAM to be present in the score file as the
# PAM may be merged from several files that involve complete genomes (no SRST2 is
# need in this case).
#
# Usage:
# This script can be either run with the following R script command
#   Rscript --vanilla maxMAFmatrix.R [allelic PAM] [score file] [prefix for output filenames]
# or run the functions in R (simply run the function definition part first).
#
# Author: Yu Wan <wanyuac@gmail.com>
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# First edition: 31 Jan 2018; the latest edition: 6 Aug 2018

# Parameter list ###############
library(optparse)

options <- list(
    make_option(c("-a", "--apam"), type = "character", dest = "apam",
                help = "An allelic presence-absence matrix"),
    make_option(c("-s", "--scores"), type = "character", dest = "scores",
                help = "A tab-delimited file of alignment scores for reliable allele calls"),
    make_option(c("-r", "--replacement"), type = "character", dest = "rep", default = NULL,
                help = "A tab-delimited file for allele name substitutions"),
    make_option(c("-o", "--output"), type = "character", dest = "output", default = "./allele",
                help = "Path and prefix for output files [default %default]")
)

# Define functions ###############
.importMappingTable <- function(tab_file) {
    # Read the table mapping allele names for replacement. It is usually named
    # allele_name_replacement.txt and produced by clustering_allele_variants.py.
    if (!is.null(tab_file)) {
        print("Reading the mapping table for allele name replacement.")
        mapping <- read.delim(tab_file, stringsAsFactors = FALSE)  # column names: Sample, Prev_allele, Var, New_allele
        mapping <- mapping[, c("Sample", "Prev_allele", "New_allele")]
    } else {
        mapping <- NULL
    }

    return(mapping)
}

.searchAlleleID <- function(strain, allele, mapping) {
    # This is a subordinate function of .replaceAlleleIDs.
    i <- which(mapping$Sample == strain & mapping$AlleleID_base == allele)
    if (length(i) > 0) {
        new_id <- mapping$Allele[i]
    } else {
        print(paste0("Warning: allele ", allele,
                    " from the score table is absent in the PAM for the sample ",
                    strain, "."))
        new_id <- NA
    }

    return(new_id)  # This vector may contain NA's when some alleles in the score file are missing in the PAM.
}

.correctAlleleNames <- function(score_tab, repl) {
    # This is a subordinate function of .replaceAlleleIDs.
    for (i in 1 : nrow(repl)) {
        r <- repl[i, ]
        strain <- r$Sample
        a_prev <- r$Prev_allele
        a_new <- r$New_allele
        j <- which(score_tab$Sample == strain & score_tab$Allele == a_prev)
        if (length(j) > 0) {  # equals one actually, as there cannot be multiple hits
            score_tab$Allele[j] <- a_new
        }
    }

    return(score_tab$Allele)
}

.replaceAlleleIDs <- function(pam, score_tab, repl) {
    # Substitute allele names in the score table with allele names in the PAM
    # For example, to replace "MphA_1663" in strain ERR134530 with "MphA_1663.141".
    # This algorithm works only when there is only a single allele call per gene in each strain.
    score_tab <- score_tab[, c("Sample", "Allele")]

    # Drop "?/*/*?" from allele names in the score file assuming all of the alleles listed in this file are reliable calls.
    # By far, the allele names in this table are original names in the SRST2 output.
    score_tab$Allele <- str_replace(string = score_tab$Allele, pattern = "([?*])+", replacement = "")

    # Replace allele names in the score table with those in the mapping table
    if (!is.null(repl)) {
        score_tab$Allele <- .correctAlleleNames(score_tab, repl)
    }

    # make a strain-allele table (mapping) from the PAM
    strains <- rownames(pam)
    alleles_all <- colnames(pam)  # They are corrected allele names with extended indices when repl != NULL.
    mapping <- data.frame(Sample = character(0), Allele = character(0), stringsAsFactors = FALSE)
    for (s in strains) {
        alleles <- alleles_all[as.logical(pam[s, ])]
        n <- length(alleles)
        if (n > 0) {
            mapping <- rbind.data.frame(mapping,
                                        data.frame(Sample = rep(s, times = length(alleles)), Allele = alleles,
                                                   stringsAsFactors = FALSE))
        }
    }

    # make a new column for allele IDs without the extended index. For instance, convert MphA_1663.141 back into MphA_1663.
    # Herein, the base name (e.g., MphA_1663 for MphA_1663.141) of allele IDs must not have any dots.
    mapping$AlleleID_base <- sapply(mapping$Allele, function(a) strsplit(x = a, split = ".", fixed = TRUE)[[1]][1])

    # Substitute allele names in the score table
	# This step introduces NA's when the sample or allele is missing in the data frame mapping (which comes from the allelic PAM).
    new_ids <- as.character(mapply(.searchAlleleID, score_tab$Sample, score_tab$Allele, MoreArgs = list(mapping = mapping)))

    return(new_ids)
}

maxMAFmatrix <- function(apam, scores, sep = "\t", mapping = NULL) {
    require(stringr)

    # Preparation of data ###############
    # A PAM is a matrix of only zeros and ones and it is expected to only
    # contain reliable allele calls.
    cls <- class(apam)
    if (cls == "character") {
        apam <- as.matrix(read.table(file = apam, header = TRUE, sep = sep, row.names = 1,
                                     check.names = FALSE, stringsAsFactors = FALSE))  # take the first column for row names
        strains <- rownames(apam)
    } else if (cls == "data.frame") {  # convert it to a binary matrix
        strains <- apam[, 1]  # assuming the first column stores sample names
        apam <- as.matrix(apam[, -1])
        rownames(apam) <- strains
    } else if (cls == "matrix") {
        strains <- rownames(apam)
    } else {
        stop("Input error: the allelic PAM must be a file name, a data frame or a matrix.")
    }
    alleles <- colnames(apam)

    # Import the score table as a data frame
    # Since column names of this matrix are replaced allele names when mapping != NULL,
    # they cannot match to allele names in the score table. It is hence necessary
    # to replace the allele name in the score table with those in the allelic PAM.
    cls <- class(scores)
    if (cls == "character") {
        scores <- read.table(file = scores, header = TRUE, sep = sep, stringsAsFactors = FALSE)  # first column: strain names
    } else if (cls == "data.frame") {
        print("The argument scores is already a data frame. No further action is needed.")
    } else if (cls == "matrix") {
        scores <- data.frame(Sample = rownames(scores), scores, stringsAsFactors = FALSE)
        rownames(scores) <- NULL
    } else {
        stop("Input error: the score table must be a file name, a data frame or a matrix.")
    }
    scores$Allele <- .replaceAlleleIDs(apam, scores, mapping)  # Some allele names become NA afterwards.

    # It is not surprising to see some alleles in the score table do not match to any alleles
    # in the PAM because the PAM may be filtered for some purpose before being applied to this
    # analysis. Nonetheless, since we only use the PAM for GeneMates, unmatched alleles in the
    # score table do not matter for expected results. As such, this function discards those NA
    # entries in the following step.
    scores <- subset(scores, !is.na(Allele))

    # Generate a maxMAF matrix ###############
    m <- matrix(NA, nrow = nrow(apam), ncol = ncol(apam))  # initialising the output matrix
    rownames(m) <- strains
    colnames(m) <- alleles

    # search maxMAF for each allele call: the following algorithm is faster than accessing pam[s, a].
    for (s in strains) {
        pa <- apam[s, ]  # a named vector of integers
        alleles_s <- names(pa)[as.logical(pa)]  # alleles in the current strain
        if (length(alleles_s) > 0) {
            for (a in alleles_s) {
                i <- which(scores$Sample == s & scores$Allele == a)
                if (length(i) > 0) {
                    m[s, a] <- round(scores$maxMAF[i], digits = 6)  # XX.XXXX%
                } else {
                    # Sometimes a strain in the PAM is not present in the score
                    # file as the latter file contains scores of unreliable alleles
                    # as well. A typical example is that this strain has a complete
                    # genome, hence the gene screen is carried out using BLAST.
                    print(paste("Note: allele", a, "of the sample", s,
                                "in the PAM is absent in the score file.", sep = " "))  # then keep the entry as NA
                }
            }
        }  # Else, keep the whole row as of all NA's (when no allele is found in the current strain).
    }

    # return the maxMAF matrix, PAM and the score table whose allele names have been replaced.
    return(list(m = m, pam = apam, scores = scores))
}

summariseMaxMAF <- function(m, pam) {
    # m: a matrix of maxMAF; pam: an allelic PAM, which is used for counting the number of occurrence
    # Initialisation
    alleles <- colnames(pam)

    # Construct the summary table
    # Count_score: number of alleles having scores and maxMAF information. Count_score <= Count.
    sumry <- data.frame(Allele = character(0), Count = integer(0), Count_score = integer(0),
                        Max = numeric(0), P75 = numeric(0), Median = numeric(0),
                        P25 = numeric(0), Min = numeric(0), stringsAsFactors = FALSE)

    # Compute summary statistics
    for (a in alleles) {
        n <- sum(pam[, a])  # number of occrrence events, given 1 for presence and 0 for absence
        maxmafs <- m[, a]
        maxmafs <- maxmafs[!is.na(maxmafs)]
        k <- length(maxmafs)
        if (k > 0) {
            maxmaf_qu <- as.numeric(quantile(maxmafs, probs = c(0, 0.25, 0.5, 0.75, 1)))
        } else {  # no maxmaf information for this allele at all
            maxmaf_qu <- rep(NA, times = 5)
        }

        sumry <- rbind.data.frame(sumry,
                                  data.frame(Allele = a, Count = n, Count_score = k,
                                             Max = maxmaf_qu[5], P75 = maxmaf_qu[4],
                                             Median = maxmaf_qu[3], P25 = maxmaf_qu[2],
                                             Min = maxmaf_qu[1], stringsAsFactors = FALSE),
                                  stringsAsFactors = FALSE)
    }
    sumry <- sumry[order(sumry$Max, decreasing = TRUE), ]

    return(sumry)
}

# Main program ###############
# Read arguments
opt.parser <- OptionParser(usage = "Rscript %prog [options]", option_list = options)
argv <- parse_args(opt.parser)

# Import the substition table for allele names when it is specified
mapping <- .importMappingTable(argv$rep)  # returns NULL when argv$rep is unspecified

# Generate summary statistics
maf <- maxMAFmatrix(apam = argv$apam, scores = argv$scores, sep = "\t", mapping = mapping)  # returns a list with elements m, pam and scores
maf_summary <- summariseMaxMAF(m = maf[["m"]], pam = maf[["pam"]])

# Save three kinds of results
# Do not need to save PAM because it is not changed in this script.
m <- cbind.data.frame(Sample = rownames(maf[["m"]]), as.data.frame(maf[["m"]]))
rownames(m) <- NULL
write.table(x = m, file = paste0(argv$output, "__maxMAFmatrix.tsv"), sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(x = maf[["scores"]], file = paste0(argv$output, "__renamedScores.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(x = maf_summary, file = paste0(argv$output, "__maxMAFsummary.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
