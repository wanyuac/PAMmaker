#!/usr/bin/env Rscript
# Assess uncertainty of allele calls with their scores.
# MLST information in the input table of allele calls (specified by the --profiles option) will not be processed by this script.
# Command format:
#   Rscript assessAlleleCallUncertainty.R --profiles [SRST2 geotyping table] --scores [merged SRST2 scores] --output [output filename without an extension]
# Example:
#   Rscript assessAlleleCallUncertainty.R --profiles Kp_srst2__compiledResults.txt --scores mergedScores__gene.scores --output Kp_srst2__reliableCalls
# Author: Yu Wan (wanyuac@gmail.com)
# Development history: 2 Aug 2016, 16 Apr 2017
# Licence: Apache-2.0

library(stringr)  # for the function "str_extract"
library(optparse)

# Options ###############
options <- list(
    make_option("--profiles", type = "character", help = "A compiled file of SRST2 allele calls across all samples"),
    make_option("--scores", type = "character", help = "A file containing merged SRST2 scores"),
    make_option("--output", type = "character", default = "srst2_reliableAlleleCalls", help = "Name of output files without filename extension [default %default]"),
    make_option("--avg_depth", type = "integer", default = 5, help = "Minimal average depth for reliable allele calls [default %default]"),
    make_option("--edge1_depth", type = "integer", default = 2, help = "Minimal edge-1 depth [default %default]"),
    make_option("--edge2_depth", type = "integer", default = 2, help = "Minimal edge-2 depth [default %default]"),
    make_option("--neighb_depth", type = "integer", default = 2, help = "Minimal read depth neighbouring a truncation [default %default]"),
    make_option("--divergence", type = "integer", default = 10, help = "Maximal nucleotide divergence [default %default]")
)

# Functions ###############
# remove MLST results from the table
rmMLST <- function(x) {
    if ("maxMAF" %in% names(x)) {  # 12 columns of MLST results are present; maxMAF is the last column.
        x <- x[, c(1, (which(names(x) == "maxMAF") + 1) : ncol(x))]
    }
    return(x)
}

# extract flags ("*" and "?") for uncertainty and novelity of allele calls
getFlags <- function(allele) {
    flag <- str_extract(allele, "[*?]+")
    if (is.na(flag)) {
        flag <- ""
    }

    return(flag)
}

# calculate the identity of each allele call
divergence <- function(scores) {
    # It is the same as the equation in SRST2's script.
    return(round(scores$Mismatches / (scores$Size - scores$Truncated_bases) * 100, digits = 8))
}

# Note that some assured allele variants of MLST genes may become unreliable after this algorithm.
# given coverage > 90% (the default parameter of SRST2)
evaluateReliability <- function(x, avg.min = 5, edge1.min = 2, edge2.min = 2, neighb.min = 2, diverg.max = 10) {
    n <- nrow(x)
    y <- rep(FALSE, times = n)
    for (i in 1 : n) {
        r <- x[i, ]  # extract one row from the data frame x
        if (is.na(r$DepthNeighbouringTruncation)) {  # DepthNeighbouringTruncation == NA <=> Truncated_bases = 0 (no truncation)
            y[i] <- (r$Avg_depth > avg.min) & (r$Edge1_depth >= edge1.min) & (r$Edge2_depth >= edge2.min) & (r$divergence <= diverg.max)
        } else {  # with at least one truncation
            y[i] <- (r$Avg_depth > avg.min) & (r$DepthNeighbouringTruncation > neighb.min) & (r$divergence <= diverg.max)
        }
    }

    return(y)  # return a logical value
}

# edit the profile of resistance genes ----------
editAlleleCall <- function(x, sample.scores) {
    selection <- which(sample.scores$Allele == x)
    n <- length(selection)
    if (n == 1) {  # expecting to have exactly one allele call
        reliability <- sample.scores$reliable[selection]  # There must be only one hit.
        if (reliability) {  # reliable
            y <- str_extract(x, "[^?]+")  # remove the question mark when it is present
        } else {
            y <- "-"  # remove this allele call if it is unreliable. This step may produce empty columns.
        }
    } else if (n > 1) {  # If there is such an error, there could be a problem in your database.
        cat("Error: sample", sample.scores$Sample[1], "yielded multiple calls for the allele", x, ".\n")
        y <- NA
    } else {  # An n = 0 indicates that the merged score file is wrong.
        cat("Error: the allele", x, "was actually not called in the sample", sample.scores$Sample[1], ".\n")
        cat("Please check you score file and the genotype profile.\n")
        y <- NA
    }
    return(y)
}

# Main program ###############
# read arguments
opt.parser <- OptionParser(usage = "Rscript %prog [options]", option_list = options)
opts <- parse_args(opt.parser)

# read allele calls
profiles <- read.delim(opts$profiles, check.names = FALSE, stringsAsFactors = FALSE)  # genetic profiles from SRST2
profiles <- rmMLST(profiles)  # detect and remove columns of MLST results
samples <- profiles$Sample
genes <- names(profiles)[-1]  # except the first column name "Sample"
rownames(profiles) <- samples  # for the ease of editing values later

# read allele scores
scores <- read.delim(opts$scores, check.names = FALSE, stringsAsFactors = FALSE)  # allele calls of resistance genes
scores$Allele <- gsub(pattern = "__", replacement = "_", x = scores$Allele)  # shorten allele names
scores$flag <- sapply(scores$Allele, getFlags)  # extract flags from every allele ID
scores$divergence <- divergence(scores)  # calculate the divergence (to respective reference sequence) for every allele call because SRST2 does not include this piece of information in score files

# evaluate reliability of every allele call
scores$reliable <- evaluateReliability(x = scores, avg.min = opts$avg_depth,
                                       edge1.min = opts$edge1_depth,
                                       edge2.min = opts$edge2_depth,
                                       neighb.min = opts$neighb_depth,
                                       diverg.max = opts$divergence)

# remove uncertainty signs from profiles in accordance with reliability ===============
for (s in samples) {
    sample.scores <- subset(scores, Sample == s)
    for (g in genes) {
        a <- profiles[s, g]  # get an allele name
        if (a != "-") {
            profiles[s, g] <- editAlleleCall(x = a, sample.scores = sample.scores)
        }
    }
}

# Notice profiles may contain empty columns after this reliability assessment.
# To-do: remove the empty columns and report genes that loses all alleles after
# this assessment.
write.table(profiles, file = paste0(opts$output, ".txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(scores, file = paste0(opts$output, ".scores"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
