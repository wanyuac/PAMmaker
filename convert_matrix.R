# This program converts an allele matrix into a presence-absence matrix of alleles (allelic PAM).
#
# Usage:
#   Rscript convert_matrix.R modified_allele_matrix.txt out_dir
#
# Notice the output directory name/path should not be ended by a forward slash.
#
# Copyright 2017 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# Development history: 5 Sep 2016, 16 Apr 2017, 30 Nov 2017

tabulateAlleles <- function(col) {
    # This function is embeded in the next function - getAllelicMatrix.
    # For each column (a gene), it generates a presence/absence matrix for all alleles present.
    tab <- table(col)
    alleles <- names(tab)
    alleles <- alleles[alleles != "-"]  # remove "-" characters from the input table. A minus sign means absence of the corresponding allele.
    n <- length(alleles)
    mat <- matrix(nrow = length(col), ncol = n)  # create an empty matrix
    colnames(mat) <- alleles
    for (i in 1 : n) {
        mat[, i] <- as.numeric(col == alleles[i])  # presence of the allele: 1, absence: 0
    }
    return(as.data.frame(mat))
}

getAllelicMatrix<- function(m) {
    # This function returns a presence/absence matrix of every allele in the genetic matrix.
    columns <- apply(m, 2, tabulateAlleles)  # for each column (a gene) in the genetic matrix, excluding the Sample column
    pam <- data.frame(Sample = rownames(m), stringsAsFactors = FALSE)
    for (i in 1 : length(columns)) {
        pam <- cbind(pam, columns[[i]])  # bind every sub-matrix in the list "columns"
    }
    return(pam)
}

# Main program ##########
params <- commandArgs(trailingOnly = TRUE)
f <- params[1]
if (length(params) > 1) {
    outdir <- params[2]
} else {
    outdir <- "."
}

m.genes <- read.delim(file = f, check.names = FALSE, stringsAsFactors = FALSE)  # The input matrix must not contain any columns that are comprised of only "-".
samples <- m.genes$Sample
m.genes <- as.matrix(m.genes[, -1])
rownames(m.genes) <- samples

# remove empty columns (a few genes may lose all allele calls due to reliability assessment)
null.counts <- as.integer(apply(m.genes, 2, function(c) sum(c != "-")))
m.genes <- m.genes[, null.counts > 0]

# make a presence/absence matrix
m.alleles <- getAllelicMatrix(m.genes)  # convert the allele matrix into a presence/absence allelic matrix
write.table(m.alleles, file = paste(outdir, "allele_paMatrix.txt", sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE)
