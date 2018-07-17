# This script merges two SRST2 gene-content tables (tab-delimited) into a single one.
# The first column of each table stores sample names.
#
# Example command:
# Rscript merge_gene_tables.R --in1 table1.txt --in2 table2.txt --out merged.txt
#
# Copyright 2018 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# First and the latest edition: 12 Apr 2018 at T3 terminal, Chongqing Jiangbei International Airport, China.

# Read arguments ###############
library(optparse)

options = list(
    make_option("--in1", type = "character", action = "store", default = NULL,
                help = "The first gene-content table", metavar="character"),
    make_option("--in2", type = "character", action = "store", default = NULL,
                help = "The second gene-content table", metavar = "character"),
    make_option("--out", type = "character", action = "store", default = FALSE,
                help = "The merged table", metavar = "character"),
    make_option("--sort_col", type = "logical", action = "store_true", default = FALSE,
                help = "Sort columns of the output table in an ascending order of column names",
                metavar = "logical")
)

opt.parser = OptionParser(option_list = options)
args = parse_args(opt.parser)

# Import tables ###############
t1 <- read.delim(file = args$in1, header = TRUE, sep = "\t", row.names = 1,
                 check.names = FALSE, stringsAsFactors = FALSE)
t2 <- read.delim(file = args$in2, header = TRUE, sep = "\t", row.names = 1,
                 check.names = FALSE, stringsAsFactors = FALSE)

print(paste("Table 1 has", nrow(t1), "rows.", sep = " "))
print(paste("Table 2 has", nrow(t2), "rows.", sep = " "))

# Merge table 2 into table 1 ###############
mergeTables <- function(t1, t2, sort_col = FALSE) {
    # Assuming row names are provided. No "Sample" column.
    samples1 <- rownames(t1)
    samples2 <- rownames(t2)

    # Analyse column names
    items1 <- names(t1)  # the rest of column names of the first table
    items2 <- names(t2)  # the rest of column names of the second table
    items_shared <- intersect(x = items1, y = items2)

    # Check if data types of both tables are the same
    c1 <- class(t1[1, items1[1]])
    c2 <- class(t2[1, items2[1]])
    if (c1 == c2) {
        if (c1 == "character") {  # gene-content table
            default_val <- "-"
        } else if (c1 == "integer") {  # allelic presence-absence matrix
            default_val <- 0
        } else {
            default_val <- NA
        }
    } else {
        stop("Error: data types of both tables must be the same.")
    }

    # Merge tables
    if (length(items_shared) == 0) {
        # The simpliest scenario: table 2 does not overlap with table 1 at all.
        # Directly append the table 2 to table 1.
        print("Table 1 and 2 do not share any items.")
        t1_add <- matrix(data = default_val,
                         nrow = nrow(t1), ncol = length(items2),
                         dimnames = list(samples1, items2))  # additional columns of t1
        t1 <- cbind.data.frame(t1, as.data.frame(t1_add, stringsAsFactors = FALSE))
        t2_add <- matrix(data = default_val,
                         nrow = nrow(t2), ncol = length(items1),
                         dimnames = list(samples2, items1))
        t2 <- cbind.data.frame(as.data.frame(t2_add, stringsAsFactors = FALSE), t2)
        t_out <- rbind.data.frame(t1, t2, stringsAsFactors = FALSE)
    } else {
        # When there are some items shared.
        items2_uniq <- setdiff(x = items2, items_shared)  # unique column names of table 2
        if (length(items2_uniq) > 0) {
            print("Table 2 has unique items.")
            t1_add <- matrix(data = default_val,
                             nrow = nrow(t1), ncol = length(items2_uniq),
                             dimnames = list(samples1, items2_uniq))
            t1 <- cbind.data.frame(t1, as.data.frame(t1_add, stringsAsFactors = FALSE))

            # Transfer values of t2 to t2_add for the shared items
            t2_add <- matrix(data = default_val,
                             nrow = nrow(t2), ncol = length(items1),
                             dimnames = list(samples2, items1))
            t2_add <- as.data.frame(t2_add, stringsAsFactors = FALSE)
            for (c in items_shared) {
                t2_add[, c] <- t2[, c]
            }

            # Concatenate the unique part of t2 to t1
            t2 <- cbind.data.frame(t2_add, t2[, items2_uniq])  # row names of new t2 equals those of t2_add
            t_out <- rbind.data.frame(t1, t2, stringsAsFactors = FALSE)
        } else {
            # When items2 is a subset of items1.
            print("Items of table 2 constitute a subset of those of table 1.")
            t2_add <- matrix(data = default_val,
                             nrow = nrow(t2), ncol = length(items1),
                             dimnames = list(samples2, items1))
            t2_add <- as.data.frame(t2_add, stringsAsFactors = FALSE)
            for (c in items_shared) {
                t2_add[, c] <- t2[, c]
            }
            t_out <- rbind.data.frame(t1, t2_add, stringsAsFactors = FALSE)
        }
    }

    # Sort columns when it is specified
    if (sort_col) {
        t_out <- t_out[, order(names(t_out), decreasing = FALSE)]
    }
    t_out_cols <- names(t_out)

    # Add the column for sample names
    t_out$Sample <- rownames(t_out)
    t_out <- t_out[, c("Sample", t_out_cols)]
    rownames(t_out) <- NULL

    return(t_out)
}

t_out <- mergeTables(t1 = t1, t2 = t2, sort_col = args$sort_col)
print(paste("The new table has", nrow(t_out), "rows.", sep = " "))
write.table(x = t_out, file = args$out, sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = TRUE)
