# Making minor allele frequency (MAF) plots (in PNG format) for genotype calls
#
# This script runs plotMAF.R for a list of samples.
#
# Author: Yu Wan (\email{wanyuac@gmail.com})
# Dependency: packages optparse, vcfR and data.table
# Copyright 2017 Yu Wan
# Licensed under the Apache License, Version 2.0
# First edition and the latest edition: 29 August 2017

library(optparse)

###### Set up options ##########
options <- list(
    make_option("--bams", type = "character", default = NULL, help = "A list of BAM files to be processed (must has a header row)"),
    make_option("--db", type = "character", default = NULL, help = "A FASTA file accompanying faidx-indexed reference sequences. It may contain multiple sequences."),
    make_option("--refseq", type = "character", default = NULL, help = "The reference sequence whose mapping results will be plotted with"),
    make_option("--script", type = "character", default = "plotMAF.R", help = "Path to plotMAF.R"),
    make_option("--Rscript", type = "character", default = "Rscript", help = "Path to the Rscript program"),
    make_option("--samtools", type = "character", default = "samtools", help = "Path to SAMtools [default %default]"),
    make_option("--bcftools", type = "character", default = "bcftools", help = "Path to BCFtools [default %default]"),
    make_option("--baseQ", type = "integer", default = 20, help = "The -Q argument for the command samtools mpileup [default %default]"),
    make_option("--output_dir", type = "character", default = ".", help = "Minimal edge-2 depth [default %default]"),
    make_option("--width", type = "integer", default = 80, help = "Figure's width [default %default]"),
    make_option("--height", type = "integer", default = 55, help = "Figure's height [default %default]"),
    make_option("--res", type = "integer", default = 300, help = "Figure's resolution [default %default]"),
    make_option("--unit", type = "character", default = "mm", help = "Unit for figure's dimensions [default %default]"),
    make_option("--override", type = "logical", action = "store_true", default = FALSE, help = "Whether to override previous outputs (incl. VCF files) [default %default]"),
    make_option("--debug", type = "logical", action = "store_true", default = FALSE, help = "Flag this option to avoid actual execution of commands.")
)

opt <- parse_args(OptionParser(usage = "Rscript %prog [options]", option_list = options))  # read options

##### Run plotMAF.R based on each BAM file ##########
bams <- read.delim(opt$bams, stringsAsFactors = FALSE)
fields <- names(bams)
n <- nrow(bams)
multi.column <- ncol(bams) > 1  # whether bams has multiple columns
run <- !opt$debug

if (!dir.exists(opt$output_dir) & run) {  # returns TRUE when opt$output_dir = "."
    dir.create(opt$output_dir)
}

print(paste0(n, " BAM files to be processed."))
for (i in 1 : n) {
    print(paste0("Processing the ", i, "th BAM."))
    
    # retrieve information from the input data frame
    rw <- bams[i, ]  # take a row from the data frame. Notice rw becomes a string (not a data frame) when ncol(bams) = 1.
    tit <- NULL
    annot <- NULL
    if (multi.column) {
        bam <- rw[1, 1]
        if ("Title" %in% fields) {
            tit <- rw$Title[1]
        }
        if ("Annotation" %in% fields) {
            annot <- rw$Annotation[1]
        }
    } else {
        bam <- rw
    }
    
    # make a MAF plot
    bname <- basename(bam)
    vcf <- paste0(opt$output_dir, "/", bname, ".vcf")
    img <- paste0(opt$output_dir, "/", vcf, ".png")
    if (file.exists(vcf) & file.exists(img) & !opt$override) {
        print(paste("Skip processing", bam, "as its outputs have been generated.", sep = " "))
        next
    } else {
        print(paste0("Generating a MAF plot from ", bam, "."))
        cmd <- paste(opt$samtools, "mpileup -L 1000 -u -f", opt$db, "-Q", opt$baseQ, "-q 1 -B", bam, "|",
                     opt$bcftools, "view -cg - >", vcf, sep = " ")
        print(paste0("Running command: ", cmd))
        if (run) {
            system(cmd, wait = TRUE)  # converting BAM into VCF
        }
        cmd <- paste(opt$Rscript, "--vanilla", opt$script, "--vcf", vcf, "--width", opt$width,
                     "--height", opt$height, "--res", opt$res, "--unit", opt$unit, sep = " ")  # Herein, vcf already contains the output directory. Do not set it for plotMAF.R.
        if (!is.null(opt$refseq)) {
            cmd <- paste(cmd, "--seq", opt$refseq, sep = " ")
        }
        if (!is.null(tit)) {
            cmd <- paste(cmd, "--title", tit, sep = " ")
        }
        if (!is.null(annot)) {
            cmd <- paste(cmd, "--annotation", annot, sep = " ")
        }
        print(paste0("Running command: ", cmd))
        if (run) {
            system(cmd, wait = TRUE)
        }
    }
}
