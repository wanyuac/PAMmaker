# This script generates a diagnostic scatter plot for inferring unicity of
# alleles in bacterial samples, namely, whether there is only a single allele in
# a sample.
#
# Inputs:
#   1. a tab-delimited file in which the first and third columns provide
# sample names and genomic mapping depths, respectively. It can be the file
# RepStats.txt in RedDog outputs (so I assume the third column stores the depth
# information).
#   2. a merged tab-delimited score file containing columns: Sample, Allele,
# Avg_depth and maxMAF with allele names renamed based on sequence clustering.
# It can be produced using maxMAFstats.R. It is not required that this table
# contain all samples in RepStats.txt.
#
# Outputs:
#   1. a scatter plot comparing maxMAF against putative copy numbers;
#   2. a tab-delimited table for data underlying the plot.
#
# Other parameters:
#   samples: a file listing sample names to be included. It is optional.
#   basename: basename and path for an output PNG image and a TSV file
#   width and height: width and height in millimeter for the image
#
# Command line:
#   Rscript --vanilla --genomicDepth [RepStats.txt] --scores [renamedScores.tsv]
# --samples [ingroup.txt] --basename unicity/polyploidy_diag --width 160
# --height 120
#
# Author: Yu Wan <wanyuac@gmail.com>
# Copyright 2018 Yu Wan
# Licensed under the Apache License, Version 2.0
# First edition: 9 Mar 2018; the latest edition: 17 July 2018

# Read arguments from the command line ###############
library(optparse)

options <- list(
    make_option("--genomicDepth", type = "character",
                help = "A tab-delimited table for sample names and genomic mapping depths"),
    make_option("--scores", type = "character",
                help = "A tab-delimited file containing merged SRST2 scores"),
    make_option("--samples", type = "character", default = "",
                help = "(optional) A newline-delimited file listing samples to be included"),
    make_option("--basename", type = "character", default = "unicity_diag",
                help = "Basename for output files [default %default]"),
    make_option("--width", type = "integer", default = 160,
                help = "Width (mm) of the image [default %default]"),
    make_option("--height", type = "integer", default = 160,
                help = "Height (mm) of the output image [default %default]"),
    make_option("--point_size", type = "double", default = 0.8,
                help = "Point size [default %default]")
)

# Define functions that can be run separately ###############
mergeTables <- function(gd, scores, ingroup = NA) {
    # This function produces the table for plotting.
    scores <- scores[, c("Sample", "Allele", "Avg_depth", "maxMAF")]
    names(scores)[3] <- "Allelic_depth"

    # filter samples when the ingroup argument is specified as a character vector
    if (! is.na(ingroup[1])) {
        scores <- subset(scores, Sample %in% ingroup)
    }

    # match read depths
    gd <- gd[, c(1, 3)]
    names(gd) <- c("Sample", "Genomic_depth")

	# The following step does not require "scores" to contain all samples in "gd".
    scores$Genomic_depth <- gd$Genomic_depth[match(scores$Sample, gd$Sample)]

	# Calculate the ratio of allelic depth over the genomic read depth
    scores$AG_ratio <- round(scores$Allelic_depth / scores$Genomic_depth,
                             digits = 4)  # measures potential of being polyploid
    scores <- scores[, c("Sample", "Allele", "Allelic_depth", "Genomic_depth",
                         "AG_ratio", "maxMAF")]
    return(scores)
}

# Main program ###############
# Import data
opt.parser <- OptionParser(usage = "Rscript %prog [options]", option_list = options)
opts <- parse_args(opt.parser)
gd <- read.delim(file = opts$genomicDepth, check.names = FALSE, header = TRUE,
                 sep = "\t", stringsAsFactors = FALSE)
scores <- read.delim(file = opts$scores, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
if (opts$samples != "") {
    ingroup <- read.delim(file = opts$samples, header = FALSE, sep = "\n",
                          stringsAsFactors = FALSE)[, 1]
} else {
    ingroup <- NA
}

# Produce and save the table to plot with
tab <- mergeTables(gd, scores, ingroup)
fn <- paste(opts$basename, "tsv", sep = ".")
print(paste0("Printing results into ", fn, "."))
write.table(x = tab, file = fn, quote = FALSE, sep = "\t", row.names = FALSE,
            col.names = TRUE)

# Generate the figure
library(grid)
library(gridExtra)
library(ggplot2)

grid_colour <- "#E0E0E0"  # colour for the major grid

# main panel: a scatter plot
# Must put theme_bw() before theme() to get a correct appearance of axes and labels.
# Set the scale and breaks for the Y axis, otherwise the ticks do not match between
# the left and main panels.
p <- ggplot(data = tab, mapping = aes(x = AG_ratio, y = maxMAF)) +
    geom_point(alpha = 0.80, colour = "red", size = opts$point_size) +
    geom_vline(xintercept = 1, colour = "grey25", size = 0.5) +
    annotate("text", x = 1.25, y = 0, label = "1", size = 4) +
    scale_x_continuous(trans = "log2") +
    scale_y_continuous(breaks = seq(0, 0.5, by = 0.1), limits = c(0, 0.5)) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line.y = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
          panel.grid.major = element_line(colour = grid_colour))

# left panel: a boxplot of maxMAF
pl <- ggplot(data = tab, mapping = aes(x = "", y = maxMAF)) +
    geom_boxplot(outlier.size = opts$point_size) +
    labs(x = NULL, y = "maxMAF") + theme_bw() +
    scale_y_continuous(breaks = seq(0, 0.5, by = 0.1), limits = c(0, 0.5)) +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 10), axis.title.y = element_text(size = 12),
          axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.line.y = element_line(colour = "black"),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_line(colour = grid_colour))

# bottom panel: a boxplot of GA_ratios
# The Y ticks must not be "expression(log[2](AG_ratio))" although the grid is
# drawn at a log scale.
pb <- ggplot(data = tab, mapping = aes(x = "", y = AG_ratio)) +
    geom_boxplot(outlier.size = opts$point_size) +
    labs(x = NULL, y = "Allele-genome depth ratio") +
    scale_y_continuous(trans = "log2") + coord_flip() + theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12),
          axis.text.y = element_blank(), axis.title.y = element_blank(),
          axis.line.x = element_line(colour = "black"),
          axis.ticks.y = element_blank(),
          panel.grid.major = element_line(colour = grid_colour))

fn <- paste(opts$basename, "png", sep = ".")
print(paste0("Making the figure ", fn, "."))

png(filename = fn, width = opts[["width"]], height = opts[["height"]], units = "mm", res = 300)
grid.arrange(pl, p, textGrob(""), pb, ncol = 2, widths = c(1, 4), heights = c(4, 1))
dev.off()

print("Done.")
