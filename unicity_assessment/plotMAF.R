# Making a minor allele frequency (MAF) plot (in PNG format) for a single genotype call.
# 
# Example:
#   Rscript --vanilla plotMAF.R --vcf example.vcf
#
# Author: Yu Wan (\email{wanyuac@gmail.com})
# Dependency: packages optparse, vcfR and data.table
# Copyright 2017 Yu Wan
# Licensed under the Apache License, Version 2.0
# First edition: 29 August 2017, the latest edition: 31 August 2017

library(optparse)
library(vcfR, quietly = TRUE)
library(data.table, quietly = TRUE)

##### Functions and constants ##########
parseDP4 <- function(dp) {  # parse the DP4 field
    DP4 <- rbindlist(lapply(dp, function(x) as.data.frame(matrix(as.integer(strsplit(x = x, split = ",", fixed = FALSE)[[1]]), nrow = 1, ncol = 4))))
    names(DP4) <- c("REF_FOR", "REF_REV", "ALT_FOR", "ALT_REV")
    
    return(DP4)
}

calcMAF <- function(dp.ref, dp.alt, dp) {  # calculate minor allele frequency from each row (r) of the data frame aln
    dps <- c(dp.ref, dp.alt)
    dp.minor <- dps[which.min(dps)]
    maf <- dp.minor / dp  # dp = dp.ref + dp.alt
    
    return(maf)
}

determineMinorAllele <- function(ref, alt, dp.ref, dp.alt) {
    ma <- ifelse(dp.alt > dp.ref, ref, alt)
    
    return(ma)
}

assignColours <- function(ref, ma) {  # ref and ma: reference and minor alleles
    if (is.na(ma)) {  # Minor allele is the alternative allele and it is unknown.
        col <- "blue"
    } else if (ma == ref) {
        col <- "red"
    } else {
        col <- "blue"
    }
    
    return(col)
}

normaliseDP <- function(dp) {
    dp.min <- min(dp)
    dp.norm <- (dp - dp.min) / (max(dp) - dp.min)  # unity-based normalisation
    
    return(dp.norm)
}

# This function converts normalised values into their un-normalised forms.
# This is an inverse function of normaliseDP and it works for generating tick labels.
reverseNormalisation <- function(x.norm, xs) {
    x.max <- max(xs)
    x.min <- min(xs)
    x <- x.norm * (x.max - x.min) + x.min
    
    return(x)
}

X.DIV <- 100  # width of breaks dividing the X axis

###### Set up options ##########
options <- list(
    make_option("--vcf", type = "character", default = NULL, help = "Path to a VCF"),
    make_option("--refseq", type = "character", default = NULL, help = "The reference sequence whose mapping results will be plotted with"),
    make_option("--output_dir", type = "character", default = ".", help = "Minimal edge-2 depth [default %default]"),
    make_option("--width", type = "integer", default = 80, help = "Figure's width [default %default]"),
    make_option("--height", type = "integer", default = 55, help = "Figure's height [default %default]"),
    make_option("--res", type = "integer", default = 300, help = "Figure's resolution [default %default]"),
    make_option("--unit", type = "character", default = "mm", help = "Unit for figure's dimensions [default %default]"),
    make_option("--title", type = "character", default = NULL, help = "Title for the figure"),
    make_option("--annotation", type = "character", default = NULL, help = "Text to be printed at the top-right corner of the output figure")
)

opt <- parse_args(OptionParser(usage = "Rscript %prog [options]", option_list = options))  # read options

##### Preparation ##########
if (!dir.exists(opt$output_dir)) {  # returns TRUE when opt$output_dir = "."
    dir.create(opt$output_dir)
}

# parse the input VCF file
vcf <- read.vcfR(opt$vcf, verbose = FALSE)
aln <- cbind.data.frame(as.data.frame(vcf@fix[, c("CHROM", "POS", "REF", "ALT", "QUAL")], stringsAsFactors = FALSE),
                        parseDP4(extract.info(x = vcf, element = "DP4", as.numeric = FALSE)),
                        stringsAsFactors = FALSE)  # read alignment information
aln$POS <- as.integer(aln$POS)  # vcf@fix is a matrix of characters, even for the POS column
print("Completed parsing the VCF file.")

# Remove indels due to tandem single-nucleotide repeats, that is, to remove the row of an indel if bases
# after the first one in the REF equal the next reference bases of the same length.
indel <- extract.indels(x = vcf, return.indels = TRUE)@fix
aln$INDEL <- aln$REF %in% indel[, "REF"]  # TRUE or FALSE
print("Indel information is incorporated.")
indels.id <- which(aln$INDEL)
aln$rm <- rep(FALSE, times = nrow(aln))
for (i in indels.id) {
    ref <- aln$REF[i]
    tail.bases <- substr(ref, 2, nchar(ref))  # eg. GTT -> TT
    ref.next <- paste(aln$REF[(i + 1) : (i + nchar(tail.bases))], collapse = "")
    aln$rm[i] <- tail.bases == ref.next
}
aln <- aln[which(!aln$rm), 1 : (ncol(aln) -1)]  # exclude those rows for indels. aln[!aln$rm, -ncol(aln)] does not work on the Linux release of R/3.3.3.
print("Indels introduced by repeats are removed.")

# The BAM file SRST2's output contains alignments of reads against all reference sequences,
# but we only use the mapping result relating to a single reference here. As such, it is
# usually necessary to specify which CHROM is of interest.
if (!is.null(opt$refseq)) {
    aln <- subset(aln, CHROM == opt$refseq)
}

##### Calculate minor allele frequency for each base call (MAF) ##########
aln$DP_REF <- aln$REF_FOR + aln$REF_REV  # depth of high-quality reads supporting the presence of the refernence base
aln$DP_ALT <- aln$ALT_FOR + aln$ALT_REV  # depth of high-quality reads supporting the actual base call (may be identical to the reference base) in the current sample
aln$DP <- aln$DP_REF + aln$DP_ALT
aln$DP_NORM <- round(normaliseDP(aln$DP), digits = 4)
aln$MAF <- round(as.numeric(mapply(calcMAF, aln$DP_REF, aln$DP_ALT, aln$DP, USE.NAMES = FALSE)), digits = 4)
print("Minor allele frequencies are calculated.")

##### Draw a MAF plot ##########
# Whenever an allele is called, its opponent allele must be the minor one.
# Notice a VCF does not record the alternative allele when a reference allele is called for each base. Hence the MINOR_ALLELE becomes NA in this scenario,
# and the only non-NA values appearing in this column are reference allele (as minor alleles).
# sum(!is.na(aln$ALT) & (aln$DP_ALT <= aln$DP_REF)) = 0
# sum(is.na(aln$ALT) & (aln$DP_ALT > aln$DP_REF)) = 0
#aln$MINOR_ALLELE <- as.character(apply(aln[, c("REF", "ALT", "DP_REF", "DP_ALT")], 1, function(r) ifelse(r[["DP_ALT"]] > r[["DP_REF"]], r[["REF"]], r[["ALT"]])))  # This criterion is the same as the SRST2.
aln$MINOR_ALLELE <- mapply(determineMinorAllele, aln$REF, aln$ALT, aln$DP_REF, aln$DP_ALT, USE.NAMES = FALSE)

# Colour points blue when minor alleles are alternative alleles (hence reference alleles are called at these positions);
# colour points red when minor alleles are reference alleles (hence alternative alleles are called at these positions - genuine variants).
aln$COLOUR <- mapply(assignColours, aln$REF, aln$MINOR_ALLELE, USE.NAMES = FALSE)

# Save the result
write.table(x = aln, file = paste0(opt$output_dir, "/", opt$vcf, "__summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = ".")

# draw the plot
print("Generating MAF plots.")
x.max <- ceiling(max(aln$POS) / X.DIV) * X.DIV
aln <- aln[order(aln$POS, decreasing = FALSE), ]
#aln <- aln[order(aln$COLOUR, aln$POS, decreasing = FALSE), ]  # to draw red points after blue ones
aln.prt <- subset(aln, !is.na(MINOR_ALLELE) & MAF >= 0.3)
n <- nrow(aln.prt)
print(paste("There are", n, "minor alleles having frequencies >= 0.3.", sep = " "))
print.title <- !is.null(opt$title)

# for debugging: opt <- data.frame(output_dir = ".", vcf = "test", width = 80, height = 55, units = "mm", res = 300, annotation = "intI1_1.1", stringsAsFactors = FALSE)
png(filename = paste0(opt$output_dir, "/", opt$vcf, ".png"), width = opt$width, height = opt$height, units = opt$unit, res = opt$res)
if (print.title) {
    par(oma = rep(0.1, times = 4), mar = c(2, 1.8, 1, 1.8), mgp = c(1, 0.35, 0))
} else {
    par(oma = rep(0.1, times = 4), mar = c(2, 1.8, 0.1, 1.8), mgp = c(1, 0.35, 0))
}
plot(x = aln$POS, y = aln$DP_NORM / 2, type = "l", col = "grey50", lwd = 0.5,
     xlim = c(0, x.max), ylim = c(0, 0.54), axes = FALSE,
     xlab = "Base position", ylab = "MAF", cex.lab = 0.5)
points(x = aln$POS, y = aln$MAF, type = "h", col = aln$COLOUR)
axis(side = 1, at = seq(0, x.max, by = X.DIV), mgp = c(1, 0, 0), cex.axis = 0.4, tcl = -0.25)
axis(side = 2, at = seq(0, 0.5, by = 0.1), las = 2, cex.axis = 0.5, tcl = -0.25)
axis(side = 4, at = seq(0, 0.5, by = 0.05),
     labels = as.character(round(reverseNormalisation(x.norm = seq(0, 1, by = 0.1), xs = aln$DP), digits = 1)),
     las = 2, cex.axis = 0.5, tcl = -0.25)
rug(side = 1, x = seq(0, x.max, by = X.DIV / 2), ticksize = -0.025)
rug(side = 2, x = seq(0, 0.5, by = 0.05), ticksize = -0.025)
rug(side = 4, x = seq(0, 0.5, by = 0.025), ticksize = -0.025)
mtext("Read depth", side = 4, line = 1, cex = 0.5)  # The line argument controls the distance between the axis and its label.
#text(x = 0, y = 0.51, labels = paste0("Min. DP = ", as.character(min(aln$DP))),
#     cex = 0.4, adj = c(0, 0))  # It is redundant when the function reverseNormalisation is applied.
if (n > 0) {
    text(x = aln.prt$POS, y = aln.prt$MAF + 0.01, labels = aln.prt$MINOR_ALLELE, col = aln.prt$COLOUR,
         adj = c(0.5, 0), cex = 0.4)
}
if (print.title) {
    title(main = opt$title, cex.main = 0.5)
}
if (!is.null(opt$annotation)) {
    text(x = x.max - 300, y = 0.51, labels = opt$annotation, adj = c(0, 0), cex = 0.4)
}
dev.off()
