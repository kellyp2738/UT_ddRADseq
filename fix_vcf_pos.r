# fix the POS column of a VCF file

options <- commandArgs(TRUE)

infile <- options[1]
outname <- options[2]

#infile <- "~/Desktop/vcf_chrom_rename_3.vcf"

vcf <- read.table(infile, comment.char="#", header=TRUE) # header lines are commented out by default
vcf[,2] <- seq(1:length(vcf[,2]))
names(vcf)[1] <- "#CHROM"

write.table(vcf, outname,  quote=FALSE, row.names=FALSE, sep="\t")
