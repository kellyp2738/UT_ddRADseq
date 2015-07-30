## calculate significance of LD R^2 and make FDR correction (matrix input)

options <- commandArgs(TRUE)

infile <- options[1]

# read data
ld.matrix <- read.table(infile)
ld.vector <- read.table("~/Desktop/UT_ddRADseq/long_plink.ld", na.strings = "nan")

# reduce data size
ld.vector <- ld.matrix[upper.tri(ld.matrix)] # change matrix to a vector to reduce its size
rm(ld.matrix) # clear large matrix out
write.table(ld.vector, "~/Desktop/UT_ddRADseq/plink_long.ld" header=FALSE, quote=FALSE, row.names=FALSE)
