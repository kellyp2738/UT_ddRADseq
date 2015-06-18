## calculate significance of LD R^2 and make FDR correction (matrix input)

library(doMC)
library(bigtabulate)

options <- commandArgs(TRUE)

infile <- options[1]

N <- 22 # number chromosomes in D. variabilis

# read data
ld.matrix <- read.table(infile)
ld.matrix <- read.big.matrix("~/Desktop/UT_ddRADseq/plink.ld", sep=" ", type='double')#, na.strings = "nan")

# reduce data size
ld.vector <- ld.matrix[upper.tri(ld.matrix)] # change matrix to a vector to reduce its size
rm(ld.matrix) # clear large matrix out

# parallelize the operations on the vector
registerDoMC(cores=4)
idx <- bigsplit(ld.vector)
# calculate p-value
ld.p.value <- dchisq((ld.vector*N), 1)

# calculate q-value
ld.fdr <- p.adjust(ld.p.value, method=c('fdr'))

# create data frame for results and subset the significant SNPs
results.df <- data.frame(cbind(ld.vector, ld.fdr))
names(results.df) <- c('r2', 'q')

# get summary stats and print results
num.pairs <- length(results.df$r2)
complete <- results.df[complete.cases(results.df),]
num.sig <- length(subset(complete, r2 > 0.8 & q < 0.05)[,1])

print(c(num.pairs, num.sig))
