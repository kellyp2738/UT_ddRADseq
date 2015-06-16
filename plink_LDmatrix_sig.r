## calculate significance of LD R^2 and make FDR correction (matrix input)

options <- commandArgs(TRUE)

infile <- options[1]

N <- 22 # number chromosomes in D. variabilis

ld.matrix <- read.table(infile)
ld.matrix <- read.table("~/Desktop/UT_ddRADseq/plink.ld", na.strings = "nan")
matrix.size <- dim(ld.matrix[1])

# calculate chi sq. p-values
chisq.mat <- matrix(data=NA, nrow=matrix.size, ncol=matrix.size)
chisq.mat[upper.tri(chisq.mat)] <- dchisq((ld.matrix[upper.tri(ld.matrix)]*N), 1) # 2 alleles, df = 1

# FDR using p.adjust() in {base}
ld.fdr <- p.adjust(chisq.mat[upper.tri(chisq.mat)], method=c('fdr'))

results.df <- data.frame(cbind(ld.matrix[upper.tri(ld.matrix)], ld.fdr))
names(results.df) <- c('r2', 'q')

complete <- results.df[complete.cases(results.df),]
num.pairs <- length(results.df$r2)
num.sig <- subset(complete, r2 > 0.8, q < 0.05)

# does not play nice with NAs: 
# library('fdrtool')
# ld.fdr <- fdrtool(chisq.mat[upper.tri(chisq.mat)], statistic=c('pvalue'), plot=FALSE)
# ld$q <- ld.fdr$qval
# sig.ld <- subset(ld, q < 0.05 & R2 > 0.8)
# num.sig <- length(sig.ld[,1])
