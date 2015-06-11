# calculate significance of LD R^2 and make FDR correction

library('fdrtool')

options <- commandArgs(TRUE)

infile <- options[1]

ld <- read.table(infile, header=TRUE)

N <- 22 # number chromosomes in D. variabilis

ld$chisq <- ld$R2*N
ld$p <- dchisq(ld$chisq, df=1) # 2 alleles, df = 1
ld$q <- fdrtool(ld$p, statistic=c('pvalue'))
sig.ld <- subset(ld, q < 0.05, R2 > 0.8)
num.sig <- length(sig.ld[,1])

print(num.sig) # this should end up putting the number of sig. LD pairs in std. error output that Python can parse