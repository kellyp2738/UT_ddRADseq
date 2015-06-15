# calculate significance of LD R^2 and make FDR correction

library('fdrtool')

options <- commandArgs(TRUE)

infile <- options[1]
#print(infile)

ld <- read.table(infile, header=TRUE)
#ld<-read.table('~/Desktop/UT_ddRADseq/plink_edit-2.ld', header=TRUE)

N <- 22 # number chromosomes in D. variabilis

ld$chisq <- ld$R2*N
ld$p <- dchisq(ld$chisq, df=1) # 2 alleles, df = 1
#print(length(ld$p))
ld.fdr <- fdrtool(ld$p, statistic=c('pvalue'), plot=FALSE)
ld$q <- ld.fdr$qval
sig.ld <- subset(ld, q < 0.05 & R2 > 0.8)
num.sig <- length(sig.ld[,1])

print(c(length(ld$p), num.sig)) # this should end up putting the number of sig. LD pairs in std. error output that Python can parse

## WHY DO I ONLY GET 84 PAIRWISE RESULTS FOR LD FROM PLINK??