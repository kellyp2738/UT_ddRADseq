## calculate significance of LD R^2 and make FDR correction (matrix input)

library(doMC)
library(bigtabulate)
library(bigmemory)
library(ff)

source("http://bioconductor.org/biocLite.R")
biocLite("qvalue")
library(qvalue)

# which correction to use? benjamini-hochberg adjusted p-value (as in p.adjust) or q-value?
#https://stat.ethz.ch/pipermail/bioconductor/attachments/20121219/00dc27b1/attachment.pl
# seems like q-values require large studies (>3000) and an expectation that a number of the values will be significant

options <- commandArgs(TRUE)

infile <- options[1]

N <- 22 # number chromosomes in D. variabilis

# simulate a small test dataset
v<-runif(15000)
test.data<-data.frame(matrix(v/10, nrow=3000, ncol=5))
test.chisq<-apply(test.data, 1:2, function(x) x*22) # convert to chi square test statistics
test.chisq.p<-apply(test.chisq, 1:2, function(x) pchisq(x, df=1, lower.tail=FALSE)) # calculate p values
test.chisq.p[!is.finite(test.chisq.p)] <- NA # convert inf to NA
qv.matrix<-data.frame(matrix(data=NA, nrow=dim(test.chisq.p)[1], ncol=dim(test.chisq.p)[2]))
fraction.sig<-c()
for(col in 1:length(test.chisq.p[1,])){
  #print(col)
  col.data<-test.chisq.p[,col]
  complete.data<-col.data[complete.cases(col.data)]
  #print(length(complete.data))
  #intermediate<-apply(test.chisq.p, 1:2, function(x) qvalue(col.data))
  qv.matrix[,col]<-qvalue(complete.data)$qvalues
  fraction.sig<-c(fraction.sig, length(qv.matrix[qv.matrix[,col]<0.028,col])/length(col.data))
}

# read data
ld.matrix <- read.table(infile)
ld.matrix <- read.table("~/Desktop/UT_ddRADseq/long_plink_pval.ld")

## OLD: TESTING BIG DATA FORMATS AND TIMING LOADING OPERATIONS ##

ld.matrix <- read.table(infile)
ld.matrix <- read.table("~/Desktop/UT_ddRADseq/long_plink_pval.ld")
ld.matrix <- read.table.ffdf(file="~/Desktop/UT_ddRADseq/long_plink_pval.ld")

ld.values<-read.table.ffdf(file='~/Desktop/UT_ddRADseq/long_plink.ld')
ld.p.vector <- dchisq((ld.values*N), 1)

ld.matrix <- read.big.matrix("~/Desktop/UT_ddRADseq/plink.ld", sep=" ", type='double')#, na.strings = "nan")

## Which way of reading in the raw LD values is fastest?

## matrix, read.table()
# Start the clock!
ptm <- proc.time()
ld.values.normal<-read.table(file='~/Desktop/UT_ddRADseq/plink.ld')
# Stop the clock
proc.time() - ptm
rm(ptm)

## long form, read.table()
# Start the clock!
ptm <- proc.time()
ld.list.normal<-read.table('~/Desktop/UT_ddRADseq/long_plink.ld')
# Stop the clock
proc.time() - ptm
rm(ptm)

## long form, ff package
# Start the clock!
ptm <- proc.time()
ld.values<-read.table.ffdf(file='~/Desktop/UT_ddRADseq/long_plink.ld')
# Stop the clock
proc.time() - ptm

## NORMAL MATRIX WITH READ.TABLE() IS THE FASTEST

ld.values.subset<-ld.values.normal[1:1000,900:800]


# calculate a chisq test statistic for each R2 in the matrix
ld.chisq<-apply(ld.values.subset, 1:2, function(x) x*22) # convert to chi square test statistics
ld.chisq.p<-apply(ld.chisq, 1:2, function(x) pchisq(x, df=1, lower.tail=FALSE)) # calculate p values
ld.chisq.p[!is.finite(ld.chisq.p)] <- NA # convert inf to NA

# calculate fdr for the chisq


# strip infinites to NA for easier removal later

test<-ld.chisq.p[complete.cases(ld.chisq.p[,1]), 1]
max(test)
head(test)

complete.q<-function(data.column){return(qvalue(data.column[complete.cases(data.column)]))}
fraction.still.sig<-c()
for(col in 1:length(chisq.p)){
  qv<-qvalue(chisq.p[,col])$qvalues
  
}
ld.qvals<-apply(ld.chisq.p, 2, function(x) qvalue(x)
q.matrix<-data.frame(matrix(unlist(ld.qvals), nrow=1000, ncol=101))
diffs<-data.frame(q.matrix[,1]-ld.chisq.p[,1]                
                
ld.qvals<-apply(ld.chisq.p, 1:2, function(x) p.adjust(x[complete.cases(x)], method=c('fdr')))


n.sig<-apply(ld.qvals, 2, )


# create a matrix of 10 rows x 2 columns
m <- matrix(c(1:10, 11:20), nrow = 10, ncol = 2)
# mean of the rows
apply(m, 1, mean)
[1]  6  7  8  9 10 11 12 13 14 15
# mean of the columns
apply(m, 2, function(x) sum(x[x<8])/2)
[1]  5.5 15.5
# divide all values by 2
apply(m, 1:2, function(x) x/2)




# reduce data size
ld.vector <- ld.matrix[upper.tri(ld.matrix)] # change matrix to a vector to reduce its size
rm(ld.matrix) # clear large matrix out

# parallelize the operations on the vector
registerDoMC(cores=4)
idx <- bigsplit(ld.vector)

split.idx<-function(vector){
  idx.list<-list()
  len<-length(vector)
  i=1
  j=1
  while(i <= len+1000){
    idx.list[[j]]=seq(i:i+999)
    i=i+1000
    j=j+1
  }
  return(idx.list)
}

test.split<-split.idx(ld.vector.normal)

acStart <- sapply ( acindices , function ( i )
  birthmonth ( x [i , c(’Year ’,’Month ’) , drop =
                    FALSE ]) )

ld.vector.normal<-ld.values.normal[upper.tri(ld.values.normal)]
ld.df<-data.frame(ld.vector.normal)
idx<-bigsplit(ld.df, ccols=c(1))
# calculate p-value
ld.p.value <- dchisq((ld.vector*N), 1)

# calculate q-value
ld.fdr <- p.adjust(ld.matrix, method=c('fdr'))

# create data frame for results and subset the significant SNPs
results.df <- data.frame(cbind(ld.vector, ld.fdr))
names(results.df) <- c('r2', 'q')

# get summary stats and print results
num.pairs <- length(results.df$r2)
complete <- results.df[complete.cases(results.df),]
num.sig <- length(subset(complete, r2 > 0.8 & q < 0.05)[,1])

print(c(num.pairs, num.sig))
