# create histogram of R^2 values (without plotting) and save to file

args<-commandArgs(trailingOnly=TRUE)

infile = args[1] #<-'plink.ld' # it will always be called this!

outfile2 = args[2] #= opts.output_R2
outfile3 = args[3] #=opts.output_hist
n = args[4] # sample size
fname = args[5] # filtered or unfiltered

#read in plink data
ld.values.normal<-read.table(infile)
#extract the upper triangle of the matrix, replace infinite values with NA, and retain only complete cases
ld.u<-ld.values.normal[upper.tri(ld.values.normal)]
ld.u[!is.finite(ld.u)] <- NA
ld.u<-ld.u[complete.cases(ld.u)]


#get distribution of non-inf R2 values
ld.hist<-hist(ld.u, breaks=seq(0,1,0.05), plot=F)
#prepare for saving
output<-data.frame(cbind(rep(n, length(ld.hist$mids)), ld.hist$mids, ld.hist$density))
names(output)<-c('n', 'mids', 'density')

#file name for representative histograms
hist.name<-file.path(paste(outfile3, '_', fname, sep=''))

# write by row
write.csv(output, file=hist.name, sep=',', row.names=F, append=T)

# separately record the fraction above certain thresholds
g.95<-length(ld.u[ld.u>0.95])/length(ld.u)
g.90<-length(ld.u[ld.u>0.9])/length(ld.u)
g.85<-length(ld.u[ld.u>0.85])/length(ld.u)
g.80<-length(ld.u[ld.u>0.8])/length(ld.u)
fs<-c(n, g.95, g.90, g.85, g.80)
#print(fs) # capture fraction sig in stdout/stderr to make one big file in python

#file name for r2 thresholds
cutoffs.name<-file.path(paste(outfile2, fname, sep=''))

write.csv(fs, cutoffs.name, sep=',', row.names=F, append=T)



# for when the time comes, the plotting script...
#plot(ld.hist$mids, ld.hist$density, type='l', axes=F, xlab=expression('R'^'2'), ylab='', cex.axis=1.5, cex.lab=1.5)
#axis(side=1, labels=c(0,1), at=c(0,1), cex.axis=1.5, cex.lab=1.5)
#polygon(x=c(min(ld.hist$mids), ld.hist$mids, max(ld.hist$mids)), y=c(0, ld.hist$density, 0), col='blue')
