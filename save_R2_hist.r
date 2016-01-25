# create histogram of R^2 values (without plotting) and save to file

args<-commandArgs(trailingOnly=TRUE)

infile<-'plink.ld' # it will always be called this!
parentDir<-args[1]
n<-args[2]
fname<-args[3] # which type of file (unfiltered/filtered) is being analyzed?
#print(n)

ld.values.normal<-read.table(infile)
ld.u<-ld.values.normal[upper.tri(ld.values.normal)]
ld.u[!is.finite(ld.u)] <- NA
ld.u<-ld.u[complete.cases(ld.u)]
ld.hist<-hist(ld.u, breaks=seq(0,1,0.05), plot=F)
output<-data.frame(cbind(rep(n, length(ld.hist$mids)), ld.hist$mids, ld.hist$density))
names(output)<-c('n', 'mids', 'density')

g.95<-length(ld.u[ld.u>0.95])/length(ld.u)
g.90<-length(ld.u[ld.u>0.9])/length(ld.u)
g.85<-length(ld.u[ld.u>0.85])/length(ld.u)
g.80<-length(ld.u[ld.u>0.8])/length(ld.u)
fs<-c(g.95, g.90, g.85, g.80)
print(fs) # capture fraction sig in stdout/stderr to make one big file in python

f.name<-file.path(parentDir paste(n, fname, 'Hist.csv', sep=''))
#g.name<-file.path('~/Desktop/UT_ddRADseq/Rsq', paste(n, 'fraction_sig.csv', sep=''))

write.csv(output, f.name, row.names=F, append=T)

# for when the time comes, the plotting script...
#plot(ld.hist$mids, ld.hist$density, type='l', axes=F, xlab=expression('R'^'2'), ylab='', cex.axis=1.5, cex.lab=1.5)
#axis(side=1, labels=c(0,1), at=c(0,1), cex.axis=1.5, cex.lab=1.5)
#polygon(x=c(min(ld.hist$mids), ld.hist$mids, max(ld.hist$mids)), y=c(0, ld.hist$density, 0), col='blue')
