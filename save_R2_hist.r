# create histogram of R^2 values (without plotting) and save to file

args<-commandArgs()

infile<-args[1]
n<-args[2]

ld.values.normal<-read.table(infile)
ld.u<-ld.values.normal[upper.tri(ld.values.normal)]
ld.u[!is.finite(ld.u)] <- NA
ld.u<-ld.u[complete.cases(ld.u)]
ld.hist<-hist(ld.u, breaks=seq(0,1,0.05), plot=F)

output<-data.frame(cbind(rep(n, length(ld.hist$mids)), hist$mids, hist$density))
names(output)<-c('n', 'mids', 'density')

f.name<-file.path('~/Desktop/UT_ddRADseq/Rsq', paste(n, '.csv', sep=''))

write.csv(output, f.name)

# for when the time comes, the plotting script...
#plot(ld.hist$mids, ld.hist$density, type='l', axes=F, xlab=expression('R'^'2'), ylab='', cex.axis=1.5, cex.lab=1.5)
#axis(side=1, labels=c(0,1), at=c(0,1), cex.axis=1.5, cex.lab=1.5)
#polygon(x=c(min(ld.hist$mids), ld.hist$mids, max(ld.hist$mids)), y=c(0, ld.hist$density, 0), col='blue')
