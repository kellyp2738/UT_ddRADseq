#!/sur/bin/env Rscript
# make a plot of the sensitivity analysis on subpopulation sample size

library(gplots)
library(plyr)
library(plotrix)
library(stringr)
library(scales)

args<-commandArgs()

infile<-args[1]

# read in the snp data file
data<-read.table(infile, header=FALSE)

data<-read.csv('~/Desktop/UT_ddRADseq/snp_bootstrap_final.txt', header=FALSE)
names(data)<-c('n', 'snps')

# create a boxplot that also displays the raw data (jittered on x-axis)
#png(file='~/Dropbox/ddRADseq/Final_Analysis/vcf_bootstrap_sampleSize.png', height=20, width=30, unit='cm', res=300)
pdf(file='~/Dropbox/ddRADseq/Final_Analysis/vcf_bootstrap_sampleSize.pdf', height=20/2.56, width=30/2.54)
par(mfrow=c(1,1), mar=c(5,5,2,1))
plot(x=NULL, y=NULL, xlim=c(0,75), ylim=c(4000,20000), axes=FALSE, xlab="", ylab="",
     main='75% Fragment Representation')
for(i in 1:((max(data$n)-min(data$n))+1)){
  n.data<-subset(data, n==unique(data$n)[i])
  points(jitter(rep(i, length(n.data$n)), amount=0.1), n.data$snps, pch=16, col=alpha('olivedrab4', 0.25))
}
boxplot(snps~n, data=data, axes=FALSE, add=TRUE, col=alpha('white', 0),
        ylim=c(round_any(min(data$snps), 1000, f=floor), 
               round_any(max(data$snps), 1000, f=ceiling)),
        ylab='SNPs Recovered (1000s)', xlab='Sample Size', cex.lab=1.5)

axis(side=1, cex.axis=1.5, at=seq(1, length(unique(data$n)), 5), labels=seq(10, 75, 5))
axis(side=2, las=1, at=seq(round_any(min(data$snps), 1000, f=floor), 
                           round_any(max(data$snps), 1000, f=floor)+1000, 1000),
     labels=seq(round_any(min(data$snps), 1000, f=floor), 
                round_any(max(data$snps), 1000, f=floor)+1000, 1000)/1000, cex.axis=1.5)
dev.off()

# determine average snp recover
ns<-c()
means<-c()
for(i in unique(data$n)){
  ns<-c(ns, i)
  means<-c(means, mean(subset(data, n==i)$snps))
}

lm.means<-lm(means~ns)
summary(lm.means)
lines(x=1:34, y=10854.95-(140.03*(1:34)), col='blue', lwd=3) # adds to scatterplot; don't run alone
points(ns, means, pch=20)
plot(residuals(lm.means))

# logistic regression on count data
log.model<-glm(snps~n, data=data, family=poisson)
summary(log.model)
plot(resid(log.model))
fit.model<-predict(log.model, type=c("response"), se.fit=TRUE)
plot(x=min(data$n):max(data$n), y=unique(fit.model$fit), type='l')
lines(x=min(data$n):max(data$n), y=unique(fit.model$fit)-(unique(fit.model$se.fit)*1.96), lty=2)
lines(x=min(data$n):max(data$n), y=unique(fit.model$fit)+(unique(fit.model$se.fit)*1.96), lty=2)
# the standard error for the parameter estimate in this model is incredibly tiny
# looking back at the distribution of the counts, the distributions seem normal
mean(data$snps[1:100])
sqrt(var(data$snps[1:100]))
# mean != variance --> logistic regression with a poisson link function does not seem to fit this data set
# linear regression on count data is okay when the counts are not near zero
# in our case, our counts are so far from zero that we needn't worry about negative predicted values
# and so... a linear regression model as an alternative

# a simple linear regression model
# not appropriate for count data, but the plot was reasonably difficult to make so its code are left intact
linear.model<-lm(snps~n, data=data)
summary(linear.model)
fit.with.ci<-predict(linear.model, data, level=0.95, interval="prediction", )
rsq.label<-bquote(R^2~.('=')~.(round(summary(linear.model)$r.squared, 2)))
pval.label<-bquote(p~.('=')~.(format(summary(linear.model)$coefficients[2,4], digits=2)))


#png(file='~/Desktop/test_plot.png', height=15, width=30, res=300, units='cm')
par(mar=c(5,5,2,2))
plot(x=data$n, y=fit.with.ci[,1], type="l", 
     ylim=c(round_any(min(data$snps), 1000, f=floor), 
            round_any(max(data$snps), 1000, f=ceiling)), 
     axes=FALSE, xlab='Sample Size', ylab='SNPs Recovered (1000s)',
     cex.lab=1.5, cex.axis=1.5)
#lines(x=data$n, y=fit.with.ci[,2], lty=2)
#lines(x=data$n, y=fit.with.ci[,3], lty=2)
lines(x=1:35, y=15240*exp(-0.14*(1:35))+6743)
points(data$n, data$snps, pch=16, col=alpha('red', 0.25))
axis(side=1, cex.axis=1.5)
axis(side=2, las=1, at=seq(round_any(min(data$snps), 1000, f=floor), 
                           round_any(max(data$snps), 1000, f=ceiling), 1000),
     labels=seq(round_any(min(data$snps), 1000, f=floor), 
                round_any(max(data$snps), 1000, f=ceiling), 1000)/1000, cex.axis=1.5)
corner.coords<-corner.label(x=1, y=1)
text(corner.coords$x-2.5, corner.coords$y-2527, labels=pval.label, adj=0)
text(corner.coords$x-2.5, corner.coords$y-3250, labels=rsq.label, adj=0)
legend(x=corner.coords$x-2.72, y=corner.coords$y, legend=c('Fit', '95% CI'), lty=c(1, 2), bty='n')
#dev.off()

plot(residuals(linear.model))
abline(h=0, col='red')
hist(residuals(linear.model))

exp.decay<-function(x, a){y=10000*exp(a*x)+2000
return(y)}
nls1<-nls(snps ~ (10000*exp(a*n)+2000), data, start=c(a=-0.1))#, weights=residuals(linear.model))
nls2<-nls(snps ~ (b*exp(a*n)+c), data, start=c(a=-0.1, b=10000, c=2000))
aov(nsl1, nls2)

### max-missing = 50%

data1<-read.csv('~/Desktop/UT_ddRADseq/snp_bootstrap_max_missing_0.5.txt', header=FALSE)
data2<-read.csv('~/Desktop/UT_ddRADseq/snp_bootstrap_max_missing_0.5_2.txt', header=FALSE)
data3<-read.csv('~/Desktop/UT_ddRADseq/snp_bootstrap_max_missing_0.5_3.txt', header=FALSE)
data4<-read.csv('~/Desktop/UT_ddRADseq/snp_bootstrap_max_missing_0.5_4.txt', header=FALSE)
full.data<-rbind(data1, data2, data3, data4)
names(full.data)<-c('n', 'snps')
boxplot(full.data$snps ~ full.data$n)

#png(file='~/Dropbox/ddRADseq/Final_Analysis/vcf_bootstrap_sampleSize_maxmissing_0.5.png', height=20, width=30, unit='cm', res=300)
png(file='~/Dropbox/ddRADseq/Final_Analysis/vcf_bootstrap_sampleSize_50_75_two_panels.png', height=40, width=30, unit='cm', res=300)

par(mfrow=c(2,1), mar=c(5,5,2,1))
plot(x=NULL, y=NULL, xlim=c(0,75), ylim=c(4000,20000), axes=FALSE, xlab="", ylab="",
     main='75% Fragment Representation', cex.main=1.5)
for(i in 1:((max(data$n)-min(data$n))+1)){
  n.data<-subset(data, n==unique(data$n)[i])
  points(jitter(rep(i, length(n.data$n)), amount=0.1), n.data$snps, pch=16, col=alpha('olivedrab4', 0.25))
}
boxplot(snps~n, data=data, axes=FALSE, add=TRUE, col=alpha('white', 0),
        ylim=c(round_any(min(data$snps), 1000, f=floor), 
               round_any(max(data$snps), 1000, f=ceiling)),
        ylab='SNPs Recovered (1000s)', xlab='Sample Size', cex.lab=1.5)

axis(side=1, cex.axis=1.5, at=seq(1, length(unique(data$n)), 5), labels=seq(10, 75, 5))
axis(side=2, las=1, at=seq(round_any(min(data$snps), 1000, f=floor), 
                           round_any(max(data$snps), 1000, f=floor)+1000, 1000),
     labels=seq(round_any(min(data$snps), 1000, f=floor), 
                round_any(max(data$snps), 1000, f=floor)+1000, 1000)/1000, cex.axis=1.5)

plot(x=NULL, y=NULL, xlim=c(0,75), ylim=c(8000, 28000), axes=FALSE, xlab="", ylab="",
     main='50% Fragment Representation', cex.main=1.5)
for(i in 1:((max(full.data$n)-min(full.data$n))+1)){
  full.n.data<-subset(full.data, n==unique(full.data$n)[i])
  points(jitter(rep(i, length(full.n.data$n)), amount=0.1), full.n.data$snps, pch=16, col=alpha('olivedrab3', 0.25))
}
boxplot(snps~n, data=full.data, axes=FALSE, add=TRUE, col=alpha('white', 0),
        ylim=c(round_any(min(full.data$snps), 1000, f=floor), 
               round_any(max(full.data$snps), 1000, f=ceiling)),
        ylab='SNPs Recovered (1000s)', xlab='Sample Size', cex.lab=1.5)

axis(side=1, cex.axis=1.5, at=seq(1, length(unique(full.data$n)), 5), labels=seq(10, 75, 5))
axis(side=2, las=1, at=seq(round_any(min(full.data$snps), 1000, f=floor), 
                           round_any(max(full.data$snps), 1000, f=floor)+1000, 1000),
     labels=seq(round_any(min(full.data$snps), 1000, f=floor), 
                round_any(max(full.data$snps), 1000, f=floor)+1000, 1000)/1000, cex.axis=1.5)

dev.off()