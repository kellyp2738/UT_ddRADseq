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

data<-read.csv('~/Desktop/snp_bootstrap_final.txt', header=FALSE)
names(data)<-c('n', 'snps')

# create a boxplot that also displays the raw data (jittered on x-axis)
par(mar=c(5,5,2,2))
boxplot(snps~n, data=data, axes=FALSE,
        ylim=c(round_any(min(data$snps), 1000, f=floor), 
               round_any(max(data$snps), 1000, f=ceiling)),
        ylab='SNPs Recovered (1000s)', xlab='Sample Size', cex.lab=1.5)
for(i in 1:((max(data$n)-min(data$n))+1)){
  n.data<-subset(data, n==unique(data$n)[i])
  points(jitter(rep(i, length(n.data$n)), amount=0.1), n.data$snps, pch=16, col=alpha('red', 0.25))
}
axis(side=1, cex.axis=1.5, at=1:length(unique(data$n)), labels=c(min(data$n):max(data$n)))
axis(side=2, las=1, at=seq(round_any(min(data$snps), 1000, f=floor), 
                           round_any(max(data$snps), 1000, f=floor)+1000, 1000),
     labels=seq(round_any(min(data$snps), 1000, f=floor), 
                round_any(max(data$snps), 1000, f=floor)+1000, 1000)/1000, cex.axis=1.5)

# determine average snp recover
ns<-c()
means<-c()
for(i in unique(data$n)){
  ns<-c(ns, i)
  means<-c(means, mean(subset(data, n==i)$snps))
}

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


png(file='~/Desktop/test_plot.png', height=15, width=30, res=300, units='cm')
par(mar=c(5,5,2,2))
plot(x=data$n, y=fit.with.ci[,1], type="l", 
     ylim=c(round_any(min(data$snps), 1000, f=floor), 
            round_any(max(data$snps), 1000, f=ceiling)), 
     axes=FALSE, xlab='Sample Size', ylab='SNPs Recovered (1000s)',
     cex.lab=1.5, cex.axis=1.5)
lines(x=data$n, y=fit.with.ci[,2], lty=2)
lines(x=data$n, y=fit.with.ci[,3], lty=2)
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
dev.off()


