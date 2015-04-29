#!/sur/bin/env Rscript
# make a plot of the sensitivity analysis on subpopulation sample size

library(gplots)
library(plyr)
library(plotrix)

args<-commandArgs()

infile<-args[1]

# read in the snp data file
data<-read.table(infile)

data<-read.csv('~/Desktop/snp_bootstrap_test.txt', header=FALSE)
names(data)<-c('n', 'snps')
boxplot(snps~n, data=data)
plot(snps~n, data=data)

ns<-c()
means<-c()
for(i in unique(data$n)){
  ns<-c(ns, i)
  means<-c(means, mean(subset(data, n==i)$snps))
}
  
plot(ns, means)


aov(model)
confint(model)

model2<-aov(snps~as.factor(n), data=data)
summary(model2)

model3<-glm(snps~n, data=data, family=poisson)
summary(model3)
anova(model3)

ci.3<-confint(model3)
lower.ci<-function(x){
  y=x*ci.3[2,2]+coef(model3)[1]
  return(y)
}

upper.ci<-function(x){
  y=x*ci.3[2,1]+coef(model3)[1]
  return(y)
}

fit<-function(x){
  y=exp(x*coef(model3)[2])+coef(model3)[1]
}

plot(data$n, fit(data$n), type='l')
lines(data$n, lower.ci(data$n), lty=2)
lines(data$n, upper.ci(data$n), lty=2)

fit3.with.ci<-predict(model3, data, level=0.95, interval="confidence")

# a simple linear regression model
# not appropriate for count data, but the plot was reasonably difficult to make so its code are left intact
model<-lm(snps~n, data=data)
fit.with.ci<-predict(model, data, level=0.95, interval="confidence")
rsq.label<-bquote(R^2~.('=')~.(round(summary(model)$r.squared, 2)))

par(mar=c(5,5,2,2))
plot(x=data$n, y=fit.with.ci[,1], type="l", 
     ylim=c(round_any(min(fit.with.ci[,2]), 100, f=floor), 
            round_any(max(fit.with.ci[,3]), 100, f=floor)+1000), 
     axes=FALSE, xlab='Sample Size', ylab='SNPs Recovered (1000s)',
     cex.lab=1.5, cex.axis=1.5)
lines(x=data$n, y=fit.with.ci[,2], lty=2)
lines(x=data$n, y=fit.with.ci[,3], lty=2)
axis(side=1, cex.axis=1.5)
axis(side=2, las=1, at=seq(round_any(min(fit.with.ci[,2]), 1000, f=floor), 
                           round_any(max(fit.with.ci[,3]), 1000, f=floor)+1000, 1000),
     labels=seq(round_any(min(fit.with.ci[,2]), 1000, f=floor), 
                round_any(max(fit.with.ci[,3]), 1000, f=floor)+1000, 1000)/1000, cex.axis=1.5)
corner.label(label=rsq.label, x=1, y=1)


