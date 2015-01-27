## tinkering with adegenet

library("ape")
library("genetics")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")

packageDescription("adegenet", fields="Version")

# for some reason the ploidy line (despite being "optional") is necessary
tick.snps<-read.snp("~/Dropbox/ddRADseq/final_tick_SNPs.snp")

# these also work
tick.snps2<-read.snp("~/Dropbox/ddRADseq/final_tick_SNPs_noPopInfo.snp")
file.show(system.file("files/exampleSnpDat.snp", package="adegenet"))
test.data<-read.snp(system.file("files/exampleSnpDat.snp", package="adegenet"))
my.test<-read.snp("~/Dropbox/adegenet_test_data_no_comments.snp")

ploidy(tick.snps)

glPlot(tick.snps)

pca.ticks<-glPca(tick.snps)

scatter(pca.ticks, posi="bottomright")
