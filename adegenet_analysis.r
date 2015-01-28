## tinkering with adegenet

library("ape")
library("genetics")
# dependency for pegas; not on CRAN
install.packages("~/Desktop//Rcpp/", repos=NULL,
                 type="source")
library("pegas")
library("seqinr")
library("ggplot2")

install.packages("~/Desktop//adegenet", repos=NULL,
                 type="source")

library("adegenet")

packageDescription("adegenet", fields="Version")
# sample data in genlight format
file.show(system.file("files/exampleSnpDat.snp",package="adegenet"))

# for some reason the ploidy line (despite being "optional") is necessary
tick.snps<-read.snp("~/Dropbox/ddRADseq/final_tick_SNPs_no_duplicates.snp")
tick.snps.geo<-read.snp("~/Dropbox//ddRADseq/final_tick_SNPs_no_duplicates_geo_pop.snp")

# check the structure of the data
ploidy(tick.snps)
glPlot(tick.snps)

# principle components analysis - population defined as host
pca.ticks<-glPca(tick.snps)
scatter(pca.ticks, posi="bottomright")
color.host=c('green', 'purple')
hosts<-tick.snps@pop
plot(pca.ticks$scores[,1], pca.ticks$scores[,2],
     col=color.host[hosts], pch=16)
host.corr<-as.data.frame(cor(hosts, pca.ticks$scores[,1]))

# principle components analysis - population defined as location
pca.geo<-glPca(tick.snps.geo)
pca.geo
scatter(pca.geo, posi="topright")
# with better formatting and color coding of locations
color=c('red','blue')
locations<-tick.snps.geo@pop
plot(pca.geo$scores[,1], pca.geo$scores[,2], 
     col=color[locations], pch=16)
points(pca.ticks$scores[,1], pca.ticks$scores[,2], col='purple')
# okay, so we're clear (and this should have been obvious)...
# the categorical data on population or host are in no way factored
# in to the principle components analysis

# distribution of allele frequencies (2nd allele)
allele.freq<-glMean(tick.snps)
hist(allele.freq)

# converting the data to a matrix allows some other operations
tick.snp.matrix<-as.matrix(tick.snps)

# host individuals
grouping.vars<-read.csv('~/Dropbox//ddRADseq/grouping_variables_no_duplicates_missingdata.csv')
host.indv<-grouping.vars$host.ID
indv.colors<-rainbow(length(unique(host.indv)))

# a quick (but not terribly rigorous) way to look at the tick relatedness
# is to make a neighbor joining tree from the SNP data
# distance matrix
tick.dist<-dist(tick.snp.matrix)
# tree
tick.nj<-nj(tick.dist)
plot(tick.nj, type='fan')
# color by collection site
tiplabels(pch=16, col=color[locations])
# color by host species
tiplabels(pch=16, col=color[hosts])
# color by host individual
tiplabels(pch=16, col=indv.colors[host.indv])
# it appears that some ticks on the same individual host are related

#tick.blocks<-seploc(tick.snps, n.block=10, parallel=FALSE)

# the clustering function works best if you give it
# the pre-calculated PCA
# there doesn't seem to be any strong support for clustering
# delta(BIC) from 1 to 2 clusters = +3 (1 cluster = better)
tick.clusters<-find.clusters(tick.snps, glPca=pca.ticks)
tick.clusters.2<-find.clusters(tick.snps, 
                               max.n.clust=length(unique(host.indv)),
                               glPca=pca.ticks)
