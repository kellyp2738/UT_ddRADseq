## tinkering with adegenet

library("ape")
library("genetics")
# dependency for pegas; not on CRAN
install.packages("~/Desktop//Rcpp/", repos=NULL,
                 type="source")
library("pegas")
library("seqinr")
library("ggplot2")
library("RCo")

install.packages("~/Desktop//adegenet", repos=NULL,
                 type="source")

library("adegenet")
library("RColorBrewer")
library("scales")

packageDescription("adegenet", fields="Version")
# sample data in genlight format
file.show(system.file("files/exampleSnpDat.snp",package="adegenet"))

# for some reason the ploidy line (despite being "optional") is necessary
#tick.snps<-read.snp("~/Dropbox/ddRADseq/final_tick_SNPs_no_duplicates.snp")
#tick.snps.geo<-read.snp("~/Dropbox//ddRADseq/final_tick_SNPs_no_duplicates_geo_pop.snp")

tick.snps2<-read.snp('~/Desktop/D_variabilis_Pseudoref/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_HB_only_maxmissing0.75_MERGED_FOR_ADEGENET.finalSNPs.snp')

# check the structure of the data
ploidy(tick.snps)
glPlot(tick.snps)

# principle components analysis


## POPULATION STRUCTURE EXISTS AT A FINE SCALE (~3km) -- COMPARISON OF HB TO SRT (COMPARE TO STRUCTURE RESULTS)
site<-read.snp('~/Desktop/UT_ddRADseq/Both_Sites_Racc_Op_Pseudoref_minmeanDP20_minGQ25_maf0.05_FOR_ADEGENET.finalSNPs_SITE.snp')
pca.site<-glPca(site)
site.name<-site@pop
site.colors<-c('red', 'blue')
plot(pca.site$scores[,1], pca.site$scores[,2], col=site.colors[site.name], pch=19)
site.dapc<-dapc(site, pop=site.name, n.pca=33, n.da=1)
scatter(site.dapc)

## WITHIN HARRISON BAYOU, HOST SPECIES CORRESPONDS TO LOW-LEVEL POPULATION DIFFERENTIATION
tick.snps2<-read.snp('~/Desktop/D_variabilis_Pseudoref/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_HB_only_maxmissing0.75_MERGED_FOR_ADEGENET.finalSNPs.snp')
dapc2<-dapc(tick.snps2, pop=tick.snps2@pop, n.pca=27, n.da=1) 
scatter(dapc2)

## WITHIN HARRISION BAYOU, TICKS CLUSTER BY HOST INDIVIDUAL
# the population assignment in the data file doesn't change the outcome of the PCA
tick.snps<-read.snp('~/Desktop/D_variabilis_Pseudoref/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_HB_only_maxmissing0.75_MERGED_FOR_ADEGENET.finalSNPs_HOST.snp')
pca.ticks<-glPca(tick.snps)
scatter(pca.ticks, posi="bottomleft")
loadingplot(pca.ticks)
color.host=c(brewer.pal(9, 'Reds')[8:9], brewer.pal(9, 'Blues')[5:9])
hosts<-tick.snps@pop
plot(pca.ticks$scores[,1], pca.ticks$scores[,2],
     col=alpha(color.host[hosts], 0.5), pch=19, xlab='PCA 1', ylab='PCA 2')
dapc1<-dapc(tick.snps, pop=hosts, n.pca=27, n.da=6)
pdf(file='~/Dropbox/ddRADseq/dapc_by_individual.pdf', height=10/2.54, width=10/2.54)
scatter(dapc1)
dev.off()

host.corr<-as.data.frame(cor(hosts, pca.ticks$scores[,1]))

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
plot(tick.nj, type='phylogram', show.tip.label=FALSE, use.edge.length=FALSE)
# color by collection site
tiplabels(pch=16, col=color.host[hosts])
# color by host species
tiplabels(pch=16, col=color[hosts])
# color by host individual
tiplabels(pch=16, col=indv.colors[host.indv])
# it appears that some ticks on the same individual host are related

#tick.blocks<-seploc(tick.snps, n.block=10, parallel=FALSE)

# the clustering function works best if you give it
# the pre-calculated PCA
# there doesn't seem to be any strong support for clustering


tick.clusters<-find.clusters(tick.snps, glPca=pca.ticks)
dapc2<-dapc(tick.clusters)
scatter(tick.clusters)
tick.clusters.2<-find.clusters(tick.snps, 
                               max.n.clust=length(unique(host.indv)),
                               glPca=pca.ticks)

dapc3<-dapc(tick.snps2, pop=tick.snps2@pop)
scatter(dapc3)

