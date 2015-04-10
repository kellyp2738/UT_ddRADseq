#############################################################################################################
##
##   ANALYSIS OF POPULATION STRUCTURE, TICK RADSEQ DATA
##
#############################################################################################################

## INSTALL & LOAD PACKAGES
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

## LOAD SUMMARY DATA FILE
extra.data<-read.csv('~/Desktop/UT_ddRADseq/ddRAD_FinalLibrary_SampleInfo_Full.csv')

## CHECK ADEGENET FUNCTIONALITY AND DATA FORMAT
packageDescription("adegenet", fields="Version")
# sample data in genlight format
file.show(system.file("files/exampleSnpDat.snp",package="adegenet"))

# for some reason the ploidy line (despite being "optional") is necessary
#tick.snps<-read.snp("~/Dropbox/ddRADseq/final_tick_SNPs_no_duplicates.snp")
#tick.snps.geo<-read.snp("~/Dropbox//ddRADseq/final_tick_SNPs_no_duplicates_geo_pop.snp")

#################################################################################################################
##
##   POPULATION STRUCTURE BY SITE -- HARRISON BAYOU VS STARR RANCH TRAIL
##
#################################################################################################################

## POPULATION STRUCTURE EXISTS AT A FINE SCALE (~3km) -- COMPARISON OF HB TO SRT (COMPARE TO STRUCTURE RESULTS)

## RELATED TICKS NOT REMOVED
# read in data
site<-read.snp('~/Desktop/D_variabilis_Pseudoref/Final_Analysis/Structure_by_Site/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_BOTH_SITES_host_filtered_only_maxmissing0.75_MERGED_FOR_ADEGENET.finalSNPs_site.snp')
# look at principal coordinates
pca.site<-glPca(site)
site.name<-site@pop
site.colors<-c('gray19', 'gray70')
plot(pca.site$scores[,1], pca.site$scores[,2], col=site.colors[site.name], pch=19)
# discriminant analysis of principal components to detect clusters
site.dapc<-dapc(site, pop=site.name, n.pca=33, n.da=1) # save as many PCs as possible
pdf(file='~/Dropbox/ddRADseq/Final_Plots_April_2015/dapc_by_site.pdf', width=10/2.54, height=10/2.54)
par(mar=c(5,5,2,2))
scatter(site.dapc, col=c('gray19', 'gray70'), axes=FALSE)
axis(side=2, las=1)
legend("topright", legend=c('HB', 'SRT'), fill=alpha(site.colors, 0.5), bty='n')
dev.off()
# dapc posterior cluster assignments (essentially a structure plot)
post.data.site<-data.frame(site.dapc$posterior)
post.df.site<-data.frame(cbind(row.names(post.data.site), post.data.site$HB, post.data.site$SRT)) #coerces to char b/c of names
names(post.df.site)<-c('ID', 'HB', 'SRT')
post.df.site.c<-merge(post.df.site, extra.data, by.x='ID', by.y='combo.label') # add site data to posteriors data frame
post.df.site.trimmed<-cbind(post.df.site.c[,1:3], post.df.site.c$coll.site)
post.df.site.trimmed.sort<-post.df.site.trimmed[order(post.df.site.trimmed$HB),]
colors.site<-c('gray19', 'gray70')
pdf(file='~/Dropbox/ddRADseq/Final_Plots_April_2015/dapc_by_site_posterior.pdf', width=25/2.54, height=10/2.54)
par(mar=c(6,4,2,8), xpd=TRUE)
barplot(t(post.df.site.trimmed.sort[2:3]), col=colors.site, axes=FALSE, ylab='Probability of Membership',
        xlim=c(0, length(post.df.site.trimmed.sort[,1])), space=0, col.axis='white') #for some reason #s still show on x-axis... make their color white
mtext(side=1, las=2, at=seq(0.5, length(post.df.site.trimmed.sort[,1])-0.5), cex=0.5, col=colors.site[post.df.site.trimmed.sort[,4]],
      text=gsub(pattern="_.*", replacement="", post.df.site.trimmed.sort[,1], perl=TRUE))
axis(side=2, las=1)
legend(x=101, y=0.55, legend=c('SRT', 'HB'), fill=c('gray19', 'gray70'), title='Cluster \nAssignment', bty="n")
dev.off()

#################################################################################################################
##
##   POPULATION WITHIN HARRISON BAYOU
##
#################################################################################################################


## WITHIN HARRISON BAYOU, HOST SPECIES CORRESPONDS TO LOW-LEVEL POPULATION DIFFERENTIATION

# read in data
tick.snps2<-read.snp('~/Desktop/D_variabilis_Pseudoref/Final_Analysis/Structure_by_Host_HB/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_HB_only_maxmissing0.75_MERGED_not_related_FOR_ADEGENET.finalSNPs.snp')
# discriminant analysis of principal components to detect clusters
dapc2<-dapc(tick.snps2, pop=tick.snps2@pop, n.pca=21, n.da=1) 
pdf(file='~/Dropbox/ddRADseq/Final_Plots_April_2015/dapc_by_host.pdf', height=10/2.54, width=10/2.54)
par(mar=c(5,5,2,2))
scatter(dapc2)
axis(side=2, las=1)
legend("topleft", legend=c('Raccoon', 'Opossum'), fill=c('red', 'blue'), bty='n')
dev.off()
# dapc posterior cluster assignments (essentially a structure plot)
post.data<-data.frame(dapc2$posterior)
post.df<-data.frame(cbind(row.names(post.data), as.numeric(post.data$O), as.numeric(post.data$R)))
names(post.df)<-c('tick', 'O', 'R')
post.df.host.type<-merge(post.df, extra.data, by.x='tick', by.y='combo.label')
post.df.ht.trimmed<-cbind(post.df.host.type[,1:3], post.df.host.type$host.species)
post.df.sort<-post.df.ht.trimmed[order(post.df.ht.trimmed[,2]),] # sort character data w/names in the same way
label.cols<-c('blue','red')
pdf(file='~/Dropbox/ddRADseq/Final_Plots_April_2015/dapc_by_host_posterior.pdf', height=10/2.54, width=25/2.54)
par(mar=c(6,4,2,8), xpd=TRUE)
barplot(t(post.df.sort[2:3]), space=0, col.axis='white', ylab='Probability of Membership',
        las=2, cex.names=0.5, xlim=c(0,length(post.df.sort[,1])), axes=FALSE, col=c('blue', 'red'))
mtext(side=1, las=2, at=seq(0.5, length(post.df.sort[,1])-0.5, by=1), text=gsub(pattern="_.*", replacement="", post.df.sort[,1], perl=TRUE),
     line=1, col=label.cols[post.df.sort[,4]], cex=0.6)
axis(side=2, las=1)
legend(x=64, y=0.55, legend=c('Raccoon', 'Opossum'), fill=c('red', 'blue'), title='Cluster \nAssignment', bty="n")
dev.off()

## WITHIN HARRISION BAYOU, TICKS DO NOT CLUSTER BY HOST INDIVIDUAL

# read in data
# switch to not-related data!
#tick.snps4<-read.snp('~/Desktop/D_variabilis_Pseudoref/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_HB_only_maxmissing0.75_MERGED_FOR_ADEGENET.finalSNPs_HOST.snp')
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

#################################################################################################################
##
##   POPULATION WITHIN STARR RANCH TRAIL
##
#################################################################################################################


## WITHIN STARR RANCH TRAIL, HOST SPECIES CORRESPONDS TO LOW-LEVEL POPULATION DIFFERENTIATION?

## read in data
tick.snps3<-read.snp('~/Desktop/D_variabilis_Pseudoref/Final_Analysis/Structure_by_Host_SRT/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_SRT_only_maxmissing0.75_MERGED_not_related_FOR_ADEGENET.finalSNPs.snp')
pca.srt<-glPca(tick.snps3)
scatter(pca.srt)

pdf(file='~/Dropbox/ddRADseq/Final_Plots_April_2015/dapc_by_host_SRT.pdf', height=10/2.54, width=10/2.54)
par(mar=c(5,5,2,2))
dapc3<-dapc(tick.snps3, pop=tick.snps3@pop, n.pca=5, n.da=1) 
scatter(dapc3)
axis(side=2, las=1)
legend("topright", legend=c('Raccoon', 'Opossum'), fill=c('red', 'blue'), bty='n')
dev.off()

post.data3<-data.frame(dapc3$posterior)
post.df3<-data.frame(cbind(row.names(post.data3), as.numeric(post.data3$O), as.numeric(post.data3$R))) #coerces to char b/c of names
names(post.df3)<-c('tick', 'O', 'R')
post.df3.sort<-post.df3[order(post.df3[,3], decreasing=TRUE),] # sort character data w/names in the same way
label.cols3<-c(rep('red', 13), rep('blue', 2))
pdf(file='~/Dropbox/ddRADseq/Final_Plots_April_2015/dapc_by_host_posterior_SRT.pdf', height=10/2.54, width=25/2.54)
par(mar=c(6,4,2,8), xpd=TRUE)
barplot(t(post.df3.sort[2:3]), space=0, col.axis='white', ylab='Probability of Membership',
        las=2, cex.names=0.5, xlim=c(0,length(post.df3.sort[,1])), axes=FALSE, col=c('blue', 'red'))
mtext(side=1, las=2, at=seq(0.5, length(post.df3.sort[,1])-0.5, by=1), text=gsub(pattern="_.*", replacement="", post.df3.sort[,1], perl=TRUE),
      line=1, col=label.cols3, cex=0.6)
axis(side=2, las=1)
legend(x=16, y=0.55, legend=c('Raccoon', 'Opossum'), fill=c('red', 'blue'), title='Cluster \nAssignment', bty="n")
dev.off()


## WITHIN STARR RANCH TRAIL, TICKS DO NOT CLUSTER BY HOST INDIVIDUAL
# the population assignment in the data file doesn't change the outcome of the PCA
tick.snps5<-read.snp('~/Desktop/D_variabilis_Pseudoref/Final_Analysis/Structure_by_Host_SRT/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_SRT_only_maxmissing0.75_MERGED_not_related_FOR_ADEGENET.finalSNPs_HOST.snp')
pca.ticks5<-glPca(tick.snps5)
scatter(pca.ticks5, posi="bottomleft")
loadingplot(pca.ticks5)
color.host=c(brewer.pal(9, 'Reds')[8:9], brewer.pal(9, 'Blues')[5:9])
hosts5<-tick.snps5@pop
plot(pca.ticks5$scores[,1], pca.ticks5$scores[,2],
     col=alpha(color.host[hosts5], 0.5), pch=19, xlab='PCA 1', ylab='PCA 2')
dapc5<-dapc(tick.snps5, pop=hosts5, n.pca=5, n.da=2)
pdf(file='~/Dropbox/ddRADseq/dapc_by_individual_srt.pdf', height=10/2.54, width=10/2.54)
scatter(dapc5)
dev.off()

dapc5$posterior

##############################################################################################################
##############################################################################################################
#####                                           OLD STUFF                                                #####
##############################################################################################################
##############################################################################################################

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

## Related individuals within HB removed
tick.nr.snps<-read.snp('~/Desktop/D_variabilis_Pseudoref/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_HB_only_maxmissing0.75_MERGED_not_related_FOR_ADEGENET.finalSNPs.snp')
pca.ticks.nr<-glPca(tick.nr.snps)
scatter(pca.ticks.nr, posi="bottomleft")
loadingplot(pca.ticks.nr)
color.host=c('red', 'blue')
hosts<-tick.nr.snps@pop
plot(pca.ticks.nr$scores[,1], pca.ticks.nr$scores[,2],
     col=alpha(color.host[hosts], 0.5), pch=19, xlab='PCA 1', ylab='PCA 2')
dapc.nr<-dapc(tick.nr.snps, pop=hosts, n.pca=21, n.da=1)
#pdf(file='~/Dropbox/ddRADseq/dapc_by_individual.pdf', height=10/2.54, width=10/2.54)
scatter(dapc.nr)
#dev.off()

## Related removed, pop = host animal
tick.nr.host.snps<-read.snp('~/Desktop/D_variabilis_Pseudoref/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_HB_only_maxmissing0.75_MERGED_not_related_FOR_ADEGENET.finalSNPs_HOST.snp')
pca.ticks.nr.host<-glPca(tick.nr.host.snps)
scatter(pca.ticks.nr.host, posi="bottomleft")
loadingplot(pca.ticks.nr.host)
color.host=c(brewer.pal(9, 'Reds')[8:9], brewer.pal(9, 'Blues')[5:9])
hosts<-tick.nr.host.snps@pop
plot(pca.ticks.nr.host$scores[,1], pca.ticks.nr.host$scores[,2],
     col=alpha(color.host[hosts], 0.5), 
     pch=19, xlab='PCA 1', ylab='PCA 2')
dapc.nr.host<-dapc(tick.nr.host.snps, pop=hosts, n.pca=21, n.da=6)
#pdf(file='~/Dropbox/ddRADseq/dapc_by_individual.pdf', height=10/2.54, width=10/2.54)
scatter(dapc.nr.host, col=c(rep('blue', 2), rep('red', 5)))
col=c(rep('blue', 23), rep('red', 40))
#dev.off()
