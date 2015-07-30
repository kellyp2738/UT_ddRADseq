#############################################################################################################
## ANALYSIS OF POPULATION STRUCTURE, TICK RADSEQ DATA
## Pierce et al. 2015 
#############################################################################################################

## INSTALL & LOAD PACKAGES
# dependency for pegas; not on CRAN for linux (on CRAN for OS X)
install.packages("~/Desktop//Rcpp/", repos=NULL,
                 type="source")
install.packages("~/Desktop//adegenet", repos=NULL,
                 type="source")
library("ape")
library("genetics")
library("pegas")
library("seqinr")
library("ggplot2")
library("RCo")
library("adegenet")
library("RColorBrewer")
library("scales")

# load the custom scatter function, adapted from scatter.dapc in adegenet
setwd() #set working directory to directory containing the analysis scripts
source("dapc_custom_scatter.r")
source('~/Dropbox/2014&Older/ModelFitting/R Scripts/PlotPrevCI.r')

## LOAD SUMMARY DATA FILE
extra.data<-read.csv('~/Desktop/UT_ddRADseq/ddRAD_FinalLibrary_SampleInfo_Full.csv')

## CHECK ADEGENET FUNCTIONALITY AND DATA FORMAT
packageDescription("adegenet", fields="Version")
# sample data in genlight format
file.show(system.file("files/exampleSnpDat.snp",package="adegenet"))


## 4 populations, infer
try.four<-find.clusters(site)
ul.try.four<-unlist(try.four$grp)
ul.df.try4<-cbind(names(ul.try.four), unlist(try.four$grp))
rownames(ul.df.try4)<-NULL
uldf<-data.frame(ul.df.try4)
names(uldf)<-c('ID', 'group')
ulm.df.try4<-merge(uldf, extra.data, by.x='ID', by.y='combo.label')#, by.x=ul.df.try4[,1], by.y='combo.label')
try4.results<-cbind(as.character(ulm.df.try4$ID), ulm.df.try4$group, 
                    paste(ulm.df.try4$host.species, ulm.df.try4$coll.site, sep='-'))
rt<-table(try4.results[,2], try4.results[,3])

## 3 populations, infer
try3<-find.clusters(site)
ul.try3<-unlist(try3$grp)
ul.df.try3<-cbind(names(ul.try3), unlist(try3$grp))
rownames(ul.df.try3)<-NULL
uldf3<-data.frame(ul.df.try3)
names(uldf3)<-c('ID', 'group')
ulm.df.try3<-merge(uldf3, extra.data, by.x='ID', by.y='combo.label')#, by.x=ul.df.try4[,1], by.y='combo.label')
try3.results<-cbind(as.character(ulm.df.try3$ID), ulm.df.try3$group, 
                    paste(ulm.df.try3$host.species, ulm.df.try3$coll.site, sep='-'))
rt3<-table(try3.results[,2], try3.results[,3])

## 2 populations, infer
try2<-find.clusters(site)
ul.try2<-unlist(try2$grp)
ul.df.try2<-cbind(names(ul.try2), unlist(try2$grp))
rownames(ul.df.try2)<-NULL
uldf2<-data.frame(ul.df.try2)
names(uldf2)<-c('ID', 'group')
ulm.df.try2<-merge(uldf2, extra.data, by.x='ID', by.y='combo.label')#, by.x=ul.df.try4[,1], by.y='combo.label')
try2.results<-cbind(as.character(ulm.df.try2$ID), ulm.df.try2$group, 
                    paste(ulm.df.try2$host.species, ulm.df.try2$coll.site, sep='-'))
rt2<-table(try2.results[,2], try2.results[,3])

barplot(t(rt))
barplot(t(rt3))
barplot(t(rt2))

par(mar=c(4,2,6,6))
table.value(rt)
table.value(rt3)
table.value(rt2)

setwd('~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Site_Host/')
site.host<-read.snp('Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_BOTH_SITES_host_filtered_only_maxmissing0.75_MERGED_FOR_ADEGENET.finalSNPs_site_and_host.snp')
site.host.names<-read.table('Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_BOTH_SITES_host_filtered_only_maxmissing0.75_MERGED_FOR_ADEGENET.finalSNPs_site_and_host.snp', 
                            skip=6, sep="\n")
tick.names<-as.character(site.host.names[seq(1, 198, 2),])
tick.names.formatted<-gsub('> ', '', tick.names, perl=TRUE)
sh.test<-dapc(site.host, n.pca=33, n.da=2)

## for EEID 2015 poster
pdf(file='~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Site_Host/2_way_dapc.pdf',
    height=10/2.54, width=10/2.54)
scatter(sh.test, legend=TRUE, posi.leg="topleft", label="", scree.da=FALSE,
        col=c('#08519c', '#a50f15', '#4292c6', '#ef3b2c'))
dev.off()

compoplot(sh.test)
find.clusters(site.host)

table.value(table(pop(site.host), sh.test$grp))

tally4<-0
assign<-c()
true<-c()
tf<-c()
for(i in 1:nInd(site.host)){
  x.rm4<-site.host[i] # remove individual i
  x.kp4<-site.host[-i] # keep all but individual i
  x.kp4.dapc<-dapc(x.kp4, n.pca=32, n.da=3)
  predict.x.rm4<-predict.dapc(x.kp4.dapc, newdata=x.rm4)
  if(as.character(predict.x.rm4$assign)==as.character(pop(x.rm4))){
    tally4=tally4+1
    tf<-c(tf, 1)
  }
  else{
    tf<-c(tf, 0)
  }
  assign<-c(assign, as.character(predict.x.rm4$assign))
  true<-c(true, as.character(pop(x.rm4)))
  #print(tally4)
}
tally4/nInd(site.host)
what.worked<-cbind(tick.names.formatted, assign, true, tf) #names vector is pulled from the snp data, so the order should be preserved
what.worked.meta<-merge(what.worked, extra.data, by.x='tick.names.formatted', by.y='combo.label')
what.worked.meta$host.sp.id<-paste(what.worked.meta$host.species, what.worked.meta$host.id, what.worked.meta$coll.site, sep='')
host.by.assign<-table(what.worked.meta$assign, what.worked.meta$host.sp.id)
barplot(host.by.assign, beside=TRUE)

barplot(what.worked)
hbo<-length(subset(what.worked, true=='HB-O')[,1])
hbr<-length(subset(what.worked, true=='HB-R')[,1])
srto<-length(subset(what.worked, true=='SRT-O')[,1])
srtr<-length(subset(what.worked, true=='SRT-R')[,1])
hbo.correct<-length(subset(what.worked, true=='HB-O' & tf=='1')[,1])
hbr.correct<-length(subset(what.worked, true=='HB-R' & tf=='1')[,1])
srto.correct<-length(subset(what.worked, true=='SRT-O' & tf=='1')[,1])
srtr.correct<-length(subset(what.worked, true=='SRT-R' & tf=='1')[,1])
hbo.pc<-hbo.correct/hbo #length(subset(what.worked, true=='HB-O' & tf=='1')[,1])/length(subset(what.worked, true=='HB-O'))
hbr.pc<-hbr.correct/hbr #length(subset(what.worked, true=='HB-R' & tf=='1')[,1])/length(subset(what.worked, true=='HB-R'))
srto.pc<-srto.correct/srto #length(subset(what.worked, true=='SRT-O' & tf=='1')[,1])/length(subset(what.worked, true=='SRT-O'))
srtr.pc<-srtr.correct/srtr #length(subset(what.worked, true=='SRT-R' & tf=='1')[,1])/length(subset(what.worked, true=='SRT-R'))

overall<-length(subset(what.worked, tf=='1')[,1])/length(what.worked[,1])

setwd('~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Site_Host/')
all.site.host<-read.snp('Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_BOTH_SITES_host_filtered_only_maxmissing0.75_MERGED_FOR_ADEGENET.finalSNPs_site_and_host.snp')
ash.test<-dapc(all.site.host, n.pca=33, n.da=2)
scatter(ash.test)
compoplot(sh.test)
find.clusters(site.host)




#################################################################################################################
##
##   POPULATION STRUCTURE BY SITE -- HARRISON BAYOU VS STARR RANCH TRAIL
##
#################################################################################################################

## POPULATION STRUCTURE EXISTS AT A FINE SCALE (~3km) -- COMPARISON OF HB TO SRT

## RELATED TICKS NOT REMOVED
# read in data
setwd('~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Site/')
site<-read.snp('Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_BOTH_SITES_host_filtered_only_maxmissing0.75_MERGED_FOR_ADEGENET.finalSNPs_site.snp')
# discriminant analysis of principal components to detect clusters
site.dapc<-dapc(site, pop=site.name, n.pca=33, n.da=1) # save as many PCs as possible
pdf(file='dapc_by_site_recolored.pdf', width=10/2.54, height=10/2.54)
par(mar=c(5,5,2,2))
site.colors2<-c('#fec44f', '#fee391')
scatter(site.dapc, col=c('#fec44f', '#fee391'), axes=FALSE)
axis(side=2, las=1)
legend("topright", legend=c('HB', 'SRT'), fill=alpha(site.colors2, 0.5), bty='n')
dev.off()

par(mfrow=c(1,1))
temp<-optim.a.score(site.dapc)

kept.id<-unlist(tapply(1:nInd(site), pop(site), function(e) sample(e, 5, replace=FALSE))) # 5 = number individuals from each population to exclude
x<-site[kept.id]
x.sup<-site[-kept.id]
x.sup.dapc<-dapc(x.sup, n.pca=29, n.da=1)
predict.x<-predict.dapc(x.sup.dapc, newdata=x)
mean(as.character(predict.x$assign)==as.character(pop(x)))

tally<-0
assign.site<-c()
true.site<-c()
tf.site<-c()
for(i in 1:nInd(site)){
  x.rm.site<-site[i] # remove individual i
  x.kp.site<-site[-i] # keep all but individual i
  x.kp.dapc.site<-dapc(x.kp.site, n.pca=32, n.da=1)
  predict.x.rm.site<-predict.dapc(x.kp.dapc.site, newdata=x.rm.site)
  if(as.character(predict.x.rm.site$assign)==as.character(pop(x.rm.site))){
    tally=tally+1
    tf.site<-c(tf.site, 1)
  }
  else{
    tf.site<-c(tf.site, 0)
  }
  print(as.character(predict.x.rm.site$assign))
  assign.site<-c(assign.site, as.character(predict.x.rm.site$assign))
  true.site<-c(true.site, as.character(pop(x.rm.site)))
  #print(tally)
}

what.worked.site<-cbind(tick.names.formatted, assign.site, true.site, tf.site) #names vector is pulled from the snp data, so the order should be preserved

hb<-length(subset(what.worked.site, true.site=='HB')[,1])
srt<-length(subset(what.worked.site, true.site=='SRT')[,1])
hb.correct<-length(subset(what.worked.site, true.site=='HB' & tf.site=='1')[,1])
srt.correct<-length(subset(what.worked.site, true.site=='SRT' & tf.site=='1')[,1])
hb.pc<-hb.correct/hb
srt.pc<-srt.correct/srt

accuracy.site<-(hb.correct+srt.correct)/length(what.worked.site[,1])

tally/nInd(site)

reg.clust<-find.clusters(site)

hier.clust<-find.clusters(site, clust=site@pop)
match<-cbind(site@)

## DAPC posterior cluster assignments
post.data.site<-data.frame(site.dapc$posterior) # extract posterior group assignments
post.df.site<-data.frame(cbind(row.names(post.data.site), post.data.site$HB, post.data.site$SRT)) # tweak the data structure to add sample names as a separate column (enables subsequent merging with data on collection site)
names(post.df.site)<-c('ID', 'HB', 'SRT')
post.df.site.c<-merge(post.df.site, extra.data, by.x='ID', by.y='combo.label') # add site data to posteriors data frame by merging on the sample names
post.df.site.trimmed<-cbind(post.df.site.c[,1:3], post.df.site.c$coll.site) # remove extraneous columns
predict.HB<-subset(post.df.site.trimmed, post.df.site.trimmed[,4]=='HB' & as.numeric(as.character(post.df.site.trimmed[,2]))>0.5) # correctly assigned to HB
predict.SRT<-subset(post.df.site.trimmed, post.df.site.trimmed[,4]=='SRT' & as.numeric(as.character(post.df.site.trimmed[,3]))>0.5) # correctly assigned to SRT
correct.HB<-length(predict.HB[,1]) # number correct HB
correct.SRT<-length(predict.SRT[,1]) # number correct SRT
accuracy.site<-(correct.HB+correct.SRT)/length(post.data.site[,1]) # accuracy of assignment to site groups

# plot posterior cluster assignments as a stacked barplot
post.df.site.trimmed.sort<-post.df.site.trimmed[order(post.df.site.trimmed$HB),]
colors.site<-c('gray19', 'gray70')
#pdf(file='dapc_by_site_posterior.pdf', width=25/2.54, height=10/2.54)
par(mar=c(6,4,2,8), xpd=TRUE)
barplot(t(post.df.site.trimmed.sort[2:3]), col=colors.site, axes=FALSE, ylab='Probability of Membership',
        xlim=c(0, length(post.df.site.trimmed.sort[,1])), space=0, col.axis='white') #for some reason #s still show on x-axis... make their color white
mtext(side=1, las=2, at=seq(0.5, length(post.df.site.trimmed.sort[,1])-0.5), cex=0.5, col=colors.site[post.df.site.trimmed.sort[,4]],
      text=gsub(pattern="_.*", replacement="", post.df.site.trimmed.sort[,1], perl=TRUE))
axis(side=2, las=1)
legend(x=101, y=0.55, legend=c('SRT', 'HB'), fill=c('gray19', 'gray70'), title='Cluster \nAssignment', bty="n")
#dev.off()

#################################################################################################################
##
##   POPULATION WITHIN HARRISON BAYOU
##
#################################################################################################################


## WITHIN HARRISON BAYOU, HOST SPECIES CORRESPONDS TO LOW-LEVEL POPULATION DIFFERENTIATION

# read in data
setwd('~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Host_HB/')
tick.snps2<-read.snp('Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_HB_only_maxmissing0.75_MERGED_not_related_FOR_ADEGENET.finalSNPs.snp')
# discriminant analysis of principal components to detect clusters
dapc2<-dapc(tick.snps2, pop=tick.snps2@pop, n.pca=21, n.da=1) 
#pdf(file='dapc_by_host.pdf', height=10/2.54, width=10/2.54)
par(mar=c(5,5,2,2))
scatter(dapc2)
axis(side=2, las=1)
legend("topleft", legend=c('Raccoon', 'Opossum'), fill=c('red', 'blue'), bty='n')
#dev.off()


tally2<-0
assign.hb<-c()
true.hb<-c()
tf.hb<-c()
for(i in 1:nInd(tick.snps2)){
  x.rm2<-tick.snps2[i] # remove individual i
  x.kp2<-tick.snps2[-i] # keep all but individual i
  x.kp2.dapc<-dapc(x.kp2, n.pca=21, n.da=1)
  predict.x.rm2<-predict.dapc(x.kp2.dapc, newdata=x.rm2)
  if(as.character(predict.x.rm2$assign)==as.character(pop(x.rm2))){
    tally2=tally2+1
    tf.hb<-c(tf.hb, 1)
  }
  else{
    tf.hb<-c(tf.hb, 0)
  }
  print(tally2)
  assign.hb<-c(assign.hb, as.character(predict.x.rm2$assign))
  true.hb<-c(true.hb, as.character(pop(x.rm2)))
}
tally2/nInd(tick.snps2)

what.worked.HB<-cbind(assign.hb, true.hb, tf.hb) #names vector is pulled from the snp data, so the order should be preserved

hb.racc<-length(subset(what.worked.HB, true.hb=='R')[,1])
hb.op<-length(subset(what.worked.HB, true.hb=='O')[,1])
hb.racc.correct<-length(subset(what.worked.HB, true.hb=='R' & tf.hb=='1')[,1])
hb.op.correct<-length(subset(what.worked.HB, true.hb=='O' & tf.hb=='1')[,1])
hb.racc.pc<-hb.racc.correct/hb.racc
hb.op.pc<-hb.op.correct/hb.op
accuracy.hb<-(hb.racc.correct+hb.op.correct)/length(what.worked.HB[,1])

# dapc posterior cluster assignments
post.data<-data.frame(dapc2$posterior)
post.df<-data.frame(cbind(row.names(post.data), as.numeric(post.data$O), as.numeric(post.data$R)))
names(post.df)<-c('tick', 'O', 'R')
post.df.host.type<-merge(post.df, extra.data, by.x='tick', by.y='combo.label')
post.df.ht.trimmed<-cbind(post.df.host.type[,1:3], post.df.host.type$host.species)
post.df.sort<-post.df.ht.trimmed[order(post.df.ht.trimmed[,2]),] # sort character data w/names in the same way
predict.HB.op<-subset(post.df.ht.trimmed, post.df.ht.trimmed[,4]=='op' & as.numeric(as.character(post.df.ht.trimmed[,2]))>0.5)
predict.HB.racc<-subset(post.df.ht.trimmed, post.df.ht.trimmed[,4]=='racc' & as.numeric(as.character(post.df.ht.trimmed[,3]))>0.5)
n.racc<-length(subset(post.df.ht.trimmed, post.df.ht.trimmed[,4]=='racc')[,1])
n.op<-length(subset(post.df.ht.trimmed, post.df.ht.trimmed[,4]=='op')[,1])
correct.HB.op<-length(predict.HB.op[,1]) #true positive opossum
correct.HB.racc<-length(predict.HB.racc[,1]) #true positive raccoon
pc.HB.op<-correct.HB.op/n.op
pc.HB.racc<-correct.HB.racc/n/racc
HB.accuracy<-(correct.HB.op+correct.HB.racc)/length(post.data[,1])

# DAPC plot
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
#tick.snps3<-read.snp('~/Desktop/D_variabilis_Pseudoref/Final_Analysis/Structure_by_Host_SRT/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_SRT_only_maxmissing0.75_MERGED_not_related_FOR_ADEGENET.finalSNPs.snp')
tick.snps3<-read.snp('~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Host_SRT/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_SRT_only_maxmissing0.75_MERGED_not_related_FOR_ADEGENET.finalSNPs.snp')
pca.srt<-glPca(tick.snps3)
scatter(pca.srt)

pdf(file='~/Dropbox/ddRADseq/Final_Plots_April_2015/dapc_by_host_SRT.pdf', height=10/2.54, width=10/2.54)
par(mar=c(5,5,2,2))
dapc3<-dapc(tick.snps3, pop=tick.snps3@pop, n.pca=5, n.da=1) 
scatter(dapc3)
axis(side=2, las=1)
legend("topright", legend=c('Raccoon', 'Opossum'), fill=c('red', 'blue'), bty='n')
dev.off()


tally.srt<-0
assign.srt<-c()
true.srt<-c()
tf.srt<-c()
for(i in 1:nInd(tick.snps3)){
  x.rm3<-tick.snps3[i] # remove individual i
  x.kp3<-tick.snps3[-i] # keep all but individual i
  x.kp3.dapc<-dapc(x.kp3, n.pca=4, n.da=1)
  predict.x.rm3<-predict.dapc(x.kp3.dapc, newdata=x.rm3)
  if(as.character(predict.x.rm3$assign)==as.character(pop(x.rm3))){
    tally.srt=tally.srt+1
    tf.srt<-c(tf.srt, 2)
  }
  else{
    tf.srt<-c(tf.srt, 0)
  }
  print(tally.srt)
  assign.ssrt<-c(assign.srt, as.character(predict.x.rm3$assign))
  true.srt<-c(true.srt, as.character(pop(x.rm3)))
}
tally.srt/nInd(tick.snps3)

# accidentally made true=2 and not true=1
what.worked.srt<-cbind(assign.srt, true.srt, tf.srt) #names vector is pulled from the snp data, so the order should be preserved
o.srt.correct<-length(subset(what.worked.srt, true.srt=='O' & tf.srt=='2')[,1])
o.srt.n<-length(subset(what.worked.srt, true.srt=='O')[,1])
r.srt.correct<-length(subset(what.worked.srt, true.srt=='R' & tf.srt=='2')[,1])
r.srt.n<-length(subset(what.worked.srt, true.srt=='R')[,1])
o.srt.pc<-o.srt.correct/o.srt.n
r.srt.pc<-r.srt.correct/r.srt.n
accuracy.srt<-(o.srt.correct+r.srt.correct)/length(what.worked.srt[,1])

# this seems wrong: accuracy.SRT<-(13+1)/15

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


########################
## Combined Plots     ##
########################

## DAPC

# both sites
pdf(file='~/Dropbox/ddRADseq/Final_Plots_April_2015/dapc_by_host_grayscale.pdf',
    height=15/2.54, width=15/2.54)
par(mfrow=c(1,1), mar=c(5,6,2,2))
scatter.custom(site.dapc, col=c('gray19', 'gray70'), cex.axis=1.5, cex.lab=1.5, ylim=c(0,0.8),
               hist.border.lty=c(3,1), hist.border.lwd=5, axes=FALSE, xlim=c(-4,6))
axis(side=1, cex.axis=1.5, at=seq(-4, 6, 2), las=1)
axis(side=2, cex.axis=1.5, las=1)
mtext(side=2, text='Density', line=4, cex=1.5)
legend.custom("topright", legend=c('HB', 'SRT'), fill=alpha(site.colors, 0.5), bty='n',
       border=site.colors, custom.border=c(3,1), custom.border.lwd=3, cex=1.5)
dev.off()
# HB
host.colors<-c('gray19', 'gray70')

pdf(file='~/Dropbox/ddRADseq/Final_Plots_April_2015/dapc_by_host_2_panel.pdf',
    height=15/2.54, width=30/2.54)
par(mfrow=c(1,2), mar=c(5,6,2,2))
scatter.custom(dapc2, col=host.colors, cex.axis=1.5, cex.lab=1.5, ylim=c(0,0.8),
               hist.border.lty=c(3,3), hist.border.lwd=5, axes=FALSE, xlim=c(-4,4))
axis(side=1, cex.axis=1.5, at=seq(-4, 4, 2))
axis(side=2, cex.axis=1.5, las=1)
mtext(side=2, text='Density', line=4, cex=1.5)
legend.custom("topleft", legend=c('Opossum', 'Raccoon'), fill=alpha(host.colors, 0.5), bty='n',
              border=host.colors, custom.border=c(3,3), custom.border.lwd=3, cex=1.5)
mtext(side=2, text="A", adj=6, las=1, padj=-12.5, cex=1.5)
# SRT
scatter.custom(dapc3, col=host.colors, cex.axis=1.5, cex.lab=1.5, ylim=c(0,0.8),
        hist.border.lty=c(1,1), hist.border.lwd=5, axes=FALSE, xlim=c(-6,4))
axis(side=1, cex.axis=1.5, at=seq(-6, 4, 2))
axis(side=2, cex.axis=1.5, las=1)
#mtext(side=2, text='Density', line=3, cex=1.5)
legend.custom("topleft", legend=c('Opossum', 'Raccoon'), fill=alpha(host.colors, 0.5), bty='n',
              border=host.colors, custom.border=c(1,1), custom.border.lwd=3, cex=1.5)
mtext(side=2, text="B", adj=6, las=1, padj=-12.5, cex=1.5)
dev.off()

## POSTERIOR CLUSTERING

pdf('~/Dropbox/ddRADseq/Final_Plots_April_2015/dapc_posterior_by_site_grayscale.pdf', 
    height=15/2.54, width=30/2.54)
par(mfrow=c(1,1), mar=c(5,5,2,8), xpd=TRUE)
# both
barplot(t(post.df.site.trimmed.sort[2:3]), col=site.colors, axes=FALSE, ylab='Probability of \nMembership',
        xlim=c(0, length(post.df.site.trimmed.sort[,1])), space=0, col.axis='white') #for some reason #s still show on x-axis... make their color white
mtext(side=1, las=2, at=seq(0.5, length(post.df.site.trimmed.sort[,1])-0.5), cex=0.5, col=site.colors[post.df.site.trimmed.sort[,4]],
      text=gsub(pattern="_.*", replacement="", post.df.site.trimmed.sort[,1], perl=TRUE))
axis(side=2, las=1)
legend(x=101, y=0.65, legend=c('SRT', 'HB'), fill=c('gray19', 'gray70'), title='Cluster \nAssignment', bty="n")
dev.off()

pdf('~/Dropbox/ddRADseq/Final_Plots_April_2015/dapc_posterior_by_host_grayscale.pdf', 
    height=30/2.54, width=30/2.54)
par(mfrow=c(2,1), mar=c(5,5,2,8), xpd=TRUE)
# HB
barplot(t(post.df.sort[2:3]), space=0, col.axis='white', ylab='Probability of \nMembership',
        las=2, cex.names=0.5, xlim=c(0,length(post.df.sort[,1])), axes=FALSE, col=site.colors)
mtext(side=1, las=2, at=seq(0.5, length(post.df.sort[,1])-0.5, by=1), text=gsub(pattern="_.*", replacement="", post.df.sort[,1], perl=TRUE),
      line=1, col=site.colors[post.df.sort[,4]], cex=0.75)
axis(side=2, las=1)
mtext(side=2, text="A", adj=5.3, las=1, padj=-12, cex=1.5)
legend(x=64, y=0.65, legend=c('Raccoon', 'Opossum'), fill=c(site.colors), title='Cluster \nAssignment', bty="n")

# SRT
label.cols3<-c(rep('gray70', 13), rep('gray19', 2))
barplot(t(post.df3.sort[2:3]), space=0, col.axis='white', ylab='Probability of \nMembership',
        las=2, cex.names=0.5, xlim=c(0,length(post.df3.sort[,1])), axes=FALSE, col=site.colors)
mtext(side=1, las=2, at=seq(0.5, length(post.df3.sort[,1])-0.5, by=1), text=gsub(pattern="_.*", replacement="", post.df3.sort[,1], perl=TRUE),
      line=1, col=label.cols3, cex=1)
axis(side=2, las=1)
mtext(side=2, text="B", adj=5.3, las=1, padj=-12, cex=1.5)
legend(x=15.2, y=0.65, legend=c('Raccoon', 'Opossum'), fill=site.colors, title='Cluster \nAssignment', bty="n")
dev.off()

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


### scatterplot of accuracy

## single model
sm.n<-c(hbo, hbr, srto, srtr)
sm.c<-c(hbo.correct, hbr.correct, srto.correct, srtr.correct)
sm.pc<-c(hbo.pc, hbr.pc, srto.pc, srtr.pc)
sm.data<-data.frame(cbind(sm.n, sm.c, sm.pc))
names(sm.data)<-c('n', 'pos', 'percent')
write.csv(sm.data, '~/Dropbox/ddRADseq/Final_Analysis/single_model_accuracy.csv')

## two part model

## by site
site.n<-c(hb, srt)
site.c<-c(hb.correct, srt.correct)
site.pc<-c(hb.pc, srt.pc)
site.data<-data.frame(cbind(site.n, site.c, site.pc))
names(site.data)<-c('n', 'pos', 'percent')
write.csv(site.data, '~/Dropbox/ddRADseq/Final_Analysis/two_part_by_site_accuracy.csv')

## HB only
hb.ro.ns<-c(hb.racc, hb.op)
hb.ro.c<-c(hb.racc.correct, hb.op.correct)
hb.ro.pc<-c(hb.racc.pc, hb.op.pc)
hb.ro.data<-data.frame(cbind(hb.ro.ns, hb.ro.c, hb.ro.pc))
names(hb.ro.data)<-c('n', 'pos', 'percent')
write.csv(hb.ro.data, '~/Dropbox/ddRADseq/Final_Analysis/two_part_hb_accuracy.csv')

## SRT only
srt.ro.ns<-c(r.srt.n, o.srt.n)
srt.ro.c<-c(r.srt.correct, o.srt.correct)
srt.ro.pc<-c(r.srt.pc, o.srt.pc)
srt.ro.data<-data.frame(cbind(srt.ro.ns, srt.ro.c, srt.ro.pc))
names(srt.ro.data)<-c('n', 'pos', 'percent')
write.csv(srt.ro.data, '~/Dropbox/ddRADseq/Final_Analysis/two_part_srt_accuracy.csv')

sm.data.ci<-KL.CI(sm.data$percent, sm.data$pos, sm.data$n, sm.data)
site.data.ci<-KL.CI(site.data$percent, site.data$pos, site.data$n, site.data)
hb.ro.data.ci<-KL.CI(hb.ro.data$percent, hb.ro.data$pos, hb.ro.data$n, hb.ro.data)
srt.ro.data.ci<-KL.CI(srt.ro.data$percent, srt.ro.data$pos, srt.ro.data$n, srt.ro.data)

plotCI(x=1:4, y=sm.data.ci$percent, ui=sm.data.ci$ci.u, li=sm.data.ci$ci.l)
plotCI(x=1:2, y=site.data.ci$percent, ui=site.data.ci$ci.u, li=site.data.ci$ci.l)
plotCI(x=1:2, y=hb.ro.data.ci$percent, ui=hb.ro.data.ci$ci.u, li=hb.ro.data.ci$ci.l)
plotCI(x=1:2, y=srt.ro.data.ci$percent, ui=srt.ro.data.ci$ci.u, li=srt.ro.data.ci$ci.l)

## Single Model

sm.percents<-cbind(c(sm.data$percent[3:4]), c(sm.data$percent[1:2]))
png(file='~/Dropbox/ddRADseq/Final_Analysis/n-1_accuracy_single_model.png',
    height=20, width=20, units='cm', res=300)
layout(mat=matrix(c(1,2,2), 1, 3 , byrow=TRUE), width=c(1,2))
par(mar=c(5, 6, 2, 2))
#layout.show()
barplot(overall*100, ylim=c(0,100), names=c('Overall'),
        las=1, ylab='', col='black', cex.axis=2, cex.lab=2, cex.names=2)
mtext(side=2, text='Percent Correct', line=4, cex=1.5)
barplot((sm.percents)*100, beside=T,
        ylim=c(0,100), col=c('blue', 'red', 'blue', 'red'),
        #col=c('#4292c6', '#ef3b2c', '#08519c', '#a50f15'), 
        names=c('Site 1', 'Site 2'), las=1, cex.axis=2, cex.lab=2, cex.names=2)
legend(x='topright', legend=c('Raccoon', 'Opossum'), fill=c('red', 'blue'), bty='n', cex=2)
dev.off()

## Two-part Model, by site
png(file='~/Dropbox/ddRADseq/Final_Analysis/n-1_accuracy_site.png',
    height=20, width=20, units='cm', res=300)
layout(mat=matrix(c(1,2), 1, 2 , byrow=TRUE), width=c(1.25,2))
barplot(accuracy.site*100, ylim=c(0,100), las=1, ylab='Percent Correct',
        cex.lab=1.5, cex.axis=1.5, cex.names=1.5, names=c('Overall'), col='black')
barplot(site.data$percent*100, col=site.colors2[1], ylim=c(0,100), ylab='', las=1,
        names=c('Site 1', 'Site 2'), cex.lab=1.5, cex.names=1.5, cex.axis=1.5)
dev.off()
combined.site<-cbind(hb.ro.data$percent, srt.ro.data$percent)

right.order<-rbind((combined.site[2,]), (combined.site[1,]))
png(file='~/Dropbox/ddRADseq/Final_Analysis/n-1_accuracy_HB.png',
    height=20, width=20, units='cm', res=300)
layout(mat=matrix(c(1,2), 1, 2 , byrow=TRUE), width=c(1.25,2))
barplot(accuracy.hb*100, ylim=c(0,100), ylab='Percent Correct', las=1, col='black',
        cex.names=1.5, cex.axis=1.5, cex.lab=1.5, names=c('Overall Site 1'))
barplot((right.order[,1])*100, beside=T, col=c('blue', 'red'),
        cex.lab=1.5, cex.axis=1.5, cex.names=1.5,
        names=c('Site 1'), space=0, ylab='', las=1, ylim=c(0,100))
dev.off()

png(file="~/Dropbox/ddRADseq/Final_Analysis/n-1_accuracy_SRT.png",
    height=20, width=20, units='cm', res=300)
layout(mat=matrix(c(1,2), 1, 2 , byrow=TRUE), width=c(1.25,2))
barplot(accuracy.srt*100, names=c('Overall'), ylim=c(0,100),
        las=1, ylab=c('Percent Correct'), cex.lab=1.5,
        cex.names=1.5, cex.axis=1.5, col='black')
barplot((right.order[,2])*100, beside=T, col=c('blue', 'red'),
        cex.lab=1.5, cex.names=1.5, cex.axis=1.5,
        names=c('Site 2'), space=0, ylab='', las=1, ylim=c(0,100))
#legend(x='topright', legend=c('Raccoon', 'Opossum'), fill=c('red', 'blue'), bty='n')
dev.off()

all.accuracy<-c(overall, accuracy.site, accuracy.hb, accuracy.srt)
sample.sizes<-c(length(what.worked[,1]), length(what.worked.site[,1]), length(what.worked.HB[,1]), length(what.worked.srt[,1]))
n.pos<-c(length(subset(what.worked, tf=='1')[,1]), (hb.correct+srt.correct), (hb.racc.correct+hb.op.correct), (o.srt.correct+r.srt.correct))
all.for.CI<-cbind(all.accuracy, sample.sizes, n.pos)
all.CIs<-KL.CI(all.for.CI[,3], all.for.CI[,2], all.for.CI)


library(gplots)

png(file='~/Dropbox/ddRADseq/Final_Analysis/n-1_accuracy.png',
    height=10, width=20, units='cm', res=300)
par(mfrow=c(1,1))
barplot2(all.accuracy*100, ylim=c(0,100), plot.ci=TRUE, ci.u=all.CIs[,4]*100, ci.l=all.CIs[,5]*100,
        names=c('Single Model', 'Site', 'Host | Site 1', 'Host | Site 2'),
        las=1, cex.lab=1.5, cex.names=1.5, cex.axis=1.5,
        ylab=c('Percent Accuracy'), col='slategray2')
dev.off()

# average accuracy by subpopulation assigment (so that it doesn't look like high accuracy when the only thing that's ever right is Raccoon and/or site 1)

mean.acc<-c(mean(sm.data$percent), mean(site.data$percent),
            mean(hb.ro.data$percent), mean(srt.ro.data$percent))
all.accuracy.mean<-cbind(all.for.CI, mean.acc)
all.CIs.mean<-KL.CI(all.accuracy.mean[,4]*all.accuracy.mean[,2], all.accuracy.mean[,2],
                    all.accuracy.mean)

png(file='~/Dropbox/ddRADseq/Final_Analysis/n-1_accuracy_average_by_category.png',
    height=10, width=22, units='cm', res=300)
par(mfrow=c(1,1), mar=c(4,6,2,1))
barplot2(mean.acc*100, ylim=c(0,100), plot.ci=TRUE, ci.u=all.CIs.mean[,5]*100, ci.l=all.CIs.mean[,6]*100,
         names=c('Single Model', 'Site', 'Host | Site 1', 'Host | Site 2'),
         las=1, cex.lab=1.5, cex.names=1.5, cex.axis=1.5,
         ylab='', col='slategray2')
mtext(side=2, line=4, text=c('Percent Accuracy'), cex=1.5)
dev.off()