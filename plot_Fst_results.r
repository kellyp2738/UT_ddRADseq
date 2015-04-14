library("VennDiagram")
library("scales")

## Analysis 3/2015

##########################
# Final, Both Sites

both<-read.table('~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Site/Final_Merged_Both_Sites_Bayescan_Input_fst.txt')
plot(both$qval, both$fst, pch=16, col='gray', cex=1.5, cex.axis=1.5, xlim=c(0,1), 
     ylim=c(0,0.15), bty='n', las=1, xlab='q-value', ylab="", cex.lab=1.5)
mean(both$fst)

##########################
# Final, SRT unrelated

SRT<-read.table('~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Host_SRT/Final_not_related_Merged_SRT_Bayescan_Input_fst.txt')
min(both$qval)
plot(SRT$qval, SRT$fst, pch=16, col='gray', cex=1.5, cex.axis=1.5, xlim=c(0,1), 
     ylim=c(0,0.15), bty='n', las=1, xlab='q-value', ylab="", cex.lab=1.5)
mean(SRT$fst)

##########################
# Final, HB unrelated

bs.not.related<-read.table('~/Dropbox/ddRADseq/D_variabilis_Pseudoref/Final_not_related_Merged_Bayescan_Input_fst.txt')
min(bs.not.related$qval)

png(file='~/Dropbox/ddRADseq/Final_Plots_April_2015/HB_fst_outlier.png', height=15, width=15, unit='cm', res=300)
par(mar=c(6,6,2,2))
plot(bs.not.related$qval, bs.not.related$fst, pch=16, col='gray', cex=1.5, cex.axis=1.5, xlim=c(0,1), 
     ylim=c(0,0.04),bty='n', las=1, xlab='q-value', ylab="", cex.lab=1.5)
mtext(text=expression('F'[ST]), side=2, line=4.3, cex=1.5)
dev.off()

mean.fst<-mean(bs.not.related$fst)
plot(density(bs.not.related$fst))
median(bs.not.related$fst)
min(bs.not.related$fst)
max(bs.not.related$fst)

##########################
# Final, All analyses

plot(both$qval, both$fst, pch=16, col='gray19', cex=1.5, cex.axis=1.5, xlim=c(0,1), 
     ylim=c(0,0.15), bty='n', las=1, xlab='q-value', ylab="", cex.lab=1.5)
points(SRT$qval, SRT$fst, pch=16, col='gray70', cex=1.5, cex.axis=1.5, xlim=c(0,1), 
       ylim=c(0,0.15), bty='n', las=1, xlab='q-value', ylab="", cex.lab=1.5)
points(bs.not.related$qval, bs.not.related$fst, pch=16, col='gray35', cex=1.5, cex.axis=1.5, xlim=c(0,1), 
       ylim=c(0,0.04),bty='n', las=1, xlab='q-value', ylab="", cex.lab=1.5)

###########################
# Overlap b/t samples

# unrelated ticks from HB... somehow this is not the right file!
HB.1<-read.table('~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Host_HB/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_HB_only_maxmissing0.75_MERGED_not_related_op_counts.frq.count', skip=1)
HB.2<-read.table('~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Host_HB/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_HB_only_maxmissing0.75_MERGED_not_related_racc_counts.frq.count', skip=1)
length(HB.1[,1])-length(HB.2[,1]) #because of how these were merged, they should automatically have the same sites... double checking that lengths match to be sure

# unrelated ticks from SRT
SRT.1<-read.table('~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Host_SRT/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_SRT_only_maxmissing0.75_MERGED_not_related_op_counts.frq.count', skip=1)
SRT.2<-read.table('~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Host_SRT/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_SRT_only_maxmissing0.75_MERGED_not_related_racc_counts.frq.count', skip=1)
length(SRT.1[,1])-length(SRT.2[,1])

# all ticks, all sites
both.1<-read.table('~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Site/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_BOTH_SITES_host_filtered_maxmissing0.75_MERGED_HB_counts.frq.count', skip=1)
both.2<-read.table('~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Site/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_BOTH_SITES_host_filtered_maxmissing0.75_MERGED_SRT_counts.frq.count', skip=1)
length(both.1[,1])-length(both.2[,1])

# all lengths match, so just use pop 1 to make Venn diagram

HB.id<-paste(HB.1[,1], HB.1[,2], sep="_")
SRT.id<-paste(SRT.1[,1], SRT.1[,2], sep="_")
both.id<-paste(both.1[,1], both.1[,2], sep="_")

# make a matrix of areas and overlaps
both.in.HB<-length(intersect(HB.id, both.id))
both.in.SRT<-length(intersect(SRT.id, both.id))
SRT.in.HB<-length(intersect(HB.id, SRT.id))
all.shared<-length(intersect(intersect(HB.id, SRT.id), both.id))

pdf(file='~/Dropbox/ddRADseq/Final_Plots_April_2015/Shared_SNPs_By_Site.pdf', height=10/2.54, width=10/2.54)
plot.new()
draw.pairwise.venn(area1=length(HB.id), area2=length(SRT.id), cross.area=2387,
                   category=c('HB', 'SRT'),
                   fill=c(alpha('gray19', 0.75), alpha('gray70', 0.75)), lwd=0 )
dev.off()