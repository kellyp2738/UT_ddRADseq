## -----------------------------------------------------------------------------------
## -- PARSE HAPLOTYPES FILE FOR BAYESCAN
## -----------------------------------------------------------------------------------

library(stringr)
library(plyr)
library(magicaxis)
source('~/BayeScan2.1/R functions/plot_R.r')

## NOTE: number of columns may need adjustment if re-running

# - all the SNPs
#pop.1<-read.table('/Users/kellypierce/Dropbox/ddRADseq/var_maf0.1_minmeanDP20_minGQ25_opossums.frq.count', skip=1)
#names(pop.1)<-c('CHROM', 'POS', 'N_ALLELES', 'N_CHR', 'ALLELE.COUNT.1', 'ALLELE.COUNT.2')
#pop.2<-read.table('/Users/kellypierce/Dropbox/ddRADseq/var_maf0.1_minmeanDP20_minGQ25_raccoons.frq.count', skip=1)
#names(pop.2)<-c('CHROM', 'POS', 'N_ALLELES', 'N_CHR', 'ALLELE.COUNT.1', 'ALLELE.COUNT.2')

# - SNPs found in at least 75% of individuals
#pop.1<-read.table('/Users/kellypierce/Dropbox/ddRADseq/var_maf0.1_minmeanDP20_minGQ25_maxmissing0.75_opossums.frq.count', skip=1)
#names(pop.1)<-c('CHROM', 'POS', 'N_ALLELES', 'N_CHR', 'ALLELE.COUNT.1', 'ALLELE.COUNT.2')
#pop.2<-read.table('/Users/kellypierce/Dropbox/ddRADseq/var_maf0.1_minmeanDP20_minGQ25_maxmissing0.75_raccoons.frq.count', skip=1)
#names(pop.2)<-c('CHROM', 'POS', 'N_ALLELES', 'N_CHR', 'ALLELE.COUNT.1', 'ALLELE.COUNT.2')

# - D. variabilis pseudoreference mapped SNPs found in at least 75% of individuals
#pop.1<-read.table('~/Dropbox/ddRADseq/pseudoref_opossum_counts_maf0.1_minmeanDP20_minGQ25_maxmissing0.75.frq.count', skip=1)
#names(pop.1)<-c('CHROM', 'POS', 'N_ALLELES', 'N_CHR', 'ALLELE.COUNT.1', 'ALLELE.COUNT.2')
#pop.2<-read.table('~/Dropbox/ddRADseq/pseudoref_raccoon_counts_maf0.1_minmeanDP20_minGQ25_maxmissing0.75.frq.count', skip=1)
#names(pop.2)<-c('CHROM', 'POS', 'N_ALLELES', 'N_CHR', 'ALLELE.COUNT.1', 'ALLELE.COUNT.2')

# - D. variabilis pseudoreference mapped SNPs found in at least 75% of individuals 
# - (filtered at whole population level first)
#pop.1<-read.table('~/Dropbox/ddRADseq/pseudoref_opossum_counts_pop_filtered_maxmissing0.75.frq.count', skip=1)
#names(pop.1)<-c('CHROM', 'POS', 'N_ALLELES', 'N_CHR', 'ALLELE.COUNT.1', 'ALLELE.COUNT.2')
#pop.2<-read.table('~/Dropbox/ddRADseq/pseudoref_raccoon_counts_pop_filtered_maxmissing0.75.frq.count', skip=1)
#names(pop.2)<-c('CHROM', 'POS', 'N_ALLELES', 'N_CHR', 'ALLELE.COUNT.1', 'ALLELE.COUNT.2')

## Update 3/2015 (yes, I'm still working on this months after graduation. don't judge me)

pop.1<-read.table('~/Desktop/UT_ddRADseq/Updated_pseudoref_counts_maxmissing_0.75_minmeanDP20_minGQ25_opHB_maf0.1_COUNTS.frq.count', skip=1)
pop.2<-read.table('~/Desktop/UT_ddRADseq/Updated_pseudoref_counts_maxmissing_0.75_minmeanDP20_minGQ25_raccHB_maf0.1_COUNTS.frq.count', skip=1)
names(pop.1)<-c('CHROM', 'POS', 'N_ALLELES', 'N_CHR', 'ALLELE.COUNT.1', 'ALLELE.COUNT.2')
names(pop.2)<-c('CHROM', 'POS', 'N_ALLELES', 'N_CHR', 'ALLELE.COUNT.1', 'ALLELE.COUNT.2')

# parse alleles, population 1
pop.1.allele.1<-unlist(strsplit(as.character(pop.1$ALLELE.COUNT.1), split=':'))
pop.1.allele.1.ID<-pop.1.allele.1[seq(1,length(pop.1.allele.1),2)]
pop.1.allele.1.count<-pop.1.allele.1[seq(2,length(pop.1.allele.1),2)]
pop.1.allele.1.parsed<-cbind(pop.1.allele.1.ID, pop.1.allele.1.count)

pop.1.allele.2<-unlist(strsplit(as.character(pop.1$ALLELE.COUNT.2), split=':'))
pop.1.allele.2.ID<-pop.1.allele.2[seq(1,length(pop.1.allele.2),2)]
pop.1.allele.2.count<-pop.1.allele.2[seq(2,length(pop.1.allele.2),2)]
pop.1.allele.2.parsed<-cbind(pop.1.allele.2.ID, pop.1.allele.2.count)

# parse alleles, population 2
pop.2.allele.1<-unlist(strsplit(as.character(pop.2$ALLELE.COUNT.1), split=':'))
pop.2.allele.1.ID<-pop.2.allele.1[seq(1,length(pop.2.allele.1),2)]
pop.2.allele.1.count<-pop.2.allele.1[seq(2,length(pop.2.allele.1),2)]
pop.2.allele.1.parsed<-cbind(pop.2.allele.1.ID, pop.2.allele.1.count)

pop.2.allele.2<-unlist(strsplit(as.character(pop.2$ALLELE.COUNT.2), split=':'))
pop.2.allele.2.ID<-pop.2.allele.2[seq(1,length(pop.2.allele.2),2)]
pop.2.allele.2.count<-pop.2.allele.2[seq(2,length(pop.2.allele.2),2)]
pop.2.allele.2.parsed<-cbind(pop.2.allele.2.ID, pop.2.allele.2.count)

# combine into a single dataframe containing all info for making a bayescan file
pop.1.parsed<-cbind(pop.1[,1:4], pop.1.allele.1.parsed, pop.1.allele.2.parsed)
pop.2.parsed<-cbind(pop.2[,1:4], pop.2.allele.1.parsed, pop.2.allele.2.parsed)
full.pop.parsed<-merge(pop.1.parsed, pop.2.parsed, by=c('CHROM', 'POS'), all=TRUE)
max(full.pop.parsed$N_ALLELES.x, na.rm=TRUE)
max(full.pop.parsed$N_ALLELES.y, na.rm=TRUE)

# add a sequence ID column
full.pop.parsed<-cbind(seq(1:length(full.pop.parsed[,1])), full.pop.parsed)
names(full.pop.parsed)<-c('seq.ID', 'pos', 'chrom', 'n.alleles.pop1',
                          'n.chrom.pop1', 'seq.pop1.p', 'count.pop1.p',
                          'seq.pop1.q', 'count.pop1.q', 'n.alleles.pop2',
                          'n.chrom.pop2', 'seq.pop2.p', 'count.pop2.p', 
                          'seq.pop2.q', 'count.pop2.q')
# generate empty data frames for storing the bayescan data
bayescan1<-data.frame(matrix(nrow=length(full.pop.parsed[,1]), ncol=10))
bayescan2<-data.frame(matrix(nrow=length(full.pop.parsed[,1]), ncol=10))

# for each row in the parsed data,
for(row in full.pop.parsed$seq.ID){
  
  # make a table of alleles in population 1
  full.pop.parsed.p1.alleles<-c(as.character(full.pop.parsed$seq.pop1.p[row]), 
                                as.character(full.pop.parsed$seq.pop1.q[row]))
  full.pop.parsed.p1.counts<-c(as.character(full.pop.parsed$count.pop1.p[row]), 
                               as.character(full.pop.parsed$count.pop1.q[row]))
  full.pop.parsed.p1.table<-data.frame(cbind(full.pop.parsed.p1.alleles, full.pop.parsed.p1.counts))
  names(full.pop.parsed.p1.table)<-c('sequence', 'count')
  #print(full.pop.parsed.p1.alleles)
  #print(full.pop.parsed.p1.counts)
  
  # make a table of alleles in population 2
  full.pop.parsed.p2.alleles<-c(as.character(full.pop.parsed$seq.pop2.p[row]), 
                                as.character(full.pop.parsed$seq.pop2.q[row]))
  full.pop.parsed.p2.counts<-c(as.character(full.pop.parsed$count.pop2.p[row]), 
                               as.character(full.pop.parsed$count.pop2.q[row]))
  full.pop.parsed.p2.table<-data.frame(cbind(full.pop.parsed.p2.alleles, full.pop.parsed.p2.counts))
  names(full.pop.parsed.p2.table)<-c('sequence', 'count')
  
  # merge the two tables to summarize the population differences
  full.pop.parsed.compare<-merge(full.pop.parsed.p1.table, full.pop.parsed.p2.table, 
                                 by='sequence', all=TRUE)
  #print(full.pop.parsed.compare)
  
  # some populations will not have any data for a locus. remove those completely empty rows
  not.missing<-full.pop.parsed.compare[rowSums(is.na(full.pop.parsed.compare)) != ncol(full.pop.parsed.compare),]
  
  # build the lines for the bayescan file
  line.for.bayescan1<-c(full.pop.parsed$seq.ID[row], sum(as.numeric(as.character(not.missing$count.x))), full.pop.parsed$n.alleles.pop1[row], as.numeric(as.character(not.missing$count.x)))
  line.for.bayescan2<-c(full.pop.parsed$seq.ID[row], sum(as.numeric(as.character(not.missing$count.y))), full.pop.parsed$n.alleles.pop2[row], as.numeric(as.character(not.missing$count.y)))
  
  # turn the NAs into 0s
  line.for.bayescan1[is.na(line.for.bayescan1)]<-0
  line.for.bayescan2[is.na(line.for.bayescan2)]<-0
  
  # add them to the file
  bayescan1[row, 1:length(line.for.bayescan1)]<-line.for.bayescan1
  bayescan2[row, 1:length(line.for.bayescan2)]<-line.for.bayescan2
  
  #print(not.missing)
  #print(full.pop.parsed[row,])
  #print(line.for.bayescan1)
  #print(line.for.bayescan2)
}

missing.1<-subset(bayescan1, X2==0)
missing.2<-subset(bayescan2, X2==0)

tu<-union(missing.1[,1], missing.2[,1]) # all the IDs that aren't shared
bs1.keep<-subset(bayescan1, bayescan1[,1] %in% tu == FALSE)
bs2.keep<-subset(bayescan2, bayescan2[,1] %in% tu == FALSE)
bs1.keep$X1<-seq(1,length(bs1.keep[,1]),1)
bs2.keep$X1<-seq(1,length(bs2.keep[,1]),1)

#write.table(bs1.keep, file='~/Dropbox/ddRADseq/var_maf0.1_minmeanDP20_minGQ25_maxmissing0.75_opossums_NotMissing_Bayescan', quote=FALSE, row.names=FALSE, col.names=FALSE)
#write.table(bs2.keep, file='~/Dropbox/ddRADseq/var_maf0.1_minmeanDP20_minGQ25_maxmissing0.75_raccoons_NotMissing_Bayescan', quote=FALSE, row.names=FALSE, col.names=FALSE)

#write.table(bs1.keep, file='~/Dropbox/ddRADseq/pseudoref_maf0.1_minmeanDP20_minGQ25_maxmissing0.75_opossums_NotMissing_Bayescan', quote=FALSE, row.names=FALSE, col.names=FALSE)
#write.table(bs2.keep, file='~/Dropbox/ddRADseq/pseudoref_maf0.1_minmeanDP20_minGQ25_maxmissing0.75_raccoons_NotMissing_Bayescan', quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(bs1.keep, file='~/Dropbox/ddRADseq/pseudoref_pop_filtered_maxmissing0.75_opossums_NotMissing_Bayescan', quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(bs2.keep, file='~/Dropbox/ddRADseq/pseudoref_pop_filtered_maxmissing0.75_raccoons_NotMissing_Bayescan', quote=FALSE, row.names=FALSE, col.names=FALSE)


## most updated Bayescan output 21/Aug/2014

# Outlier loci, genotypes, and position in the Stacks catalog
pd<-read.table('/Users/kellypierce/Dropbox/ddRADseq/D_variabilis_Pseudoref/Bayescan_Input_pseudoref_pop_filtered_maxmissing0.75_opossums_NotMis_fst.txt')
plot(pd$qval, pd$fst)
plot_bayescan(pd) # the outlier is #1581 in the bayescan input file
mean(pd$fst)

# pick a random non-significant locus
med.sig<-subset(pd, qval<0.75)
rand.locus<-sample(x=seq(1, length(med.sig[,1]),1), size=1) # in this case, 229
med.sig[rand.locus,]
med.locus<-full.pop.parsed[229,]

# make non-outlier locus allele frequency barplots
non.outlier<-matrix(c(52/62,10/62,118/134,16/134), nrow=2, ncol=2) #Fst=0.004

png(file='~/Dropbox/ddRADseq/D_variabilis_Pseudoref/NonOutlierLocus_PopFiltered.png', height=10, width=15, unit='cm', res=300)
#par(bg="black", col.axis="white",col.lab="white", col.main="white", col.sub="white", fg="white")
par(mfrow=c(1,1), mar=c(5,6,2,4))
barplot(as.matrix(non.outlier), col=c('#636363', '#bdbdbd'),
        xlim=c(0.2,2.7), width=0.8, border=NA,
        cex.axis=1.5, cex.names=1.5, cex.lab=1.5, las=1)
mtext(text=c('Allele Frequency'), side=2, line=4, cex=1.5)
mtext(text=c('Opossum\nTicks', 'Raccoon\nTicks'), side=1, line=3, 
      at=c(0.55,1.54), cex=1.5)
legend(x='topright', legend=c('G', 'T'), border=NA,
       fill=c('#636363', '#bdbdbd'), bty='n', cex=1.5)
dev.off()

# outlier at q<0.05
bs1.keep[1581,] #1951 in the parsed haplotype dataframe
full.pop.parsed[1951,] #catalog locus 5394
cl.5394<-subset(full.pop.parsed, pos=="5394_pseudoreference_pe_concatenated_without_rev_complement")
cl.5394 #only 1 SNP for this fragment
pd[1581,]

# back to the original data
setwd('~/Dropbox/ddRADseq/paired_denovo_outputs_m10_Aug2014')

# Import the data and check the structure
tags<-read.table('batch_1.catalog.tags.tsv', header=FALSE)
outlier<-subset(tags, V3==5394)
out.seq<-unlist(outlier)
seq<-c('C','A','T','G','C','A','G','G','T','C','A','C','G','G','A','T','C','G','T','G','A','T','G','C','G','A','G','G','T','T','C','A','T','A','C','T','G','A','A','G','G',
       'C','G','A','C','A','C','C','G','A','T','A','C','C','C','A','G','C','A','G','T','A','G','T','G','C','T','A','G','C','A','C','G','G','G','T','C','C','G','G','C','G',
       'A','A','G','C','A','A','G','C','T','G','C','C','A','T','A','A','T','T','C','T','G','A','C','C','G','C','A','C','G','C','A','T','T','G','G','C','G','A','C','C','T',
       'G','G','T','G','C','C','A','G','C','G','T','G','G','T','G','C','T','T','G','A','C','C','A','T','T','A','C','A','T','T','A','G','G','C','A','T','T','G','C','A','T',
       'T','G','G','C','C','G','G','G','C','T','C','G','T','C','T','G','G','T','T','C','A','C','A','T','C','G','T','C','G','G','C','G','C')
length(seq)
seq[1:100]
read1<-c()
for(i in 1:100){
  #print(seq[i])
  read1<-paste(read1, seq[i], sep='')
}
marginal.outlier<-subset(tags, V3==6275)

seq<-outlier$V9
# outlier at q<0.1
low.q<-subset(pd, qval<0.1) #the only other one is #1657 in bayescan file
bs1.keep[1657,] #2059
full.pop.parsed[2059,] #catalog locus 6275
cl.6275<-subset(full.pop.parsed, pos=="6275_pseudoreference_pe_concatenated_without_rev_complement")
cl.6275 #only 1 SNP for this fragment

# plot overall fst results
png(file='~/Dropbox/ddRADseq/D_variabilis_Pseudoref/BayescanFst_PopFiltered_Output_BlackBG.png', height=10, width=15, unit='cm', res=300)
par(bg="black", col.axis="white",col.lab="white", col.main="white", col.sub="white", fg="white")
par(mfrow=c(1,1), mar=c(5,7,4,2))
plot(pd$qval, pd$fst, pch=16, col='gray', cex=1.5, cex.axis=1.5, xlim=c(0,1), ylim=c(0,0.06),
     bty='n', las=1, xlab='q-value', ylab="", cex.lab=1.5)
#abline(v=0.05)
points(x=subset(pd, qval==min(pd$qval))$qval,
       y=subset(pd, qval==min(pd$qval))$fst,
       col=c('#74a9cf'), pch=16, cex=1.5)
points(x=pd$qval[1657],y=pd$fst[1657],
       col=c('#d95f0e'), pch=16, cex=1.5)
mtext(text='Fst', side=2, line=4.3, cex=1.5)
legend(x='topright', legend=c('Significant, q<0.05', 'Significant, q<0.1'), pch=16,
       border=NA, col=c('#74a9cf', '#d95f0e'), bty='n', cex=1)
dev.off()

# plot overal alpha results
png(file='~/Dropbox/ddRADseq/D_variabilis_Pseudoref/BayescanAlpha_PopFiltered_Output.png', height=10, width=15, unit='cm', res=300)
#par(bg="black", col.axis="white",col.lab="white", col.main="white", col.sub="white", fg="white")
par(mfrow=c(1,1), mar=c(5,7,4,2))
plot(pd$qval, pd$alpha, pch=16, col='gray', cex=1.5, cex.axis=1.5, xlim=c(0,1), ylim=c(0,2.5),
     bty='n', las=1, xlab='q-value', ylab="", cex.lab=1.5)
#abline(v=0.05)
points(x=subset(pd, qval==min(pd$qval))$qval,
       y=subset(pd, qval==min(pd$qval))$alpha,
       col=c('#74a9cf'), pch=16, cex=1.5)
points(x=pd$qval[1657],y=pd$alpha[1657],
       col=c('#d95f0e'), pch=16, cex=1.5)
mtext(text='alpha', side=2, line=4.3, cex=1.5)
legend(x='topright', legend=c('Significant, q<0.05', 'Significant, q<0.1'), pch=16,
       border=NA, col=c('#74a9cf', '#d95f0e'), bty='n', cex=1)
dev.off()

# make outlier locus allele frequency barplots
outlier<-matrix(c(1,0,85/114,29/114), nrow=2, ncol=2) #Fst=0.05
marginal.outlier<-matrix(c(16/58, 42/58, 5/132, 127/132), nrow=2, ncol=2)

png(file='~/Dropbox/ddRADseq/D_variabilis_Pseudoref/OutlierLocus_PopFiltered_BlackBG.png', height=10, width=15, unit='cm', res=300)
par(bg="black", col.axis="white",col.lab="white", col.main="white", col.sub="white", fg="white")
par(mfrow=c(1,1), mar=c(5,6,2,4))
barplot(as.matrix(outlier), col=c('#74a9cf', '#bdc9e1'),
        xlim=c(0.2,2.7), width=0.8, border=NA,
        cex.axis=1.5, cex.names=1.5, cex.lab=1.5, las=1)
mtext(text=c('Allele Frequency'), side=2, line=4, cex=1.5)
mtext(text=c('Opossum\nTicks', 'Raccoon\nTicks'), side=1, line=3, 
      at=c(0.55,1.54), cex=1.5)
legend(x='topright', legend=c('A', 'C'), border=NA,
       fill=c('#74a9cf', '#bdc9e1'), bty='n', cex=1.5)
dev.off()

png(file='~/Dropbox/ddRADseq/D_variabilis_Pseudoref/MarginalOutlierLocus_PopFiltered_BlackBG.png', height=10, width=15, unit='cm', res=300)
par(bg="black", col.axis="white",col.lab="white", col.main="white", col.sub="white", fg="white")
par(mfrow=c(1,1), mar=c(5,6,2,4))
barplot(as.matrix(marginal.outlier), col=c('#d95f0e', '#fec44f'),
        xlim=c(0.2,2.7), width=0.8, border=NA,
        cex.axis=1.5, cex.names=1.5, cex.lab=1.5, las=1)
mtext(text=c('Allele Frequency'), side=2, line=4, cex=1.5)
mtext(text=c('Opossum\nTicks', 'Raccoon\nTicks'), side=1, line=3, 
      at=c(0.55,1.54), cex=1.5)
legend(x='topright', legend=c('G', 'T'), border=NA,
       fill=c('#d95f0e', '#fec44f'), bty='n', cex=1.5)
dev.off()


## Output of Bayescan

# no significant outliers for Fst
pseudoref.data<-read.table('/Users/kellypierce/Dropbox/ddRADseq/D_variabilis_Pseudoref/Bayescan_Input_pseudoref_maf0.1_minmeanDP20_minGQ25_maxmissing0.75_opossums_NotMis_fst.txt')

png(file='~/Dropbox/ddRADseq/D_variabilis_Pseudoref/BayescanFst_Output_BlackBG.png', height=10, width=15, unit='cm', res=300)
par(mfrow=c(1,1), mar=c(5,7,4,2))
par(bg="black", col.axis="white",col.lab="white", col.main="white", col.sub="white", fg="white")
plot(pseudoref.data$qval, pseudoref.data$fst, pch=16, col='gray', cex=1.2, cex.axis=1.5,
     xlim=c(0.65, 0.95), ylim=c(0,0.008), bty="n", las=1, xlab='q-value', ylab="", cex.lab=1.5)
#points(x=subset(pseudoref.data, qval==min(pseudoref.data$qval))$qval,
#       y=subset(pseudoref.data, qval==min(pseudoref.data$qval))$fst,
#       col=c('#74a9cf'), pch=16, cex=1.5)
mtext(text='Fst', side=2, line=5, cex=1.5)
dev.off()

#plot_bayescan(pseudoref.data)

# the most differentiated locus (which is nowhere near significant)

best<-subset(pseudoref.data, qval==min(pseudoref.data$qval)) #237
best.op<-bs1.keep[237,4:5] #locus ID 422
best.racc<-bs2.keep[237,4:5] #locus ID 422

full.pop.parsed[422,] # get the allele identities (catalog locus 1392)
# in opossums, 20 A and 36 G
# in raccoons, 90 A and 42 G

op.total<-sum(best.op)
racc.total<-sum(best.racc)

best.freqs<-matrix(data=c(best.op[1]/op.total, best.op[2]/op.total, best.racc[1]/racc.total, best.racc[2]/racc.total),
                   nrow=2)

png(file='~/Dropbox/ddRADseq/D_variabilis_Pseudoref/SampleLocus_BlackBG.png', height=10, width=15, unit='cm', res=300)
par(bg="black", col.axis="white",col.lab="white", col.main="white", col.sub="white", fg="white")
par(mfrow=c(1,1), mar=c(5,6,2,4))
barplot(as.matrix(best.freqs), col=c('#74a9cf', '#bdc9e1'),
        xlim=c(0.2,2.7), width=0.8, border=NA,
        cex.axis=1.5, cex.names=1.5, cex.lab=1.5, las=1)
mtext(text=c('Allele Frequency'), side=2, line=4, cex=1.5)
mtext(text=c('Opossum\nTicks', 'Raccoon\nTicks'), side=1, line=3, 
      at=c(0.55,1.54), cex=1.5)
legend(x='topright', legend=c('A', 'G'), border=NA,
       fill=c('#74a9cf', '#bdc9e1'), bty='n', cex=1.5)
dev.off()

#col=c('#253494','#2c7fb8') alt colors
#col=c('#41b6c4', '#a1dab4')



## ------------------------------------------------------------------------------------------
## Old outputs from Stacks filtering

#better.snp<-plot_bayescan('/Users/kellypierce/Dropbox/ddRADseq/var_maf0.1_minmeanDP20_minGQ25_maxmissing0.75_NotMissing_Bayescan_I_fst.txt', FDR=0.05)
better.snp.data<-read.table('/Users/kellypierce/Dropbox/ddRADseq/var_maf0.1_minmeanDP20_minGQ25_maxmissing0.75_NotMissing_Bayescan_I_fst.txt')
plot(better.snp.data$qval, better.snp.data$fst)

png(file='Fst_Alpha_mpileup_denovo_reference.png', width=15, height=15, unit='cm', res=300)
par(mfrow=c(2,1),mar=c(4,6,1,2))
plot((better.snp.data$qval), better.snp.data$fst, las=1, pch=16, cex.lab=1.5, cex.axis=1.2, 
     xlab='', ylab="", col='gray')#, axes=FALSE)
#magaxis(side=1, unlog=TRUE, box=TRUE, cex.axis=1.2)
#magaxis(side=1, unlog=FALSE, box=TRUE, cex.axis=1.2)
#magaxis(side=2, las=1, cex.axis=1.2)
mtext('Fst', side=2, line=4, cex=1.2)
mtext('q-value', side=1, line=2.5, cex=1.2)
# positive value of alpha suggests diversifying selection
plot((better.snp.data$qval), better.snp.data$alpha, las=1, ylab='', xlab='', pch=16,
     cex.lab=1.5, cex.axis=1.2, col='gray')#, axes=FALSE)
#magaxis(side=1, unlog=TRUE, box=TRUE, cex.axis=1.2)
#magaxis(side=2, las=1, cex.axis=1.2)
mtext('alpha', side=2, line=4, cex=1.2)
mtext('q-value', side=1, line=2.5, cex=1.2)
dev.off()

snp.dens<-read.table('~/Dropbox/ddRADseq/var_maf0.1_minmeanDP20_minGQ25_maxmissing0.75_SNPdensity.snpden', header=TRUE)
het<-read.table('~/Dropbox/ddRADseq/var_maf0.1_minmeanDP20_minGQ25_maxmissing0.75_het.het', header=TRUE)
site.qual<-read.table('~/Dropbox/ddRADseq/var_maf0.1_minmeanDP20_minGQ25_maxmissing0.75_sitequal.lqual', header=TRUE)

par(mfrow=c(1,1))
hist(snp.dens$VARIANTS.KB, breaks=seq(0,120,1))
hist(het$F)
hist(site.qual$QUAL)
plot(site.qual$POS, site.qual$QUAL) #all scores are the same

length(unique(site.qual$CHROM))
length(site.qual$CHROM)
