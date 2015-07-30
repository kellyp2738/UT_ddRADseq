# trying to estimate linkage disequilibrium..


library('genetics')
library('adegenet')
library('fdrtool')

# data in genlight format
site<-read.snp('~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Site/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_BOTH_SITES_host_filtered_only_maxmissing0.75_MERGED_FOR_ADEGENET.finalSNPs_site.snp')

# the original VCF 012 output file is better suited for genind format
# it needs some translation though in bash, as follows:

# from before, remove the first column and change -1 to - as the NA character (see prep_for_adegenet.md for details)
# sed 's|0|AA|g' Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_BOTH_SITES_host_filtered_only_maxmissing0.75_MERGED_FOR_ADEGENET.rmCol1.012 > Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_BOTH_SITES_host_filtered_only_maxmissing0.75_MERGED_FOR_ADEGENET.AA.012
# sed 's|1|AB|g' Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_BOTH_SITES_host_filtered_only_maxmissing0.75_MERGED_FOR_ADEGENET.AA.012 > Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_BOTH_SITES_host_filtered_only_maxmissing0.75_MERGED_FOR_ADEGENET.AB.012
# sed 's|2|BB|g' Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_BOTH_SITES_host_filtered_only_maxmissing0.75_MERGED_FOR_ADEGENET.AB.012 > Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_BOTH_SITES_host_filtered_only_maxmissing0.75_MERGED_FOR_ADEGENET.BB.012

site.012<-data.frame(read.table('~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Site/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_BOTH_SITES_host_filtered_only_maxmissing0.75_MERGED_FOR_ADEGENET.BB.012',
                     na.strings="-"))

# the original VCF 012 output and the genlight data have the same sample order
# use the genlight sample IDs to inform creation of the genind dataset
site.genind<-df2genind(site.012, ind.names=site@ind.names,
                       missing='<NA>', ploidy=2)

# convert to genotype format for LD calculation in genetics package
site.geno<-genind2genotype(site.genind)

site.LD<-LD(site.geno)
save(site.LD, file="~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Site/LD_estimate.txt")
write(site.LD, file="~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Site/LD_extimate_text.txt")
load("~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Site/LD_estimate.txt")

plot(site.LD)

summary(site.LD)

# calculate FDR for p-values
hist(site.LD$'R^2')
hist(site.LD$'P-value')
plot(site.LD$'R^2', site.LD$'P-value')
sig.high.R2<-subset

# p-values are in the upper triangle of a matrix
p.only<-upper.tri(site.LD$'P-value')
fdr.corrected<-fdrtool(site.LD$'P-value'[p.only], statistic=c('pvalue'))

# column and row names are SNP IDs
#pvals<-fdr.corrected$pval
qvals<-fdr.corrected$qval
fraction.sig<-length(qvals[qvals<0.05])/length(qvals)

# put the qval data back into a matrix
q.mat<-matrix(data=NA, nrow=site$n.loc, ncol=site$n.loc)
q.mat[upper.tri(q.mat)]<-qvals

plot(site.LD$'R^2', q.mat)

# comparing the upper triangle of this newly created matrix with the upper
# triangle of the site.LD$'P-value' shows that the data are in the same order
# therefore the column and row names for the new q-value matrix accurately reflect
# the IDs of the SNP pairs and can be used to match q-values with actual genotype
# data from the original "site" dataset.
#test.mat<-matrix(data=NA, nrow=2387, ncol=2387)
#test.mat[upper.tri(test.mat)]<-pvals

rm.IDs<-c()
rm.IDs.smart<-c()
for(i in 1:site$n.loc){
  for(j in 1:site$n.loc){
    if(!is.na(q.mat[i,j])){
      if(q.mat[i,j]<0.05){
        rm.IDs<-c(rm.IDs, i, j) # IDs may be added multiple times if they are in LD with several other IDs
        if(!(i %in% rm.IDs.smart)){
          if(!(j %in% rm.IDs.smart)){
            rm.IDs.smart<-c(rm.IDs.smart, sample(c(i,j), 1)) # randomly choose one member of the pair to remove
            # if either i or j is already in the list, our work here is done...
          }
        }
      }
    }
  }
}

length(unique(rm.IDs)) # we see that ever locus is in significant LD with at least one other locus
rm.IDs.table<-table(rm.IDs)

length(rm.IDs.smart)
rm.IDs.smart.table<-table(rm.IDs.smart)

length(unique(rm.IDs))-length(rm.IDs.smart)

# if we were to calculate the significance of R^2 by hand...

r2<-seq(0,1,0.01) #possible values of R^2
n<-2387 #let's say N=number of SNPs, not N=number of chromosomes
chisq.stat<-r2*n
sig<-dchisq(chisq.stat, df=1) #2 alleles - 1 = 1 df
plot(r2, sig) # all of these values are significant... the N here is incorrect