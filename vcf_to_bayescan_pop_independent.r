## -----------------------------------------------------------------------------------
## -- PARSE ALLELE COUNTS FILE FOR BAYESCAN
## -----------------------------------------------------------------------------------

# The algorithm here is much simpler than the one in vcf_to_bayescan.r because it
# assumes that the input data is a merged VCF file made from SNP sites that are highly
# represented in the dataset (no more than 25% missing data per site, and sites present
# in both raccoon and opossum datasets). This filtering makes it highly likely that data
# subsetted by host individual will also meet the same requirement that all sites are 
# present in each subpopulation (but this is not necessarily true). The vcf_to_baysecan.r
# file doesn't assume all sites are present in subpopulation allele count files, and so 
# checks to make sure each population allele count tabulation contains the same number of 
# sites and removes sites that aren't present in both populations. 

# Interactively run the script and CHECK THAT EACH POPULATION ALLELE TABULATION contains the
# same number of sites, and that the number of sites is equal to that from the original VCF
# file. If this is the case, all is well. If not, this script needs to be extended to return
# only sites shared across all populations.

library(stringr)
library(plyr)
library(magicaxis)
library(gplots)
#source('~/BayeScan2.1/R functions/plot_R.r') # not installed locally

setwd('~/Desktop/UT_ddRADseq/By_Indv_Counts/') #path to allele count files, housed in their own unique directory
pop.files<-list.files(getwd())

alleles.list<-list() #empty container for parsed alleles files; each entry is the parsed allele file for a different subpopulation
alleles.save<-list()

for(n in 1:length(pop.files)){ #iteratively process files
  pop<-read.table(file.path(getwd(), pop.files[n]), skip=1)
  names(pop)<-c('CHROM', 'POS', 'N_ALLELES', 'N_CHR', 'ALLELE.COUNT.1', 'ALLELE.COUNT.2')
  
  # parse alleles
  pop.1.allele.1<-unlist(strsplit(as.character(pop$ALLELE.COUNT.1), split=':'))
  pop.1.allele.1.ID<-pop.1.allele.1[seq(1,length(pop.1.allele.1),2)]
  pop.1.allele.1.count<-pop.1.allele.1[seq(2,length(pop.1.allele.1),2)]
  pop.1.allele.1.parsed<-cbind(pop.1.allele.1.ID, pop.1.allele.1.count)
  
  pop.1.allele.2<-unlist(strsplit(as.character(pop$ALLELE.COUNT.2), split=':'))
  pop.1.allele.2.ID<-pop.1.allele.2[seq(1,length(pop.1.allele.2),2)]
  pop.1.allele.2.count<-pop.1.allele.2[seq(2,length(pop.1.allele.2),2)]
  pop.1.allele.2.parsed<-cbind(pop.1.allele.2.ID, pop.1.allele.2.count)
  
  pop.1.parsed<-cbind(pop[,1:4], pop.1.allele.1.parsed, pop.1.allele.2.parsed)
  pop.save<-cbind(seq(1:length(pop.1.parsed[,1])), pop.1.parsed[,4], pop.1.parsed[,3], pop.1.parsed[,6], pop.1.parsed[,8])
  # store the parsed alleles in a list
  alleles.list[[n]]<-pop.1.parsed
  alleles.save[[n]]<-pop.save
  print(c('Number of Sites', pop.files[n], length(pop.save[,1])))
  write.table(pop.save, file=file.path('~/Dropbox/ddRADseq/', paste(pop.files[n], '_Bayescan', sep="")), quote=FALSE, row.names=FALSE, col.names=FALSE) 
}

by.ind<-read.table('~/Dropbox/ddRADseq/D_variabilis_Pseudoref/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_HB_only_maxmissing0.75_MERGED_BY_INDIVIDUAL_count_Baye_fst.txt')
plot(by.ind$qval, by.ind$fst)