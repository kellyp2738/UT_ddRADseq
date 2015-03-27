HB<-read.table('~/Desktop/UT_ddRADseq/Updated_pseudoref_maxmissing_0.75_minmeanDP20_minGQ25_HB_maf0.1.recode.vcf')
SRT<-read.table('~/Desktop/UT_ddRADseq/Updated_pseudoref_maxmissing_0.75_minmeanDP20_minGQ25_SRT_maf0.1.recode.vcf')

HB.sites<-paste(HB[,1], HB[,2], sep="_")
HB.update<-cbind(HB, HB.sites)

SRT.sites<-paste(SRT[,1], SRT[,2], sep="_")
SRT.update<-cbind(SRT, SRT.sites)

shared.sites<-intersect(HB.update$HB.sites, SRT.update$SRT.sites)

HB.SRT.merge<-merge(HB.update, SRT.update, by.x="HB.sites", by.y="SRT.sites")
