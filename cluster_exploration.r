merged<-read.table('~/Desktop/UT_ddRADseq/Updated_pseudoref_maxmissing_0.75_minmeanDP20_minGQ25_HB_maf0.1_STRUCTURE')
merged.infile<-read.table('~/Desktop/console/sixth_try/infile', na.strings='-9')
hb<-read.table('~/Desktop/UT_ddRADseq/Updated_pseudoref_HB_SRT_merge_maxmissing_0.75_minmeanDP20_minGQ25_STRUCTURE')
merged.raw<-read.table('~/Desktop/UT_ddRADseq/Updated_pseudoref_HB_SRT_merge_maxmissing_0.75_minmeanDP20_minGQ25.vcf')
dim(merged.raw)
head(unique(merged.raw[,1]))
head(merged.raw)

head(merged.infile)
complete<-merged.infile[complete.cases(merged.infile),]

hb.infile<-subset(merged.infile, merged.infile[,2]=='1')
srt.infile<-subset(merged.infile, merged.infile[,2]=='2')

K3.data<-read.table('~/Desktop/distruct1.1/seventh_try_tick_K3.indivq')
K3.names<-c()
for(n in 1:length(K3.data[,2])){
  K3.names<-c(K3.names, strsplit(as.character(K3.data[n,2]), '_')[[1]][1])
}
K3.data<-cbind(K3.data, K3.names)

sample.info<-read.csv('~/Dropbox/Raw Data/Vertical_Movement_and_other_CLNWR_2012/ddRAD_FinalLibrary_SampleInfo_Full.csv')

augmented.K3<-merge(K3.data, sample.info, by.x='K3.names', by.y='tick.id')
names(augmented.K3)<-c('tick.id', 'K3.row.num', 'truncated.id', 'missing', 'a.priori.pop', '', 
                       'cluster1', 'cluster2', 'cluster3', 'inline.barcode', 'Illlumina.barcode', 'host.sp',
                       'host.id', 'coll.site', 'trap.id', 'hb.roadside', 'combinatoric.id')
plot(augmented.K3$trap.id, augmented.K3$cluster1)
plot(augmented.K3$host.id, augmented.K3$cluster1)
