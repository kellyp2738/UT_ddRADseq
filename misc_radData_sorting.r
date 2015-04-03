HB_SRT<-read.table('~/Desktop/distruct1.1/Hierarchical_HB_SRT.indivq')
HB_SRT_sorted<-HB_SRT[order(HB_SRT[,6]),]
write.table(HB_SRT_sorted, file="~/Desktop/distruct1.1/Hierarchical_HB_SRT_sorted.indivq", sep=" ",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

racc_op<-read.table('~/Desktop/distruct1.1/Hierarchical_racc_op.indivq')
racc_op_sorted<-racc_op[order(racc_op[,6]),]
write.table(racc_op_sorted, file="~/Desktop/distruct1.1/Hierarchical_racc_op_sorted.indivq", sep=" ",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

racc_op_both_sites<-read.table('~/Desktop/distruct1.1/ninth_try.indivq')
racc_op_both_sites_sorted<-racc_op_both_sites[order(racc_op_both_sites[,6]),]
write.table(racc_op_both_sites_sorted, file="~/Desktop/distruct1.1/ninth_try_sorted.indivq", sep=" ",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

HB.popfilter<-read.table('~/Desktop/distruct1.1/tenth_try_by_sp.indivq')
popfilter.sorted<-HB.popfilter[order(HB.popfilter[,6]),]
write.table(popfilter.sorted, file="~/Desktop/distruct1.1/tenth_try_by_sp_sorted.indivq", sep=" ",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

HB.popfilter.vcf<-read.table('~/Desktop/console/tenth_try/infile')
HB.vcf.sort<-HB.popfilter.vcf[order(HB.popfilter.vcf[,1]),]
write.table(HB.vcf.sort, file='~/Desktop/console/tenth_try/infile_sorted', sep=" ",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

final<-read.table('~/Desktop/D_variabilis_Pseudoref/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_HB_only_maxmissing0.75_MERGED_STRUCTURE')
final.sort<-final[order(final[,1]),]
write.table(final.sort, file='~/Desktop/D_variabilis_Pseudoref/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_HB_only_maxmissing0.75_MERGED_STRUCTURE_sorted',
            quote=FALSE, row.names=FALSE, col.names=FALSE)
  
indivq<-read.table('~/Desktop/distruct1.1/final_by_sp.indivq')
indivq.sorted<-indivq[order(indivq[,4], indivq[,6]),]
write.table(indivq.sorted, file='~/Desktop/distruct1.1/final_by_sp_sorted.indivq', sep=" ", quote=FALSE, row.names=FALSE, col.names=FALSE)
