########################################################################################################
## Merge the VCF files created by filtering at the lowest level of data organization (site by species) #
########################################################################################################

#################################
# Load and process the VCF data #
#################################

# load in the filtered VCF data
# the first several hundred thousand lines (hahahaha) are VCF header stuff and not data
# R will load the VCF automatically and skip those lines because they are preceeded by "#", which R takes as a comment symbol
# But... then you lose the header line which has the sample IDs, because that's preceeded by "#" as well
# Instead of making hard edit to the VCF file (a dangerous proposition), use the skip lines argument and set the header as the first non-skipped line
#opHB<-read.table('~/Desktop/UT_ddRADseq/Hierarchical_pseudoref_maxmissing_0.75_minmeanDP20_minGQ25_op_HB_maf0.05.recode.vcf', comment.char="", skip=417187, header=TRUE)
#opSRT<-read.table('~/Desktop/UT_ddRADseq/Hierarchical_pseudoref_maxmissing_0.75_minmeanDP20_minGQ25_op_SRT_maf0.05.recode.vcf', comment.char="", skip=417187, header=TRUE)
#raccHB<-read.table('~/Desktop/UT_ddRADseq/Hierarchical_pseudoref_maxmissing_0.75_minmeanDP20_minGQ25_racc_HB_maf0.05.recode.vcf', comment.char="", skip=417187, header=TRUE)
#raccSRT<-read.table('~/Desktop/UT_ddRADseq/Hierarchical_pseudoref_maxmissing_0.75_minmeanDP20_minGQ25_racc_SRT_maf0.05.recode.vcf', comment.char="", skip=417187, header=TRUE)

# non-hierarchical data (misnamed -- not filtered by site as well)
#opAll<-read.table('~/Desktop/UT_ddRADseq/Hierarchical_pseudoref_counts_maxmissing_0.75_minmeanDP20_minGQ25_op_all_maf0.05.recode.vcf', comment.char="", skip=417187, header=TRUE)
#raccAll<-read.table('~/Desktop/UT_ddRADseq/Hierarchical_pseudoref_counts_maxmissing_0.75_minmeanDP20_minGQ25_racc_all_maf0.05.recode.vcf', comment.char="", skip=417187, header=TRUE)

# non-hierarchical data, HB only, filtered first by whole population (minmeanDP20, minGQ25, maxmissing0.75), 
# then by raccoon/opossum subsets (maf0.05 each)

#opHB<-read.table('~/Desktop/D_variabilis_Pseudoref/opossum_pseudoref_maf0.05_minmeanDP20_minGQ25_maxmissing0.75.recode.vcf', comment.char="", skip=417187, header=TRUE)
#raccHB<-read.table('~/Desktop/D_variabilis_Pseudoref/raccoon_pseudoref_maf0.05_minmeanDP20_minGQ25_maxmissing0.75.recode.vcf', comment.char="", skip=417187, header=TRUE)

## Final HB Only Analysis
#opHB<-read.table('~/Desktop/D_variabilis_Pseudoref/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_HB_only_maxmissing0.75_OPOSSUM.recode.vcf', comment.char="", skip=417187, header=TRUE)
#raccHB<-read.table('~/Desktop/D_variabilis_Pseudoref/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_HB_only_maxmissing0.75_RACCOON.recode.vcf', comment.char="", skip=417187, header=TRUE)

## Final SRT Only Analysis
opSRT<-read.table('~/Desktop/D_variabilis_Pseudoref/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_SRT_only_maxmissing0.75_OPOSSUM.recode.vcf', comment.char="", skip=417187, header=TRUE)
raccSRT<-read.table('~/Desktop/D_variabilis_Pseudoref/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_SRT_only_maxmissing0.75_RACCOON.recode.vcf', comment.char="", skip=417187, header=TRUE)



# Unique identifiers are the "chromosome" (pseudoreference sequence) and the site of the SNP.
# Pseudoreference sequence IDs are non-unique, SNP sites are non-unique, but the combination of the two is unique
# Create columns of Pseudoreference ID + SNP site to uniquely identify each SNP in each dataset
#opHB.sites<-paste(opHB[,1], opHB[,2], sep="_")
#opHB.update<-cbind(opHB, opHB.sites)

opSRT.sites<-paste(opSRT[,1], opSRT[,2], sep="_")
opSRT.update<-cbind(opSRT, opSRT.sites)

#raccHB.sites<-paste(raccHB[,1], raccHB[,2], sep="_")
#raccHB.update<-cbind(raccHB, raccHB.sites)

raccSRT.sites<-paste(raccSRT[,1], raccSRT[,2], sep="_")
raccSRT.update<-cbind(raccSRT, raccSRT.sites)

#opAll.sites<-paste(opAll[,1], opAll[,2], sep='_')
#opAll.update<-cbind(opAll, opAll.sites)
#raccAll.sites<-paste(raccAll[,1], raccAll[,2], sep='_')
#raccAll.update<-cbind(raccAll, raccAll.sites)

# Identify sites shared by all samples. (unnecessary, but a quick check to see how things are working... the merge() function will do this as well)
shared.sites.hierarchical<-intersect(intersect(opHB.update$opHB.sites, opSRT.update$opSRT.sites), intersect(raccHB.update$raccHB.sites, raccSRT.update$raccSRT.sites))
length(shared.sites.hierarchical)
write.table(shared.sites.hierarchical, file='~/Desktop/UT_ddRADseq/hierarchical_all_shared_sites', quote=FALSE, col.names=FALSE, row.names=FALSE)

# Identify sites shared only by samples from HB site. (also unnecessary)
shared.sites.HB<-intersect(opHB.update$opHB.sites, raccHB.update$raccHB.sites)
length(shared.sites.HB)

# Merge data by site
HB.merge<-merge(opHB.update, raccHB.update, by.x="opHB.sites", by.y="raccHB.sites")
SRT.merge<-merge(opSRT.update, raccSRT.update, by.x="opSRT.sites", by.y="raccSRT.sites")

all.data.merge<-merge(HB.merge, SRT.merge, by.x="opHB.sites", by.y="opSRT.sites")

op.racc.merge<-merge(opAll.update, raccAll.update, by.x='opAll.sites', by.y='raccAll.sites')

##########################################
# Prep the data for export as a VCF file #
##########################################

## Full data set
# The merge process leaves us with a few redundant columns. Check column names to identify which columns can be removed
names(all.data.merge)
# Bind together only the necessary columns.
all.data.vcf<-cbind(all.data.merge[,2:39], all.data.merge[,49:102], all.data.merge[,112:113], all.data.merge[,123:136])
# Check that the number of columns with genotype data matches expectation
length(all.data.vcf[1,]) #99 individuals in full data set -- checks out (108-9 leading columns)
# The first 10 columns got renamed... return them to their original state
final.names<-c("#X.CHROM", names(opHB[2:9]), names(all.data.vcf[10:length(names(all.data.vcf))]))
names(all.data.vcf)<-final.names
# How'd that work out?
head(all.data.vcf) #seems ok!
## write it file (output must be tab-separated)
write.table(all.data.vcf, file="~/Desktop/UT_ddRADseq/Hierarchical_Merged_Pseudoref_minmeanDP20_minGQ25_maf0.05.vcf", quote=FALSE,
            row.names=FALSE, sep="\t")

## HB only, same process as above
names(HB.merge)
HB.vcf<-cbind(HB.merge[,2:39], HB.merge[,49:102])
length(HB.vcf[1,]) #83 individuals in HB dataset -- checks out (92-9 leading columns)
HB.names<-sub('\\.', '-', c(names(opHB[2:9]),names(HB.vcf[10:length(names(HB.vcf))]))) #fix the formatting of the names to be consistent with keep_files
final.HB.names<-c("#CHROM", HB.names)
names(HB.vcf)<-final.HB.names
head(HB.vcf)
## write it to file
write.table(HB.vcf, file="~/Desktop/D_variabilis_Pseudoref//Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_HB_only_maxmissing0.75_MERGED.vcf", quote=FALSE,
            row.names=FALSE, sep="\t")

## SRT only, same process as above
names(SRT.merge)
SRT.vcf<-cbind(SRT.merge[,2:12], SRT.merge[,22:35])
length(SRT.vcf[1,]) #16 individuals in SRT dataset -- checks out (92-9 leading columns)
SRT.names<-sub('\\.', '-', c(names(opSRT[2:9]),names(SRT.vcf[10:length(names(SRT.vcf))]))) #fix the formatting of the names to be consistent with keep_files
final.SRT.names<-c("#CHROM", SRT.names)
names(SRT.vcf)<-final.SRT.names
head(SRT.vcf)
## write it to file
write.table(SRT.vcf, file="~/Desktop/D_variabilis_Pseudoref/Final_Analysis/Structure_by_Host_SRT/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_SRT_only_maxmissing0.75_MERGED.vcf", quote=FALSE,
            row.names=FALSE, sep="\t")

## Take the final analysis SRT and HB and merge them together...

## Final HB Only Analysis
HB<-read.table('~/Desktop/UT_ddRADseq/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_HB_only_maxmissing0.75_MERGED.vcf', comment.char="", header=TRUE)

## Final SRT Only Analysis
SRT<-read.table('~/Desktop/D_variabilis_Pseudoref/Final_Analysis/Structure_by_Host_SRT/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_SRT_only_maxmissing0.75_MERGED.vcf', comment.char="", header=TRUE)

HB.sites<-paste(HB[,1], HB[,2], sep="_")
HB.update<-cbind(HB, HB.sites)

SRT.sites<-paste(SRT[,1], SRT[,2], sep="_")
SRT.update<-cbind(SRT, SRT.sites)

all.data.merge<-merge(HB.update, SRT.update, by.x="HB.sites", by.y="SRT.sites")

# The merge process leaves us with a few redundant columns. Check column names to identify which columns can be removed
names(all.data.merge)
# Bind together only the necessary columns.
all.data.vcf<-cbind(all.data.merge[,2:93], all.data.merge[,103:118])
# Check that the number of columns with genotype data matches expectation
length(all.data.vcf[1,]) #99 individuals in full data set -- checks out (108-9 leading columns)
# The first 10 columns got renamed... return them to their original state
final.names<-c("#CHROM", names(HB[2:9]), gsub("[.]", "-", names(all.data.vcf[10:length(names(all.data.vcf))])))
names(all.data.vcf)<-final.names
# How'd that work out?
head(all.data.vcf) #seems ok!
## write it file (output must be tab-separated)
write.table(all.data.vcf, file="~/Desktop/D_variabilis_Pseudoref/Final_Analysis/Structure_by_Site/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_BOTH_SITES_host_filtered_only_maxmissing0.75_MERGED.vcf", quote=FALSE,
            row.names=FALSE, sep="\t")

#---------------------------------------------------------------------
## Final for geographic analysis
## Raccoons and opossums from both sites (not hierarchically filtered)
names(op.racc.merge)
op.racc.vcf<-cbind(op.racc.merge[,2:41], op.racc.merge[,51:118])
length(op.racc.vcf[1,]) #108-9=99 -- checks out
final.op.racc.names<-c("#X.CHROM", names(opAll[2:9]), names(op.racc.vcf[10:length(names(op.racc.vcf))]))
names(op.racc.vcf)<-final.op.racc.names
write.table(op.racc.vcf, file="~/Desktop/UT_ddRADseq/Both_Sites_Racc_Op_Pseudoref_minmeanDP20_minGQ25_maf0.05.vcf", quote=FALSE,
            row.names=FALSE, sep="\t")
