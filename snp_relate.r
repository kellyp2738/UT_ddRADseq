## SNPRelate Analysis on HB Only Data

#####################################################################
## Install packages

source("http://bioconductor.org/biocLite.R")
biocLite("gdsfmt")
biocLite("SNPRelate")

# netcdf is a dependency that R won't be able to automatically download
# install Ubuntu dependencies first, then install with argument to point to the correct directory for the header file
#sudo apt-get install libnetcdf-dev
#sudo apt-get install udunits-bin
#sudo apt-get install libudunits2-dev
install.packages("ncdf", type = "source", configure.args="--with-netcdf-include=/usr/include")
install.packages("ncdf") # works just fine like this for OSX
biocLite("GWASTools")
#####################################################################

#####################################################################
## Load packages

library(GWASTools)
library(SNPRelate)
library(gplots)
library(fields)
library(RColorBrewer)
#####################################################################

#####################################################################
## BOTH SITES 
##

extra.data<-read.csv('~/Desktop/UT_ddRADseq/ddRAD_FinalLibrary_SampleInfo_Full.csv', header=TRUE)

## Load VCF file and convert to GDS object
## following instructions found here: http://www.bioconductor.org/packages/release/bioc/vignettes/GWASTools/inst/doc/Formats.pdf
## this site also has great documentation on getting your data into SNPRelate: http://corearray.sourceforge.net/tutorials/SNPRelate/#installation-of-the-package-snprelate

## format chromosome IDs properly
vcf.all.fix<-read.table('~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Site/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_BOTH_SITES_host_filtered_only_maxmissing0.75_MERGED.vcf', 
                    comment.char="", header=TRUE)
vcf.all.fix[,1]<-gsub(pattern="_*[a-z]*", replacement="", vcf.all.fix[,1], perl=TRUE)
write.table(vcf.all.fix, file='~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Site/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_BOTH_SITES_host_filtered_maxmissing0.75_MERGED_chromFix.vcf',
            quote=FALSE, row.names=FALSE, sep="\t")
## manual edits to the VCF file so it can be loaded:
##    - change X.CHROM to #CHROM
##    - replace "." with "-" in the sample names (prob unnecessary, but it was done)

## read in edited file
gds.all.ticks<-"snps_all_again2.gds"
snpgdsVCF2GDS('~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Site/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_BOTH_SITES_host_filtered_maxmissing0.75_MERGED_chromFix.vcf', 
              gds.all.ticks, verbose=TRUE)
(gds<-GdsGenotypeReader(gds.all.ticks))
getScanID(gds) #check IDs

snpID <- getSnpID(gds)
chromosome <- as.integer(getChromosome(gds))
position <- getPosition(gds)
alleleA <- getAlleleA(gds)
alleleB <- getAlleleB(gds)
rsID <- getVariable(gds, "snp.rs.id")
qual <- getVariable(gds, "snp.annot/qual")
filter <- getVariable(gds, "snp.annot/filter")
snpAnnot <- SnpAnnotationDataFrame(data.frame(snpID, chromosome, position,
                                              rsID, alleleA, alleleB,
                                              qual, filter, stringsAsFactors=FALSE))

#scanID<-as.character(extra.data$combo.label) doesn't work... doesn't stay as character
does.it.match<-cbind(sort(getScanID(gds)), as.character(extra.data$combo.label)) #if extra.data and getScanID(gds) are in the same order, just pull the scanIDs directly from gds
pop.group<-as.character(extra.data$coll.site)
sex<-as.character(extra.data$sex)
scanAnnot<-ScanAnnotationDataFrame(data.frame(scanID=getScanID(gds), pop.group, sex))
genoData<-GenotypeData(gds, snpAnnot=snpAnnot, scanAnnot=scanAnnot)
getGenotype(genoData)
getSex(genoData)
getScanVariable(genoData, "pop.group")

createDataFile(path=".", filename="fixed_all_snps", file.type="gds",
               variables="genotype", snp.annotation=snpAnnot, scan.annotation=scanAnnot)

snpgdsCombineGeno(c(genoData), "snp_output")
snpgdsSummary(genoData)

close(genoData)
open(genoData)


## reopen data
snpgdsClose(genofile)
genofile<-snpgdsOpen('snps_all_again2.gds', allow.duplicate=TRUE) #in /Users/kelly

#genofile <- snpgdsOpen(snpgdsExampleFileName())

get.attr.gdsn(index.gdsn(genofile))
genofile<-GdsGenotypeReader('snps_all.gds')

## Fst
getSex(genofile)
pop_code<-read.gdsn(index.gdsn(genofile, path="sample.annot/pop.group"))

genofile$scanAnnot

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pop_code <- read.gdsn(index.gdsn(genofile, "sampleAnnot"))
attributes(genofile)
genofile$names

snpgdsFst(genofile, population=getScanVariable(genoData, "pop.group"), method="W&C84")

## Identity by descent (IBD)

#not LD pruned
ibd<-snpgdsIBDMLE(genofile, autosome.only=FALSE, kinship=TRUE) 
ibd.coeff<-snpgdsIBDSelection(ibd)
write.table(ibd.coeff, file='~/Desktop/D_variabilis_Pseudoref/identity_by_descent_HB.txt',
            quote=FALSE, row.names=FALSE) #since moved; see file path below
ibd.coeff<-read.table(file='~/Desktop/')
plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1))
plot(ibd$kinship, pch=16, col=alpha('blue', 0.25), xlim=c(0,1), ylim=c(0,1))
related<-ibd.coeff[which(ibd.coeff$kinship>=1/16),]
related.unique<-unique(c(related$ID1, related$ID2))
ibd.thresholds<-c(0,1/16,1/8,1/4,1/2) #color by relatedness (0-1/16=unrelated, 1/16-1/8=first cousins, 1/8-1/4=half sibs or dbl first cousisn, 1/4-1/2=parent/child or full sibs, 1/2=identical)
tick.names<-gsub(pattern="_.*", replacement="", unique(c(ibd.coeff$ID1, ibd.coeff$ID2)), perl=TRUE)
#pdf(file='~/Desktop/D_variabilis_Pseudoref/kinship_estimates_IBD_HB.pdf', height=30/2.54, width=30/2.54)
par(mar=c(6,6,4,10), xpd=TRUE)
image(x=1:ncol(ibd$kinship), y=1:ncol(ibd$kinship), z=ibd$kinship, col=brewer.pal(4, 'Blues'), 
      breaks=ibd.thresholds, axes=FALSE, xlab="", ylab="")
mtext(side=1, line=1, text=tick.names, at=1:ncol(ibd$kinship), las=2, cex=0.5)
mtext(side=2, line=1, text=tick.names, at=1:ncol(ibd$kinship), las=2, cex=0.5)
legend(x=85, y=50, , fill=brewer.pal(4, 'Blues'),
       legend=c('0-1/16', '1/16-1/8', '1/8-1/4', '1/4-1/2'),
       #legend=c('unrelated', 'first cousins', 'half-siblings', 'parent-child or full siblings', 'identical'),
       bty='n', border='white', title='Kinship \nCoefficient')

lines(x=c(0.5,8.5), y=c(0.5, 0.5))
lines(y=c(0.5,8.5), x=c(0.5, 0.5))
lines(x=c(0.5,8.5), y=c(8.5, 8.5))
lines(y=c(0.5,8.5), x=c(8.5, 8.5))

lines(x=c(29.5,29.5), y=c(8.5, 29.5))
lines(y=c(29.5,29.5), x=c(8.5, 29.5))
lines(y=c(8.5,8.5), x=c(8.5, 29.5))
lines(x=c(8.5,8.5), y=c(8.5, 29.5))

lines(x=c(29.5,29.5), y=c(29.5, 39.5))
lines(y=c(29.5,29.5), x=c(29.5, 39.5))
lines(x=c(39.5,39.5), y=c(29.5, 39.5))
lines(y=c(39.5,39.5), x=c(29.5, 39.5))

lines(x=c(39.5,39.5), y=c(39.5, 43.5))
lines(y=c(39.5,39.5), x=c(39.5, 43.5))
lines(x=c(43.5,43.5), y=c(39.5, 43.5))
lines(y=c(43.5,43.5), x=c(39.5, 43.5))

lines(x=c(43.5,43.5), y=c(43.5, 61.5))
lines(y=c(43.5,43.5), x=c(43.5, 61.5))
lines(x=c(61.5,61.5), y=c(43.5, 61.5))
lines(y=c(61.5,61.5), x=c(43.5, 61.5))

lines(x=c(61.5,61.5), y=c(61.5, 65.5))
lines(y=c(61.5,61.5), x=c(61.5, 65.5))
lines(x=c(65.5,65.5), y=c(61.5, 65.5))
lines(y=c(65.5,65.5), x=c(61.5, 65.5))

lines(x=c(83.5,83.5), y=c(65.5, 83.5))
lines(y=c(83.5,83.5), x=c(65.5, 83.5))
lines(x=c(65.5,65.5), y=c(65.5, 83.5))
lines(y=c(65.5,65.5), x=c(65.5, 83.5))

dev.off()

# write out a list of related individuals to exclude from future analyses
#write.table(related.unique, file='~/Desktop/UT_ddRADseq/keep_files/exclude_related.txt', 
#            quote=FALSE, row.names=FALSE, col.names=FALSE)




#####################################################################
## HARRISON BAYOU ONLY 
##
## Load VCF file and convert to GDS object
## following instructions found here: http://www.bioconductor.org/packages/release/bioc/vignettes/GWASTools/inst/doc/Formats.pdf
## this site also has great documentation on getting your data into SNPRelate: http://corearray.sourceforge.net/tutorials/SNPRelate/#installation-of-the-package-snprelate

## format chromosome IDs properly
vcf.fix<-read.table('/home/antolinlab/Desktop/D_variabilis_Pseudoref/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_HB_only_maxmissing0.75_MERGED.vcf', 
                    comment.char="", header=TRUE)
vcf.fix[,1]<-gsub(pattern="_*[a-z]*", replacement="", vcf.fix[,1], perl=TRUE)
write.table(vcf.fix, file='/home/antolinlab/Desktop/D_variabilis_Pseudoref/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_HB_only_maxmissing0.75_MERGED_chromFix.vcf',
            quote=FALSE, row.names=FALSE, sep="\t")
## manual edits to the VCF file so it can be loaded:
##    - change X.CHROM to #CHROM
##    - replace "." with "-" in the sample names (prob unnecessary, but it was done)

## read in edited file
gds.ticks<-"snps_fixed.gds"
#snpgdsVCF2GDS('/home/antolinlab/Desktop/D_variabilis_Pseudoref/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_HB_only_maxmissing0.75_MERGED_chromFix.vcf', 
#              gds.ticks, verbose=FALSE)
snpgdsVCF2GDS('~/Dropbox/ddRADseq/Final_Analysis/Structure_by_Site/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_HB_only_maxmissing0.75_MERGED_chromFix.vcf', 
              gds.ticks, verbose=FALSE)
(gds<-GdsGenotypeReader(gds.ticks))
getScanID(gds) #check IDs

snpID <- getSnpID(gds)
chromosome <- as.integer(getChromosome(gds))
position <- getPosition(gds)
alleleA <- getAlleleA(gds)
alleleB <- getAlleleB(gds)
rsID <- getVariable(gds, "snp.rs.id")
qual <- getVariable(gds, "snp.annot/qual")
filter <- getVariable(gds, "snp.annot/filter")
snpAnnot <- SnpAnnotationDataFrame(data.frame(snpID, chromosome, position,
                                              rsID, alleleA, alleleB,
                                              qual, filter, stringsAsFactors=FALSE))
genoData<-GenotypeData(gds, snpAnnot=snpAnnot)
getGenotype(genoData)
close(genoData)

## reopen data
genofile<-snpgdsOpen('snps_fixed.gds', allow.duplicate=TRUE) #in /home/antolinlab

## LD pruning
snpset<-snpgdsLDpruning(genofile, ld.threshold=0.1, autosome.only=FALSE)

## PCA results are really consistent with adegenet analysis on full dataset
snpset.id<-unlist(snpset)
pca<-snpgdsPCA(genofile, snp.id=snpset.id, autosome.only=FALSE)
pc.percent<-pca$varprop*100
head(round(pc.percent, 2))
tab<-data.frame(sample.id=pca$sample.id, EV1=pca$eigenvect[,1],
                EV2=pca$eigenvect[,2], stringsAsFactors=FALSE)
plot(tab$EV2, tab$EV1)

## Identity by descent (IBD)

#not LD pruned
ibd<-snpgdsIBDMLE(genofile, autosome.only=FALSE, kinship=TRUE) 
ibd.coeff<-snpgdsIBDSelection(ibd)
write.table(ibd.coeff, file='~/Desktop/D_variabilis_Pseudoref/identity_by_descent_HB.txt',
            quote=FALSE, row.names=FALSE) #since moved; see file path below
ibd.coeff<-read.table(file='~/Desktop/')
plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1))
plot(ibd$kinship, pch=16, col=alpha('blue', 0.25), xlim=c(0,1), ylim=c(0,1))
related<-ibd.coeff[which(ibd.coeff$kinship>=1/16),]
related.unique<-unique(c(related$ID1, related$ID2))
ibd.thresholds<-c(0,1/16,1/8,1/4,) #color by relatedness (0-1/16=unrelated, 1/16-1/8=first cousins, 1/8-1/4=half sibs or dbl first cousisn, 1/4-1/2=parent/child or full sibs, 1/2=identical)
tick.names<-gsub(pattern="_.*", replacement="", unique(c(ibd.coeff$ID1, ibd.coeff$ID2)), perl=TRUE)
pdf(file='~/Desktop/D_variabilis_Pseudoref/kinship_estimates_IBD_HB.pdf', height=30/2.54, width=30/2.54)
par(mar=c(6,6,4,10), xpd=TRUE)
image(x=1:ncol(ibd$kinship), y=1:ncol(ibd$kinship), z=ibd$kinship, col=brewer.pal(4, 'Blues'), 
      breaks=ibd.thresholds, axes=FALSE, xlab="", ylab="")
mtext(side=1, line=1, text=tick.names, at=1:ncol(ibd$kinship), las=2, cex=0.5)
mtext(side=2, line=1, text=tick.names, at=1:ncol(ibd$kinship), las=2, cex=0.5)
legend(x=85, y=50, , fill=brewer.pal(4, 'Blues'),
       legend=c('0-1/16', '1/16-1/8', '1/8-1/4', '1/4-1/2'),
       #legend=c('unrelated', 'first cousins', 'half-siblings', 'parent-child or full siblings', 'identical'),
       bty='n', border='white', title='Kinship \nCoefficient')

lines(x=c(0.5,8.5), y=c(0.5, 0.5))
lines(y=c(0.5,8.5), x=c(0.5, 0.5))
lines(x=c(0.5,8.5), y=c(8.5, 8.5))
lines(y=c(0.5,8.5), x=c(8.5, 8.5))

lines(x=c(29.5,29.5), y=c(8.5, 29.5))
lines(y=c(29.5,29.5), x=c(8.5, 29.5))
lines(y=c(8.5,8.5), x=c(8.5, 29.5))
lines(x=c(8.5,8.5), y=c(8.5, 29.5))

lines(x=c(29.5,29.5), y=c(29.5, 39.5))
lines(y=c(29.5,29.5), x=c(29.5, 39.5))
lines(x=c(39.5,39.5), y=c(29.5, 39.5))
lines(y=c(39.5,39.5), x=c(29.5, 39.5))

lines(x=c(39.5,39.5), y=c(39.5, 43.5))
lines(y=c(39.5,39.5), x=c(39.5, 43.5))
lines(x=c(43.5,43.5), y=c(39.5, 43.5))
lines(y=c(43.5,43.5), x=c(39.5, 43.5))

lines(x=c(43.5,43.5), y=c(43.5, 61.5))
lines(y=c(43.5,43.5), x=c(43.5, 61.5))
lines(x=c(61.5,61.5), y=c(43.5, 61.5))
lines(y=c(61.5,61.5), x=c(43.5, 61.5))

lines(x=c(61.5,61.5), y=c(61.5, 65.5))
lines(y=c(61.5,61.5), x=c(61.5, 65.5))
lines(x=c(65.5,65.5), y=c(61.5, 65.5))
lines(y=c(65.5,65.5), x=c(61.5, 65.5))

lines(x=c(83.5,83.5), y=c(65.5, 83.5))
lines(y=c(83.5,83.5), x=c(65.5, 83.5))
lines(x=c(65.5,65.5), y=c(65.5, 83.5))
lines(y=c(65.5,65.5), x=c(65.5, 83.5))

dev.off()

# write out a list of related individuals to exclude from future analyses
write.table(related.unique, file='~/Desktop/UT_ddRADseq/keep_files/exclude_related.txt', 
            quote=FALSE, row.names=FALSE, col.names=FALSE)

#LD pruned
ibd.noLD<-snpgdsIBDMLE(genofile, snp.id=snpset.id, autosome.only=FALSE, kinship=TRUE)
ibd.noLD.coeff<-snpgdsIBDSelection(ibd.noLD)
plot(ibd.noLD.coeff$k0, ibd.noLD.coeff$k1, xlim=c(0,1), ylim=c(0,1))
plot(ibd.noLD$kinship, pch=16, col=alpha('blue', 0.5), xlim=c(0,1), ylim=c(0,1))


#####################################################################
## STARR RANCH TRAIL ONLY 
##
## Load VCF file and convert to GDS object
## following instructions found here: http://www.bioconductor.org/packages/release/bioc/vignettes/GWASTools/inst/doc/Formats.pdf
## this site also has great documentation on getting your data into SNPRelate: http://corearray.sourceforge.net/tutorials/SNPRelate/#installation-of-the-package-snprelate

## format chromosome IDs properly
vcf.SRT.fix<-read.table('/home/antolinlab/Desktop/D_variabilis_Pseudoref/Final_Analysis/Structure_by_Host_SRT/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_SRT_only_maxmissing0.75_MERGED.vcf', 
                    comment.char="", header=TRUE)
vcf.SRT.fix[,1]<-gsub(pattern="_*[a-z]*", replacement="", vcf.SRT.fix[,1], perl=TRUE)
write.table(vcf.SRT.fix, file='/home/antolinlab/Desktop/D_variabilis_Pseudoref/Final_Analysis/Structure_by_Host_SRT/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_SRT_only_maxmissing0.75_MERGED_chromFix.vcf',
            quote=FALSE, row.names=FALSE, sep="\t")
## manual edits to the VCF file so it can be loaded:
##    - change X.CHROM to #CHROM
##    - replace "." with "-" in the sample names (prob unnecessary, but it was done)

## read in edited file
gds.SRT.ticks<-"snps_SRT.gds"
snpgdsVCF2GDS('/home/antolinlab/Desktop/D_variabilis_Pseudoref/Final_Analysis/Structure_by_Host_SRT/Final_Pseudoref_minmeanDP20_minGQ25_maf0.05_SRT_only_maxmissing0.75_MERGED_chromFix.vcf', 
              gds.SRT.ticks, verbose=FALSE)
(gds<-GdsGenotypeReader(gds.SRT.ticks))
getScanID(gds) #check IDs

snpID <- getSnpID(gds)
chromosome <- as.integer(getChromosome(gds))
position <- getPosition(gds)
alleleA <- getAlleleA(gds)
alleleB <- getAlleleB(gds)
rsID <- getVariable(gds, "snp.rs.id")
qual <- getVariable(gds, "snp.annot/qual")
filter <- getVariable(gds, "snp.annot/filter")
snpAnnot <- SnpAnnotationDataFrame(data.frame(snpID, chromosome, position,
                                              rsID, alleleA, alleleB,
                                              qual, filter, stringsAsFactors=FALSE))
genoData<-GenotypeData(gds, snpAnnot=snpAnnot)
getGenotype(genoData)
close(genoData)

## reopen data
genofile<-snpgdsOpen('snps_SRT.gds', allow.duplicate=TRUE) #in /home/antolinlab

## LD pruning
snpset<-snpgdsLDpruning(genofile, ld.threshold=0.1, autosome.only=FALSE)

## PCA results are really consistent with adegenet analysis on full dataset
snpset.id<-unlist(snpset)
pca<-snpgdsPCA(genofile, snp.id=snpset.id, autosome.only=FALSE)
pc.percent<-pca$varprop*100
head(round(pc.percent, 2))
tab<-data.frame(sample.id=pca$sample.id, EV1=pca$eigenvect[,1],
                EV2=pca$eigenvect[,2], stringsAsFactors=FALSE)
plot(tab$EV2, tab$EV1)

## Identity by descent (IBD)

#not LD pruned
ibd<-snpgdsIBDMLE(genofile, autosome.only=FALSE, kinship=TRUE) 
ibd.coeff<-snpgdsIBDSelection(ibd)
write.table(ibd.coeff, file='~/Desktop/D_variabilis_Pseudoref/Final_Analysis/Structure_by_Host_SRT/identity_by_descent_SRT.txt',
            quote=FALSE, row.names=FALSE)
ibd.coeff<-read.table(file='~/Desktop/D_variabilis_Pseudoref/Final_Analysis/Structure_by_Host_SRT/identity_by_descent_SRT.txt', header=TRUE)
plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1))
plot(ibd$kinship, pch=16, col=alpha('blue', 0.25), xlim=c(0,1), ylim=c(0,1))
related<-ibd.coeff[which(ibd.coeff$kinship>=1/16),]
related.unique<-unique(c(related$ID1, related$ID2))
ibd.thresholds<-c(0,1/16,1/8,1/4,1/2) #color by relatedness (0-1/16=unrelated, 1/16-1/8=first cousins, 1/8-1/4=half sibs or dbl first cousisn, 1/4-1/2=parent/child or full sibs, 1/2=identical)
tick.names<-gsub(pattern="_.*", replacement="", ibd$sample.id, perl=TRUE) #unique(c(as.character(ibd.coeff$ID1, ibd.coeff$ID2)))
pdf(file='~/Dropbox/ddRADseq/Final_Plots_April_2015/kinship_estimates_IBD_SRT.pdf', height=30/2.54, width=30/2.54)
par(mar=c(6,6,4,10), xpd=TRUE)
image(x=1:ncol(ibd$kinship), y=1:ncol(ibd$kinship), z=ibd$kinship, col=brewer.pal(4, 'Blues'), 
      breaks=ibd.thresholds, axes=FALSE, xlab="", ylab="")
mtext(side=1, line=1, text=tick.names, at=1:ncol(ibd$kinship), las=2)
mtext(side=2, line=1, text=tick.names, at=1:ncol(ibd$kinship), las=2)
legend(x=17, y=10, , fill=brewer.pal(4, 'Blues'),
       legend=c('0-1/16', '1/16-1/8', '1/8-1/4', '1/4-1/2'),
       #legend=c('unrelated', 'first cousins', 'half-siblings', 'parent-child or full siblings', 'identical'),
       bty='n', border='white', title='Kinship \nCoefficient')

lines(x=c(0.5,2.5), y=c(0.5, 0.5))
lines(y=c(0.5,2.5), x=c(0.5, 0.5))
lines(x=c(0.5,2.5), y=c(2.5, 2.5))
lines(y=c(0.5,2.5), x=c(2.5, 2.5))

lines(x=c(2.5,2.5), y=c(2.5, 3.5))
lines(y=c(2.5,2.5), x=c(2.5, 3.5))
lines(y=c(3.5,3.5), x=c(2.5, 3.5))
lines(x=c(3.5,3.5), y=c(2.5, 3.5))

lines(x=c(4.5,4.5), y=c(4.5, 3.5))
lines(y=c(4.5,4.5), x=c(4.5, 3.5))
lines(y=c(3.5,3.5), x=c(4.5, 3.5))
lines(x=c(3.5,3.5), y=c(4.5, 3.5))

lines(x=c(4.5,4.5), y=c(4.5, 6.5))
lines(y=c(4.5,4.5), x=c(4.5, 6.5))
lines(y=c(6.5,6.5), x=c(4.5, 6.5))
lines(x=c(6.5,6.5), y=c(4.5, 6.5))

lines(x=c(9.5,9.5), y=c(9.5, 6.5))
lines(y=c(9.5,9.5), x=c(9.5, 6.5))
lines(y=c(6.5,6.5), x=c(9.5, 6.5))
lines(x=c(6.5,6.5), y=c(9.5, 6.5))

lines(x=c(9.5,9.5), y=c(9.5, 11.5))
lines(y=c(9.5,9.5), x=c(9.5, 11.5))
lines(y=c(11.5,11.5), x=c(9.5, 11.5))
lines(x=c(11.5,11.5), y=c(9.5, 11.5))

lines(x=c(14.5,14.5), y=c(14.5, 11.5))
lines(y=c(14.5,14.5), x=c(14.5, 11.5))
lines(y=c(11.5,11.5), x=c(14.5, 11.5))
lines(x=c(11.5,11.5), y=c(14.5, 11.5))

lines(x=c(14.5,14.5), y=c(14.5, 16.5))
lines(y=c(14.5,14.5), x=c(14.5, 16.5))
lines(y=c(16.5,16.5), x=c(14.5, 16.5))
lines(x=c(16.5,16.5), y=c(14.5, 16.5))

dev.off()

# write out a list of related individuals to exclude from future analyses
write.table(related.unique, file='~/Desktop/UT_ddRADseq/keep_files/exclude_related.txt', 
            quote=FALSE, row.names=FALSE, col.names=FALSE)

#LD pruned
ibd.noLD<-snpgdsIBDMLE(genofile, snp.id=snpset.id, autosome.only=FALSE, kinship=TRUE)
ibd.noLD.coeff<-snpgdsIBDSelection(ibd.noLD)
plot(ibd.noLD.coeff$k0, ibd.noLD.coeff$k1, xlim=c(0,1), ylim=c(0,1))
plot(ibd.noLD$kinship, pch=16, col=alpha('blue', 0.5), xlim=c(0,1), ylim=c(0,1))

