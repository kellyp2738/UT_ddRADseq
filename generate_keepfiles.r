################################################################################################################
# Short script for parsing the tick sample data into different "keep" files that can be used for VCF filtering #
################################################################################################################

## RUN ONLY ONCE ##
# read in the data
sample.data<-read.csv("~/Desktop/UT_ddRADseq/ddRAD_FinalLibrary_SampleInfo.csv")

# reconstruct the combinatoric label used in the VCF file
combo.label<-paste(sample.data$tick.id, '_', sample.data$Inline.Barcode, '-', sample.data$Illumina.Barcode, sep='')
full.data<-cbind(sample.data, combo.label)
write.csv(full.data, file='~/Desktop/UT_ddRADseq/ddRAD_FinalLibrary_SampleInfo_Full.csv', row.names=FALSE)

## RUN AS NEEDED

setwd('~/Desktop/UT_ddRADseq/')
full.data<-read.csv('ddRAD_FinalLibrary_SampleInfo_Full.csv')

# define a function that subsets the data
keep.subset<-function(data, column.idx){
  characteristics<-unique(data[column.idx])
  for(i in characteristics[[1]]){
    print(i)
    keep.name=file.path(getwd(), paste('keep_', i, '.txt', sep=''))
    keep.data=subset(data, data[column.idx]==i)$combo.label
    write.table(keep.data, file=keep.name, row.names=FALSE, col.names=FALSE, quote=FALSE)
  }
}

# separate Harrison Bayou and Starr Ranch Trail
keep.subset(full.data, 6)

# separate raccoon and opossum hosts
keep.subset(full.data, 4)

# separate Harrison Bayou north and south
keep.subset(subset(full.data, coll.site=='HB'), 8)

# separate HB raccoon and opossum hosts
site.sp<-paste(full.data$host.species, full.data$coll.site, sep='_') #make new column to not overwrite other racc and op files
full.data.update<-cbind(full.data, site.sp)
keep.subset(full.data.update, 10)

# separate HB individual hosts
keep.subset(subset(full.data, coll.site=='HB'), 5)
