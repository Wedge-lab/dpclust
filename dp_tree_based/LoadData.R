load.data <- function(datpath, samplename, list_of_data_files, Chromosome, WT.count, mut.count, subclonal.CN, no.chrs.bearing.mut, data_file_suffix) {
  # Takes:
  # - An array of data file names to be read in
  # - WT.count, the colname of the column in the data files that contains the number of WT reads
  # - mut.count, the name of the column in the data files that contains the number of mutation bearing reads
  # - subclonal.CN, the name of the total copynumber column
  # - no.chrs.bearing.mut, the column that contains the number of chromosomes not bearing any mutations
  data=list()
  for(s in 1:length(list_of_data_files)){
    data[[s]] = read.table(paste(datpath,samplename,list_of_data_files[s],data_file_suffix,sep=""),header=T,sep="\t")
    data[[s]] = data[[s]][data[[s]][,Chromosome] %in% 1:22,]
  }
  
  WTCount = array(0,c(nrow(data[[1]]),no.subsamples))
  mutCount = array(0,c(nrow(data[[1]]),no.subsamples))
  totalCopyNumber = array(0,dim(WTCount))
  copyNumberAdjustment = array(0,dim(WTCount))
  non.deleted.muts = vector(mode="logical",length=nrow(data[[1]]))
  for(s in 1:length(subsamples[[run]])){
    WTCount[,s] = as.numeric(data[[s]][,WT.count])
    mutCount[,s] = as.numeric(data[[s]][,mut.count])
    totalCopyNumber[,s] = as.numeric(data[[s]][,subclonal.CN])
    copyNumberAdjustment[,s] = as.numeric(data[[s]][,no.chrs.bearing.mut])
    non.deleted.muts[data[[s]]$no.chrs.bearing.mut>0]=T
  }
  
  kappa = array(1,dim(mutCount))
  for(i in 1:length(subsamples[[run]])){
    #multiply by no.chrs.bearing.mut, so that kappa is the fraction of reads required for fully clonal mutations, rather than mutation at MCN = 1
    kappa[,i] = mutationCopyNumberToMutationBurden(1,data[[i]]$subclonal.CN,cellularity[i]) * data[[i]]$no.chrs.bearing.mut
  }
  
  # Remove those mutations that have missing values
  no.muts = nrow(totalCopyNumber)
  not.there.wt = sapply(1:no.muts, FUN=function(i, dat1) { sum(is.na(dat1[i,]))>0}, dat1=WTCount)
  not.there.mut = sapply(1:no.muts, FUN=function(i, dat1) { sum(is.na(dat1[i,]))>0}, dat1=mutCount)
  not.there.cn = sapply(1:no.muts, FUN=function(i, dat1) { sum(is.na(dat1[i,]))>0}, dat1=totalCopyNumber)
  not.there.cna = sapply(1:no.muts, FUN=function(i, dat1) { sum(is.na(dat1[i,]))>0}, dat1=copyNumberAdjustment)
  not.there.kappa = sapply(1:no.muts, FUN=function(i, dat1) { sum(is.na(dat1[i,]))>0}, dat1=kappa)
  select = !(not.there.wt | not.there.mut | not.there.cn | not.there.cna | not.there.kappa)
  WTCount = WTCount[select,]
  mutCount = mutCount[select,]
  totalCopyNumber = totalCopyNumber[select,]
  copyNumberAdjustment = copyNumberAdjustment[select,]
  non.deleted.muts = non.deleted.muts[select]
  kappa = kappa[select,]
  print(paste("Removed",no.muts-nrow(WTCount), "mutations with missing data"))
  
  return(list(WTCount=WTCount, mutCount=mutCount, totalCopyNumber=totalCopyNumber, copyNumberAdjustment=copyNumberAdjustment, non.deleted.muts=non.deleted.muts, kappa=kappa))
}