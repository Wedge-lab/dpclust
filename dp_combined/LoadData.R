load.data <- function(datpath, samplename, list_of_data_files, cellularity, Chromosome, WT.count, mut.count, subclonal.CN, no.chrs.bearing.mut, mutation.copy.number, subclonal.fraction, data_file_suffix) {
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

  no.subsamples = length(list_of_data_files)
  no.muts = nrow(data[[1]])
  
  chromosome = matrix(0,no.muts,no.subsamples)
  WTCount = matrix(0,no.muts,no.subsamples)
  mutCount = matrix(0,no.muts,no.subsamples)
  totalCopyNumber = matrix(0,no.muts,no.subsamples)
  copyNumberAdjustment = matrix(0,no.muts,no.subsamples)
  non.deleted.muts = vector(mode="logical",length=nrow(data[[1]]))
  mutationCopyNumber = matrix(NA,no.muts,no.subsamples)
  subclonalFraction = matrix(NA,no.muts,no.subsamples)
  for(s in 1:length(list_of_data_files)){
    chromosome[,s] = as.numeric(data[[s]][,Chromosome])
    WTCount[,s] = as.numeric(data[[s]][,WT.count])
    mutCount[,s] = as.numeric(data[[s]][,mut.count])
    totalCopyNumber[,s] = as.numeric(data[[s]][,subclonal.CN])
    copyNumberAdjustment[,s] = as.numeric(data[[s]][,no.chrs.bearing.mut])
    non.deleted.muts[data[[s]][,no.chrs.bearing.mut]>0]=T
    mutationCopyNumber[,s] = as.numeric(data[[s]][,mutation.copy.number])
    subclonalFraction[,s] = as.numeric(data[[s]][,subclonal.fraction])
  }
  
  kappa = matrix(1,no.muts,no.subsamples)
  for(i in 1:length(list_of_data_files)){
    #multiply by no.chrs.bearing.mut, so that kappa is the fraction of reads required for fully clonal mutations, rather than mutation at MCN = 1
    kappa[,i] = mutationCopyNumberToMutationBurden(1,data[[i]][,subclonal.CN],cellularity[i]) * data[[i]][,no.chrs.bearing.mut]
  }

  # Remove those mutations that have missing values
  not.there.wt = apply(WTCount, 1, function(x) { sum(is.na(x))>0 })
  not.there.mut = apply(mutCount, 1, function(x) { sum(is.na(x))>0 })
  not.there.cn = apply(totalCopyNumber, 1, function(x) { sum(is.na(x))>0 })
  not.there.cna = apply(copyNumberAdjustment, 1, function(x) { sum(is.na(x))>0 })
  not.there.kappa = apply(kappa, 1, function(x) { sum(is.na(x))>0 })
  # Remove those mutations that have no coverage. These cause for trouble lateron.
  not.coverage = apply(WTCount+mutCount, 1, function(x) { any(x==0) })
  not.cna = apply(copyNumberAdjustment, 1, function(x) { any(x==0) })
  
  print(paste("Removed", sum(not.there.wt),"with missing WTCount", sep=" "))
  print(paste("Removed", sum(not.there.mut),"with missing mutCount", sep=" "))
  print(paste("Removed", sum(not.there.cn),"with missing totalCopyNumber", sep=" "))
  print(paste("Removed", sum(not.there.cna),"with missing copyNumberAdjustment", sep=" "))
  print(paste("Removed", sum(not.there.kappa),"with missing kappa", sep=" "))
  print(paste("Removed", sum(not.coverage),"with no coverage", sep=" "))
  print(paste("Removed", sum(not.cna),"with zero copyNumberAdjustment", sep=" "))

  select = !(not.there.wt | not.there.mut | not.there.cn | not.there.cna | not.there.kappa | not.coverage | not.cna)
  chromosome = as.matrix(chromosome[select,])
  WTCount = as.matrix(WTCount[select,])
  mutCount = as.matrix(mutCount[select,])
  totalCopyNumber = as.matrix(totalCopyNumber[select,])
  copyNumberAdjustment = as.matrix(copyNumberAdjustment[select,])
  non.deleted.muts = as.matrix(non.deleted.muts[select])
  kappa = as.matrix(kappa[select,])
  mutationCopyNumber = as.matrix(mutationCopyNumber[select,])
  subclonalFraction = as.matrix(subclonalFraction[select,])
  print(paste("Removed",no.muts-nrow(WTCount), "mutations with missing data"))

  return(list(chromosome=chromosome, WTCount=WTCount, mutCount=mutCount, totalCopyNumber=totalCopyNumber, copyNumberAdjustment=copyNumberAdjustment, non.deleted.muts=non.deleted.muts, kappa=kappa, mutation.copy.number=mutationCopyNumber, subclonal.fraction=subclonalFraction))
}
