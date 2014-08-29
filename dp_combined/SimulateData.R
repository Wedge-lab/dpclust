source("interconvertMutationBurdens.R")

getSimpleSubclones <- function(n, size, prob, cellularity) {
  mutCount = c()
  WTCount = c()
  no.chrs.bearing.mut = c()
  totalCopyNumber = c()
  for (i in 1:length(n)) {
    new.muts = rbinom(n[i], size[i], prob[i])
    mutCount = c(mutCount,new.muts)
    WTCount = c(WTCount, rep(size[i],n[i])-new.muts)
    no.chrs.bearing.mut = c(no.chrs.bearing.mut, rep(1,n[i]))
    totalCopyNumber = c(totalCopyNumber, rep(2,n[i]))
  }
  mutation.copy.number = mutCount/(mutCount+WTCount) * totalCopyNumber
  mutation.copy.number[is.na(mutation.copy.number)] = 0
  
  kappa = mutationCopyNumberToMutationBurden(1,totalCopyNumber,cellularity) * no.chrs.bearing.mut
  
  return(list(mutCount=matrix(mutCount), WTCount=matrix(WTCount), kappa=matrix(kappa), copyNumberAdjustment=matrix(no.chrs.bearing.mut), totalCopyNumber=matrix(totalCopyNumber), mutation.copy.number=matrix(mutation.copy.number)))
}

createCNA <- function(n, lambda=rep(1, length(n)), round_numbers=T) {
  copyNumber = c()
  for (i in 1:length(n)) {
    copyNumber = c(copyNumber, rpois(n[i],lambda))
  }

  if (!round_numbers) {
    # Add some noise only to those where there is a CNA
    # TODO: perhaps it would be better to use the same random number for stretches of mutations
    inv = copyNumber>0
    copyNumber[inv] = copyNumber[inv] + runif(sum(inv),0,1)
  }
  # TODO: possibly add large (>10) number of copies with a very small chance
  
  return(copyNumber)
}

applyCNA <- function(mutCount, WTCount, no.chrs.bearing.mut, totalCopyNumber, cellularity, new_copynumber) {
  for (i in 1:length(new_copynumber)) {
    cna_on_mutation = runif(1) < 0.5
    if (cna_on_mutation) {
      mutCount[i] = mutCount[i] * new_copynumber[i]
      if (new_copynumber[i] > 0) {
        no.chrs.bearing.mut[i] = no.chrs.bearing.mut[i] + new_copynumber[i]
      } else {
        no.chrs.bearing.mut[i] = 0
      }
    } else {
      WTCount[i] = WTCount[i] * new_copynumber[i]
    }
    totalCopyNumber[i] = totalCopyNumber[i] + new_copynumber[i]
  }
  mutation.copy.number = mutCount/(mutCount+WTCount) * totalCopyNumber
  mutation.copy.number[is.na(mutation.copy.number)] = 0
  
  kappa = mutationCopyNumberToMutationBurden(1,totalCopyNumber,cellularity) * no.chrs.bearing.mut
  return(list(mutCount=matrix(mutCount), WTCount=matrix(WTCount), kappa=matrix(kappa), copyNumberAdjustment=matrix(no.chrs.bearing.mut), totalCopyNumber=matrix(totalCopyNumber), mutation.copy.number=matrix(mutation.copy.number)))
}

# createCopyNumberAdjustment <- function(mutCount) {
# #   alleleAffected = runif(n,0,1) <= 0.5
#   # is 1, unless there is no mutCount  
#   copyNumberAdjustment = rep(1,length(mutCount))
#   copyNumberAdjustment[mutCount == 0] = 0
#   return(copyNumberAdjustment)
# }
# 
# createDepthProfile <- function() {
#   function that creates varying depth
# }

appendSubcloneData <- function(list_of_subclone_datasets) {
  combined = do.call(cbind, list_of_subclone_datasets)
  mutCount = do.call(rbind, combined['mutCount',])
  WTCount = do.call(rbind, combined['WTCount',])
  copyNumberAdjustment = do.call(rbind, combined['copyNumberAdjustment',])
  totalCopyNumber = do.call(rbind, combined['totalCopyNumber',])
  mutation.copy.number = do.call(rbind, combined['mutation.copy.number',])
  kappa = do.call(rbind, combined['kappa',])
  return(list(mutCount=mutCount, WTCount=WTCount, kappa=kappa, copyNumberAdjustment=copyNumberAdjustment, totalCopyNumber=totalCopyNumber, mutation.copy.number=mutation.copy.number))
}

mergeData <- function(list_of_datasets) {
  combined = do.call(cbind, list_of_datasets)
  mutCount = do.call(cbind, combined['mutCount',])
  WTCount = do.call(cbind, combined['WTCount',])
  copyNumberAdjustment = do.call(cbind, combined['copyNumberAdjustment',])
  totalCopyNumber = do.call(cbind, combined['totalCopyNumber',])
  mutation.copy.number = do.call(cbind, combined['mutation.copy.number',])
  kappa = do.call(cbind, combined['kappa',])
  return(list(mutCount=mutCount, WTCount=WTCount, kappa=kappa, copyNumberAdjustment=copyNumberAdjustment, totalCopyNumber=totalCopyNumber, mutation.copy.number=mutation.copy.number))
}

simulatedDataToInputFile <- function(dataset, outfile_names) {
  no.muts = nrow(dataset$mutCount)
  no.subsamples = ncol(dataset$mutCount)
  column_names = c("Chromosome","WT.count", "mut.count", "subclonal.CN","mutation.copy.number","no.chrs.bearing.mut")
  
  for (i in 1:no.subsamples) {
    dat = cbind(rep(1,no.muts),dataset$WTCount[,i],dataset$mutCount[,i],dataset$totalCopyNumber[,i], dataset$mutation.copy.number[,i],dataset$copyNumberAdjustment[,i])
    colnames(dat) = column_names
    write.table(dat,filename=outfile_names[i],quote=F,row.names=F, sep="\t")
  }
}

writeDataset <- function(outpath, samplename, subsamplenames, datafile_suffix, dataset, cellularity) {
#   sample  subsample       datafile        cellularity
#   simulated_2     sample1 simulated_2_sample1_3clusters_150muts_no_cna.txt        1
#   simulated_2     sample2 simulated_2_sample2_3clusters_150muts_no_cna.txt        1
#   simulated_2     sample3 simulated_2_sample3_3clusters_150muts_no_cna.txt        1
  no.subsamples = length(subsamplenames)
  
  outfile_names = c()
  for (i in 1:no.subsamples) {
    outfile_names = c(outfile_names, paste(outpath,samplename,"_",subsamplenames[i],"_",datafile_suffix, sep=""))
  }
  
  dataset_index = data.frame(sample=rep(samplename, no.subsamples), subsample=subsamplenames, datafile=outfile_names, cellularity=cellularity)
  write.table(dataset_index,filename=paste(outpath,samplename,".txt",sep=""),quote=F,row.names=F, sep="\t")
  simulatedDataToInputFile(dataset, outfile_names)
}