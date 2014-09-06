source("interconvertMutationBurdens.R")

#
# Testcases
# - getSimpleSubclones with 0.5, 0.3 and 0.1 should yield 3 clusters with allele frequencies around those numbers
# 
#
#
#
#

# Datasets generated:
# cellularity = c(1,1,1)
# no.muts = 50
# sim = simulate.muts(100, 50, 3, mut.frac.of.cells=c(1,1,1))
# ds = sim.muts2dataset(sim, 100, cellularity)
# writeDataset("", "simulated_1", c("1", "2", "3"), "dp_input.txt", ds, cellularity)
# 
# ds2 = sim.muts2dataset(simulate.muts(25, 50, 3, mut.frac.of.cells=c(1,1,0)), 25, cellularity)
# ds3 = sim.muts2dataset(simulate.muts(25, 50, 3, mut.frac.of.cells=c(1,0,1)), 25, cellularity)
# combined = appendSubcloneData(list(ds, ds2, ds3))
# writeDataset("", "simulated_2", c("1", "2", "3"), "dp_input.txt", combined, cellularity)


simulate.muts = function(no.muts, depth.per.cn, no.subsamples, mut.cn=rep(1, no.subsamples), mut.subcl.cn=rep(0, no.subsamples), wt.cn=rep(1, no.subsamples), wt.subcl.cn=rep(0, no.subsamples), mut.frac.of.cells=rep(1, no.subsamples), cellularity=rep(1, no.subsamples), mut.align.bias=1) {
  #
  # Function that simulates mutations. It calculates the exact number of reads we expect to see carying the mutation as well as
  # the depth. It will then put these numbers into a poisson (depth) and binomial (reads for mutation) to simulate readcounts.
  # returned are mutCount, depth, the true.mutCount and the true.WTcount.
  #
  mutCount = array(NA, no.muts)
  WTcount = array(NA, no.muts)
  nmut = matrix(NA, nrow=no.muts, ncol=no.subsamples)
  ndepth = matrix(NA, nrow=no.muts, ncol=no.subsamples)
  for (i in 1:no.subsamples) { # for each subsample
    mutCount.i = cellularity[i]*(mut.frac.of.cells[i]*((mut.subcl.cn[i]+mut.cn[i])*(depth.per.cn*mut.align.bias)))
    WTcount.1.i = cellularity[i]*((1-mut.frac.of.cells[i])*((mut.subcl.cn[i]+mut.cn[i])*(depth.per.cn*(2-mut.align.bias))))
    WTcount.2.i = (1+1-cellularity[i])*(1*((wt.subcl.cn[i]+wt.cn[i])*(depth.per.cn*(2-mut.align.bias))))
    WTcount.i = WTcount.1.i + WTcount.2.i
    
    # Depth takes on a poisson distribution
    ndepth.i = rpois(no.muts,mutCount.i+WTcount.i)
    # mutCount a binomial
    nmut.i = rbinom(no.muts, round(mutCount.i+WTcount.i), mutCount.i/(mutCount.i+WTcount.i))
    print(nmut.i)    
    # Store this subsample
    mutCount[i] = mutCount.i
    WTcount[i] = WTcount.i
    nmut[,i] = nmut.i
    ndepth[,i] = ndepth.i
  }
  return(list(mutCount=nmut, depth=ndepth, true.mutCount=mutCount, true.WTcount=WTcount))#, WTcount.1=WTcount.1, WTcount.2=WTcount.2))
}
sim.muts2dataset = function(sim, no.muts, cellularity) {
  #
  # Converts simulated mutations into a dataset.
  #
  no.subsamples = ncol(sim$mutCount)
  dataset = list()
  dataset$mutCount = sim$mutCount
  dataset$WTCount = sim$depth-sim$mutCount
  dataset$copyNumberAdjustment = matrix(rep(1, no.subsamples*no.muts), ncol=no.subsamples)
  dataset$totalCopyNumber = matrix(rep(2, no.subsamples*no.muts), ncol=no.subsamples)
  dataset$mutation.copy.number = sim$mutCount/sim$depth * dataset$totalCopyNumber
  dataset$kappa = mutationCopyNumberToMutationBurden(1,dataset$totalCopyNumber,cellularity) * dataset$copyNumberAdjustment
  dataset$subclonal.fraction = sim$mutCount / (sim$depth*dataset$kappa)
  dataset$subclonal.fraction[dataset$kappa == 0] = NA
  
  return(dataset)
}

appendSubcloneData <- function(list_of_subclone_datasets) {
  #
  # rbinds a list of datasets together into one
  #
  combined = do.call(cbind, list_of_subclone_datasets)
  mutCount = do.call(rbind, combined['mutCount',])
  WTCount = do.call(rbind, combined['WTCount',])
  copyNumberAdjustment = do.call(rbind, combined['copyNumberAdjustment',])
  totalCopyNumber = do.call(rbind, combined['totalCopyNumber',])
  mutation.copy.number = do.call(rbind, combined['mutation.copy.number',])
  kappa = do.call(rbind, combined['kappa',])
  subclonal.fraction = do.call(rbind, combined['subclonal.fraction',])
  return(list(mutCount=mutCount, WTCount=WTCount, kappa=kappa, copyNumberAdjustment=copyNumberAdjustment, totalCopyNumber=totalCopyNumber, mutation.copy.number=mutation.copy.number, subclonal.fraction=subclonal.fraction))
}

mergeColumns <- function(list_of_datasets) {
  #
  # cbinds a list of datasets together into one
  #
  combined = do.call(cbind, list_of_datasets)
  mutCount = do.call(cbind, combined['mutCount',])
  WTCount = do.call(cbind, combined['WTCount',])
  copyNumberAdjustment = do.call(cbind, combined['copyNumberAdjustment',])
  totalCopyNumber = do.call(cbind, combined['totalCopyNumber',])
  mutation.copy.number = do.call(cbind, combined['mutation.copy.number',])
  kappa = do.call(cbind, combined['kappa',])
  subclonal.fraction = do.call(cbind, combined['subclonal.fraction'])
  return(list(mutCount=mutCount, WTCount=WTCount, kappa=kappa, copyNumberAdjustment=copyNumberAdjustment, totalCopyNumber=totalCopyNumber, mutation.copy.number=mutation.copy.number, subclonal.fraction=subclonal.fraction))
}

simulatedDataToInputFile <- function(dataset, outfile_names) {
  no.muts = nrow(dataset$mutCount)
  no.subsamples = ncol(dataset$mutCount)
  column_names = c("chr","pos","WT.count", "mut.count", "subclonal.CN","mutation.copy.number","subclonal.fraction","no.chrs.bearing.mut")
  
  for (i in 1:no.subsamples) {
    dat = cbind(array(rep(1,no.muts)),
                array(rep(1,no.muts)),
                dataset$WTCount[,i],
                dataset$mutCount[,i],
                dataset$totalCopyNumber[,i],
                dataset$mutation.copy.number[,i],
                dataset$subclonal.fraction[,i],
                dataset$copyNumberAdjustment[,i])

    colnames(dat) = column_names
    write.table(dat,file=outfile_names[i],quote=F,row.names=F, sep="\t")
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
  write.table(dataset_index,file=paste(outpath,samplename,".txt",sep=""),quote=F,row.names=F, sep="\t")
  simulatedDataToInputFile(dataset, outfile_names)
}