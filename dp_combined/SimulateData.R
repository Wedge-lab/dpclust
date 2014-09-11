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

# Not used
# sim1 = simulate.muts(5, 50, 3, mut.frac.of.cells=c(1,1,1), mut.cn=c(0,2,3))
# sim2 = simulate.muts(5, 50, 3, mut.frac.of.cells=c(1,1,1), mut.cn=c(3,0,2))
# sim3 = simulate.muts(5, 50, 3, mut.frac.of.cells=c(1,1,1), mut.cn=c(2,3,0))
# sim4 = simulate.muts(5, 50, 3, mut.frac.of.cells=c(1,1,1), wt.cn=c(2,3,0))
# sim5 = simulate.muts(5, 50, 3, mut.frac.of.cells=c(1,1,1), wt.cn=c(0,2,3))
# sim6 = simulate.muts(5, 50, 3, mut.frac.of.cells=c(1,1,1), wt.cn=c(3,0,2))
# sim7 = simulate.muts(5, 50, 3, mut.frac.of.cells=c(1,1,1), wt.subcl.cn=c(0.3,0,0))
# sim8 = simulate.muts(5, 50, 3, mut.frac.of.cells=c(1,1,1), wt.subcl.cn=c(0,0.1,0.1))
# sim9 = simulate.muts(5, 50, 3, mut.frac.of.cells=c(1,1,1), mut.cn=c(0,0.6,0.2))
# sim10 = simulate.muts(5, 50, 3, mut.frac.of.cells=c(1,1,1), mut.cn=c(0.1,0,0.4))
# 
# mut.cn = matrix(c(rpois(15*3,1.5), rep(1, 35*3)), ncol=3, byrow=T)
# mut.cn[,1] = mut.cn[,1]+1
# 
# mut.cn = matrix(c(0,2,3, 3,0,2, 2,3,0, rep(1, ))
# 
# ds2 = simulate.subclone(c(1,1,1), c(50,50,50), 3, mut.cn=mut.cn



# simulate.subclone = function(no.muts, depth.per.cn, no.subsamples, mut.cn=matrix(1, ncol=no.subsamples, nrow=sum(no.muts)), mut.subcl.cn=matrix(1, ncol=no.subsamples, nrow=sum(no.muts)), wt.cn=matrix(1, ncol=no.subsamples, nrow=sum(no.muts)), wt.subcl.cn=matrix(1, ncol=no.subsamples, nrow=sum(no.muts)), mut.frac.of.cells=matrix(1, ncol=no.subsamples, nrow=sum(no.muts)), cellularity=rep(1, no.subsamples), mut.align.bias=1) {
#   simmed.ds = list()
#   for (i in 1:sum(no.muts)) {
#     # simulate a mutation with the given specs and transform it into a dataset
#     simmed.ds[[i]] = sim.muts2dataset(simulate.muts(1, 
#                                                     depth.per.cn[i], 
#                                                     no.subsamples,
#                                                     mut.cn=mut.cn[i,],
#                                                     mut.subcl.cn=mut.subcl.cn[i,],
#                                                     wt.cn=wt.cn[i,],
#                                                     wt.subcl.cn=wt.subcl.cn[i,],
#                                                     mut.frac.of.cells=mut.frac.of.cells[i,],
#                                                     cellularity=cellularity,
#                                                     mut.align.bias=mut.align.bias), 
#                                       1, 
#                                       cellularity)
#   }
#   # Concatenate all the single mutation datasets
#   ds = appendSubcloneData(simmed.ds)
#   return(ds)
# }
create.testing.samples = function() {
  no.subsamples = 3
  #if(0) {
  # Affected by nothing, 3 simple subclones across 3 samples
  ds = generate.dataset(no.subclones=3, 
                        no.muts=c(100,25,25), 
                        no.subsamples=no.subsamples, 
                        cov=c(50,50,50), 
                        mut.frac.of.cells=matrix(c(c(1,0.5,0.3), c(1,0,1), c(1,1,0)), ncol=no.subsamples),
                        sample.cn.param=c(-1,-1,-1), # pois parameter that affects CN values per sample
                        sample.cn.frac=c(0,0,0), # frac of muts affected by CN per sample
                        sample.cn.del.frac=c(0,0,0), # frac of muts that are affected by a deletion
                        sample.subcl.cn.frac=c(0,0,0), # frac of muts affected by subclonal CN per sample
                        cellularity=rep(1, no.subsamples)) #rnorm(3,mean=0.7,sd=0.2))
  writeDataset("", "simulated_001", c("01", "02", "03"), "dp_input.txt", ds, ds$cellularity)

  # Affected by clonal copynumber
  ds = generate.dataset(no.subclones=3, 
                        no.muts=c(100,25,25), 
                        no.subsamples=no.subsamples, 
                        cov=c(50,50,50), 
                        mut.frac.of.cells=matrix(c(c(1,0.5,0.3), c(1,0,1), c(1,1,0)), ncol=no.subsamples),
                        sample.cn.param=c(1.45,1.55,1.5), # pois parameter that affects CN values per sample
                        sample.cn.frac=c(0.25,0.25,0.25), # frac of muts affected by CN per sample
                        sample.cn.del.frac=c(0,0,0), # frac of muts that are affected by a deletion
                        sample.subcl.cn.frac=c(0,0,0), # frac of muts affected by subclonal CN per sample
                        cellularity=rep(1, no.subsamples)) #rnorm(3,mean=0.7,sd=0.2))
  writeDataset("", "simulated_002", c("01", "02", "03"), "dp_input.txt", ds, ds$cellularity)
  
  # Affected by clonal delections
  ds = generate.dataset(no.subclones=3, 
                        no.muts=c(100,25,25), 
                        no.subsamples=no.subsamples, 
                        cov=c(50,50,50), 
                        mut.frac.of.cells=matrix(c(c(1,0.5,0.3), c(1,0,1), c(1,1,0)), ncol=no.subsamples),
                        sample.cn.param=c(-1,-1,-1), # pois parameter that affects CN values per sample
                        sample.cn.frac=c(0,0,0), # frac of muts affected by CN per sample
                        sample.cn.del.frac=c(0.09,0.11,0.1), # frac of muts that are affected by a deletion
                        sample.subcl.cn.frac=c(0,0,0), # frac of muts affected by subclonal CN per sample
                        cellularity=rep(1, no.subsamples)) #rnorm(3,mean=0.7,sd=0.2))
  writeDataset("", "simulated_003", c("01", "02", "03"), "dp_input.txt", ds, ds$cellularity)
  
  # Affected by sublonal copynumber
  ds = generate.dataset(no.subclones=3, 
                        no.muts=c(100,25,25), 
                        no.subsamples=no.subsamples, 
                        cov=c(50,50,50), 
                        mut.frac.of.cells=matrix(c(c(1,0.5,0.3), c(1,0,1), c(1,1,0)), ncol=no.subsamples),
                        sample.cn.param=c(-1,-1,-1), # pois parameter that affects CN values per sample
                        sample.cn.frac=c(0,0,0), # frac of muts affected by CN per sample
                        sample.cn.del.frac=c(0,0,0), # frac of muts that are affected by a deletion
                        sample.subcl.cn.frac=c(0.11,0.09,0.1), # frac of muts affected by subclonal CN per sample
                        cellularity=rep(1, no.subsamples)) #rnorm(3,mean=0.7,sd=0.2))
  writeDataset("", "simulated_004", c("01", "02", "03"), "dp_input.txt", ds, ds$cellularity)
  
  # Affected by lower cellularity
  cellularity = rnorm(no.subsamples,mean=0.7,sd=0.2)
  cellularity[cellularity > 1] = 1
  ds = generate.dataset(no.subclones=3, 
                        no.muts=c(100,25,25), 
                        no.subsamples=no.subsamples, 
                        cov=c(50,50,50), 
                        mut.frac.of.cells=matrix(c(c(1,0.5,0.3), c(1,0,1), c(1,1,0)), ncol=no.subsamples),
                        sample.cn.param=c(-1,-1,-1), # pois parameter that affects CN values per sample
                        sample.cn.frac=c(0,0,0), # frac of muts affected by CN per sample
                        sample.cn.del.frac=c(0,0,0), # frac of muts that are affected by a deletion
                        sample.subcl.cn.frac=c(0,0,0), # frac of muts affected by subclonal CN per sample
                        cellularity=cellularity)
  writeDataset("", "simulated_005", c("01", "02", "03"), "dp_input.txt", ds, ds$cellularity)
  
  # All above combined
  cellularity = rnorm(no.subsamples,mean=0.7,sd=0.2)
  cellularity[cellularity > 1] = 1
  ds = generate.dataset(no.subclones=3, 
                        no.muts=c(100,25,25), 
                        no.subsamples=no.subsamples, 
                        cov=c(50,50,50), 
                        mut.frac.of.cells=matrix(c(c(1,0.5,0.3), c(1,0,1), c(1,1,0)), ncol=no.subsamples),
                        sample.cn.param=c(1.45,1.55,1.5), # pois parameter that affects CN values per sample
                        sample.cn.frac=c(0.25,0.25,0.25), # frac of muts affected by CN per sample
                        sample.cn.del.frac=c(0.09,0.11,0.1), # frac of muts that are affected by a deletion
                        sample.subcl.cn.frac=c(0.11,0.09,0.1), # frac of muts affected by subclonal CN per sample
                        cellularity=cellularity)
  writeDataset("", "simulated_006", c("01", "02", "03"), "dp_input.txt", ds, ds$cellularity)

  # This sample was created to test a few things, but it turned out to be a fundamental bug that affected all samples.
  # They therefore needed to be regenerated anyway. This sample is the same as 006.
  cellularity = rnorm(no.subsamples,mean=0.7,sd=0.2)
  cellularity[cellularity > 1] = 1
  ds = generate.dataset(no.subclones=3, 
                        no.muts=c(100,25,25), 
                        no.subsamples=no.subsamples, 
                        cov=c(50,50,50), 
                        mut.frac.of.cells=matrix(c(c(1,0.5,0.3), c(1,0,1), c(1,1,0)), ncol=no.subsamples),
                        sample.cn.param=c(1.45,1.55,1.5), # pois parameter that affects CN values per sample
                        sample.cn.frac=c(0.25,0.25,0.25), # frac of muts affected by CN per sample
                        sample.cn.del.frac=c(0.09,0.11,0.1), # frac of muts that are affected by a deletion
                        sample.subcl.cn.frac=c(0.11,0.09,0.1), # frac of muts affected by subclonal CN per sample
                        cellularity=cellularity)
  writeDataset("", "simulated_007", c("01", "02", "03"), "dp_input.txt", ds, ds$cellularity)
#   }
}

generate.dataset = function(no.subclones, no.muts, no.subsamples, cov=rep(50, no.subsamples), mut.frac.of.cells=matrix(rep(1, no.subclones*no.subsamples), ncol=no.subsamples), sample.cn.param=rep(-1, no.subsamples), sample.cn.frac=rep(0, no.subsamples), sample.cn.del.frac=rep(0, no.subsamples), sample.subcl.cn.frac=rep(0, no.subsamples), cellularity=rep(1, no.subsamples)) {
  dataset = list()
  for (i in 1:no.subclones) {
#     print(paste("",i))
    subclone = simulate.subclone(no.muts[i], 
                                 no.subsamples=no.subsamples, 
                                 cov=cov, 
                                 mut.frac.of.cells=mut.frac.of.cells[i,],
                                 mut.cn.lambda=sample.cn.param,
                                 wt.cn.lambda=sample.cn.param,
                                 sample.cn.frac=sample.cn.frac, 
                                 sample.cn.del.frac=sample.cn.del.frac, 
                                 sample.subcl.cn.frac=sample.subcl.cn.frac,
                                 cellularity=cellularity)

    dataset[[i]] = mergeColumns(subclone)
    dataset[[i]]$subcloneid = matrix(rep(i, no.muts[i]), ncol=1)
    
#     print(cbind(dataset[[i]]$mutCount, dataset[[i]]$WTCount, dataset[[i]]$subclonal.fraction))
#     print()
  }
  
  dataset = appendSubcloneData(dataset)
  dataset$cellularity = cellularity
  
  return(dataset)
}

simulate.subclone = function(no.muts, no.subsamples, cov, mut.frac.of.cells, mut.cn.lambda=rep(-1, no.subsamples), wt.cn.lambda=rep(-1, no.subsamples), cellularity=rep(1, no.subsamples), sample.cn.frac=rep(0, no.subsamples), sample.cn.del.frac=rep(0, no.subsamples), sample.subcl.cn.frac=rep(0, no.subsamples)) { #.sd=0.2
  
  create.cn = function(no.muts, cn.lambda, cn.frac, del.frac) {
    if (cn.lambda > -1) { 
      # Draw copynumber from pois distribution+1
      cn = rpois(no.muts, cn.lambda)+1
      
      # CN should affect a fixed percentage of these mutations
      sel = rbinom(no.muts, 1, 1-cn.frac)==1
      cn[sel] = 1 # Reset to 1
      
      # A set percentage should be affected by a loss
      sel = rbinom(no.muts, 1, del.frac)==1
      cn[sel] = 0 # Reset to 0
      
    } else {
      # No CN changes
      cn = rep(1, no.muts)
    }
    return(cn)
  }
  
  create.subcl.cn = function(no.muts, subcl.cn.frac) {
    if (subcl.cn.frac > 0) {
      # Draw subclonal copynumber from uniform [0,1]
      subcl.cn = runif(no.muts)
      
      # subclonal CN should affect a fixed percentage of these mutations
      sel = rbinom(no.muts, 1, 1-subcl.cn.frac)==1
      subcl.cn[sel] = 0 # Reset to 0
      
    } else {
      # No subclonal CN
      subcl.cn = rep(0, no.muts)
    }
    return(subcl.cn)
  }
  
  subclone = list()
  for (i in 1:no.subsamples) {
#     print(paste("\t",i))
    mut.cn = create.cn(no.muts, mut.cn.lambda[i], sample.cn.frac[i], sample.cn.del.frac[i])
    mut.subcl.cn = create.subcl.cn(no.muts, sample.subcl.cn.frac[i])
    wt.cn = create.cn(no.muts, wt.cn.lambda[i], sample.cn.frac[i], sample.cn.del.frac[i])
    wt.subcl.cn = create.subcl.cn(no.muts, sample.subcl.cn.frac[i])
    
    muts = generate.mut(no.muts=no.muts, depth.per.cn=cov[i], mut.cn=mut.cn, mut.subcl.cn=mut.subcl.cn, wt.cn=wt.cn, wt.subcl.cn=wt.subcl.cn, mut.frac.of.cells=mut.frac.of.cells[i], cellularity=cellularity[i], mut.align.bias=1)
    subclone[[i]] = sim.muts2dataset(muts, no.muts, cellularity[i])
  }
  return(subclone)
}

# simulate.muts = function(no.muts, depth.per.cn, no.subsamples, mut.cn=rep(1, no.subsamples), mut.subcl.cn=rep(0, no.subsamples), wt.cn=rep(1, no.subsamples), wt.subcl.cn=rep(0, no.subsamples), mut.frac.of.cells=rep(1, no.subsamples), cellularity=rep(1, no.subsamples), mut.align.bias=1) {
#   #
#   # Function that simulates mutations. It calculates the exact number of reads we expect to see carying the mutation as well as
#   # the depth. It will then put these numbers into a poisson (depth) and binomial (reads for mutation) to simulate readcounts.
#   # returned are mutCount, depth, the true.mutCount and the true.WTcount.
#   #
#   mutCount = array(NA, no.muts)
#   WTcount = array(NA, no.muts)
#   nmut = matrix(NA, nrow=no.muts, ncol=no.subsamples)
#   ndepth = matrix(NA, nrow=no.muts, ncol=no.subsamples)
#   for (i in 1:no.subsamples) { # for each subsample
#     
#     mut = generate.mut(no.muts, depth.per.cn, mut.cn[i], mut.subcl.cn[i], wt.cn[i], wt.subcl.cn[i], mut.frac.of.cells[i], cellularity[i], mut.align.bias)
#     
# #     mutCount.i = cellularity[i]*(mut.frac.of.cells[i]*((mut.subcl.cn[i]+mut.cn[i])*(depth.per.cn*mut.align.bias)))
# #     WTcount.1.i = cellularity[i]*((1-mut.frac.of.cells[i])*((mut.subcl.cn[i]+mut.cn[i])*(depth.per.cn*(2-mut.align.bias))))
# #     WTcount.2.i = (1+1-cellularity[i])*(1*((wt.subcl.cn[i]+wt.cn[i])*(depth.per.cn*(2-mut.align.bias))))
# #     WTcount.i = WTcount.1.i + WTcount.2.i
# #     
# #     # Depth takes on a poisson distribution
# #     ndepth.i = rpois(no.muts,mutCount.i+WTcount.i)
# #     # mutCount a binomial
# #     nmut.i = rbinom(no.muts, round(mutCount.i+WTcount.i), mutCount.i/(mutCount.i+WTcount.i))
# #     print(nmut.i)    
#     # Store this subsample
#     mutCount[i] = mut$true.mutCount #.i
#     WTcount[i] = mut$true.WTcount #.i
#     nmut[,i] = mut$nmut #.i
#     ndepth[,i] = mut$ndepth #.i
#   }
#   return(list(mutCount=nmut, depth=ndepth, true.mutCount=mutCount, true.WTcount=WTcount))#, WTcount.1=WTcount.1, WTcount.2=WTcount.2))
# }

generate.mut = function(no.muts, depth.per.cn, mut.cn, mut.subcl.cn, wt.cn, wt.subcl.cn, mut.frac.of.cells, cellularity, mut.align.bias) {
#   mutCount = cellularity*(mut.frac.of.cells*((mut.subcl.cn+mut.cn)*(depth.per.cn*mut.align.bias)))
#   WTcount.1 = cellularity*((1-mut.frac.of.cells)*((mut.subcl.cn+mut.cn)*(depth.per.cn*(2-mut.align.bias))))
#   WTcount.2 = (1+1-cellularity)*(1*((wt.subcl.cn+wt.cn)*(depth.per.cn*(2-mut.align.bias))))
#   print(paste("\t\t",mut.frac.of.cells))
  # Tumour cells carying mut
  mutCount = cellularity*(mut.cn+mut.subcl.cn)*mut.frac.of.cells*depth.per.cn*mut.align.bias
  # Tumour cells carying WT - 1+1-mut.frac.of.cells here accounts for 1*tumour cells WT and (1-mut.frac.of.cells)*tumour cells not mut
  WTcount.1 = cellularity*(wt.cn+wt.subcl.cn)*(1+1-mut.frac.of.cells)*depth.per.cn*(2-mut.align.bias)
  # Stromal cells
  WTcount.2 = rep((1-cellularity)*2*depth.per.cn*(2-mut.align.bias), no.muts) #(1-mut.frac.of.cells)*
  
  WTcount = WTcount.1 + WTcount.2
  WTcount[WTcount==0] = 1 # Setting WTcount to 1 in order to obtain a finite, but very small AF below 
  
  totalCopyNumber = mut.cn+mut.subcl.cn+wt.cn+wt.subcl.cn #rep(mut.cn+mut.subcl.cn+wt.cn+wt.subcl.cn, no.muts)
  copyNumberAdjustment = mut.cn #+mut.subcl.cn #add mut.subcl.cn here to achieve cleaner data. without it some mutations are in too high fraction of cells
  # Depth takes on a poisson distribution
  ndepth = rpois(no.muts,mutCount+WTcount)
  # mutCount a binomial
  #nmut = rbinom(no.muts, round(mutCount+WTcount), mutCount/(mutCount+WTcount))
  nmut = rbinom(no.muts, round(ndepth), mutCount/(mutCount+WTcount))
  nmut[ndepth==0] = 0
  
  return(list(mutCount=nmut, depth=ndepth, true.mutCount=mutCount, true.WTcount=WTcount, totalCopyNumber=totalCopyNumber, copyNumberAdjustment=copyNumberAdjustment))
}


sim.muts2dataset = function(sim, no.muts, cellularity) {
  #
  # Converts simulated mutations into a dataset.
  #
  no.subsamples = ncol(sim$mutCount)
  dataset = list()
  dataset$mutCount = sim$mutCount
  dataset$WTCount = sim$depth-sim$mutCount
  dataset$WTCount[dataset$WTCount < 0] = 0
  dataset$copyNumberAdjustment = sim$copyNumberAdjustment #matrix(rep(1, no.subsamples*no.muts), ncol=no.subsamples)
  dataset$totalCopyNumber = sim$totalCopyNumber #matrix(rep(2, no.subsamples*no.muts), ncol=no.subsamples)
#   dataset$mutation.copy.number = sim$mutCount/sim$depth * dataset$totalCopyNumber
  dataset$mutation.copy.number = mutationBurdenToMutationCopyNumber(sim$mutCount/(sim$depth), sim$totalCopyNumber, cellularity, rep(2,length(sim$mutCount)))
  dataset$kappa = mutationCopyNumberToMutationBurden(1,dataset$totalCopyNumber,cellularity) * dataset$copyNumberAdjustment
  dataset$subclonal.fraction = sim$mutCount / (sim$depth*dataset$kappa)
  dataset$subclonal.fraction[dataset$kappa == 0] = 0
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
  subcloneid = do.call(rbind, combined['subcloneid',])
  return(list(mutCount=mutCount, WTCount=WTCount, kappa=kappa, copyNumberAdjustment=copyNumberAdjustment, totalCopyNumber=totalCopyNumber, mutation.copy.number=mutation.copy.number, subclonal.fraction=subclonal.fraction, subcloneid=subcloneid))
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
  subclonal.fraction = do.call(cbind, combined['subclonal.fraction',])
  return(list(mutCount=mutCount, WTCount=WTCount, kappa=kappa, copyNumberAdjustment=copyNumberAdjustment, totalCopyNumber=totalCopyNumber, mutation.copy.number=mutation.copy.number, subclonal.fraction=subclonal.fraction, subcloneid=list_of_datasets[[1]]$subcloneid))
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
