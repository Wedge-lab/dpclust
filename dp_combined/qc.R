library(ggplot2)
library(reshape2)
library(gtools)

args=commandArgs(TRUE)
infile = toString(args[1])
datpath = toString(args[2])
outpath = toString(args[3])

# wd = getwd()
# setwd("~/repo/dirichlet/dp_combined")
# source('interconvertMutationBurdens.R')
# source('LoadData.R')
# setwd(wd)
library(DPClust)

# infile = "_allDirichletProcessInfo.txt"
# datpath = "/lustre/scratch110/sanger/sd11/epitax/dirichlet_input/"
# outpath = "~/epi/"

createPng <- function(p, filename, width, height) {
  png(filename=filename, width=width, height=height)
  print(p)
  dev.off()
}

meltFacetPlotData <- function(data, subsamplenames) {
  d = as.data.frame(data)
  colnames(d) = subsamplenames
  data.m = melt(d)
  return(data.m)
}

createHistFacetPlot <- function(data, title, xlab, ylab, binwidth) {
  p = ggplot(data) + aes(x=value, y=..count..) + geom_histogram(binwidth=binwidth, colour="black", fill="gray") + facet_grid(variable ~ .)
  p = p + theme_bw(base_size=35) + ggtitle(title) + xlab(xlab) + ylab(ylab)
  return(p)
}

createBoxFacetPlot <- function(data, title, xlab, ylab) {
  p = ggplot(data) + aes(x=variable, y=value) + geom_boxplot() + facet_grid(subsample ~ .)
  p = p + theme_bw() + ggtitle(title) + xlab(xlab) + ylab(ylab)
  return(p)
}

createViolinFacetPlot <- function(data, title, xlab, ylab) {
  p = ggplot(data) + aes(x=variable, y=value) + geom_violin() + facet_grid(subsample ~ .)
  p = p + theme_bw() + ggtitle(title) + xlab(xlab) + ylab(ylab)
  return(p)
}

createQCDocument <- function(res, samplename, subsamplenames, outpath, cellularity) {
  p = createHistFacetPlot(meltFacetPlotData(res$totalCopyNumber, subsamplenames), paste(samplename, "totalCopyNumber"), "Copynumber", "Count", binwidth=1)
  createPng(p, paste(outpath, samplename, "_totalCopyNumber.png", sep=""), width=1500, height=500*length(subsamplenames))
  
  p = createHistFacetPlot(meltFacetPlotData(res$mutation.copy.number, subsamplenames), paste(samplename, "mutation.copy.number"), "mutation.copy.number", "Count", binwidth=0.1)
  p = p + xlim(0,5)
  createPng(p, paste(outpath, samplename, "_mutation.copy.number.png", sep=""), width=1500, height=500*length(subsamplenames))

  p = createHistFacetPlot(meltFacetPlotData(res$copyNumberAdjustment, subsamplenames), paste(samplename, "copyNumberAdjustment"), "copyNumberAdjustment", "Count (log10)", binwidth=1)
  p = p + scale_y_log10()
  createPng(p, paste(outpath, samplename, "_copyNumberAdjustment.png", sep=""), width=1500, height=500*length(subsamplenames))
  
  p = createHistFacetPlot(meltFacetPlotData(res$mutCount/(res$mutCount+res$WTCount), subsamplenames), paste(samplename, "alleleFrequency"), "Allele Frequency", "Count", binwidth=0.01)
  p = p + xlim(0,1)
  createPng(p, paste(outpath, samplename, "_alleleFrequency.png", sep=""), width=1500, height=500*length(subsamplenames))
  
  p = createHistFacetPlot(meltFacetPlotData(res$kappa, subsamplenames), paste(samplename, "kappa"), "Kappa", "Count", binwidth=0.01)
  createPng(p, paste(outpath, samplename, "_kappa.png", sep=""), width=1500, height=500*length(subsamplenames))
  
#   p = createHistFacetPlot(meltFacetPlotData(res$subclonal.fraction, subsamplenames), paste(samplename, "subclonal.fraction"), "Subclonal fraction", "Count", binwidth=0.01)
#   createPng(p, paste(outpath, samplename, "_subclonal.fraction.png", sep=""), width=1500, height=500*length(subsamplenames))

  p = createHistFacetPlot(meltFacetPlotData((res$mutCount+res$WTCount), subsamplenames), paste(samplename, "depth"), "Depth", "Count", binwidth=5)
  createPng(p, paste(outpath, samplename, "_depth.png", sep=""), width=1500, height=500*length(subsamplenames))

#   fractionOfCells = res$mutation.copy.number / res$copyNumberAdjustment
#   meltFacetPlotData(fractionOfCells, subsamplenames)
  p = createHistFacetPlot(meltFacetPlotData(res$subclonal.fraction, subsamplenames), paste(samplename, "Fraction Of Cells"), "Fraction of Cells", "Count", binwidth=0.05)
  p = p + geom_vline(xintercept=0.5, colour="red", linetype="longdash", size=2) + geom_vline(xintercept=1.5, colour="red", linetype="longdash", size=2) + xlim(0,3) #scale_y_log10()
  createPng(p, paste(outpath, samplename, "_fractionOfCells.png", sep=""), width=1500, height=500*length(subsamplenames))
  
#   manualMutCopyNum = mutationBurdenToMutationCopyNumber(res$mutCount/(res$mutCount+res$WTCount),res$totalCopyNumber ,cellularity)
#   p = createHistFacetPlot(meltFacetPlotData(manualMutCopyNum, subsamplenames), paste(samplename, "manualMutCopyNum"), "Mutation copy number", "Count", binwidth=0.1)
#   createPng(p, paste(outpath, samplename, "_manualMutCopyNum.png", sep=""), width=500, height=250*length(subsamplenames))
#   
#   manualFractionOfCells = mutationBurdenToMutationCopyNumber(res$mutCount/(res$mutCount+res$WTCount),res$totalCopyNumber ,cellularity) / res$copyNumberAdjustment
#   p = createHistFacetPlot(meltFacetPlotData(manualFractionOfCells, subsamplenames), paste(samplename, "manualFractionOfCells"), "Fraction of Cells", "Count", binwidth=0.01)
#   createPng(p, paste(outpath, samplename, "_manualFractionOfCells.png", sep=""), width=1500, height=500*length(subsamplenames))

  # Manually melt the data
  d.m = data.frame()
  for (i in 1:length(subsamplenames)) {
    d.loc = data.frame(subsample=subsamplenames[i],
                       variable=factor(res$chromosome[,i], levels=gtools::mixedsort(unique(res$chromosome))), 
                       value=res$subclonal.fraction[,i])
    d.m = rbind(d.m, d.loc)
  }
  p = createBoxFacetPlot(d.m, paste(samplename, "subclonal fraction per chrom"), "Chromosome", "Subclonal fraction")
  p = p + ylim(0,3)
  createPng(p, paste(outpath, samplename, "_subclonalFractionPerChromosome.png", sep=""), width=1500, height=500*length(subsamplenames))
  
  # Manually melt the data
  d = as.data.frame(res$chromosome)
  colnames(d) = subsamplenames
  d.m = data.frame()
  for (i in 1:ncol(d)) {
    d.loc = d[res$subclonal.fraction > 1.5,i]
    if(length(d.loc) > 0) { 
      d.m = rbind(d.m, data.frame(variable=colnames(d)[i], 
                                  value=factor(d.loc, levels=gtools::mixedsort(unique(res$chromosome)))))
    }
  }
  
  if (nrow(d.m) > 0) {
    p = createHistFacetPlot(d.m, paste(samplename, "subclonal fraction > 1.5"), "Chromosome", "Count", binwidth=1)
    p = p + ylim(0,3)
    createPng(p, paste(outpath, samplename, "_large.subclonal.fraction.by.chrom.png", sep=""), width=1500, height=500*length(subsamplenames))
  }

  # Plot mutation copy number per chromosome
  d.m = data.frame()
  for (i in 1:length(subsamplenames)) {
    d.loc = data.frame(subsample=subsamplenames[i],
                       variable=factor(res$chromosome[,i], levels=gtools::mixedsort(unique(res$chromosome))), 
                       value=res$mutation.copy.number[,i])
    d.m = rbind(d.m, d.loc)
  }
  p = createViolinFacetPlot(d.m, paste(samplename, "mutation copy number per chrom"), "Chromosome", "Mutation Copy Number")
  p = p + ylim(0,3)
  createPng(p, paste(outpath, samplename, "_MutationCopyNumberPerChromosome.png", sep=""), width=1500, height=500*length(subsamplenames))

  
  # Plot frac reads reporting mutation vs total depth
  d.m = data.frame()
  for (i in 1:ncol(res$mutCount)) {
    d.m = rbind(d.m, data.frame(x=res$mutCount/(res$mutCount+res$WTCount), y=(res$mutCount+dataset$WTCount), samplename=rep(subsamplenames[i], nrow(res$mutCount))))
  }  
    
  if (nrow(d.m) > 0) {
    p = ggplot(d.m) + aes(x=x, y=y) + geom_point() + facet_grid(. ~ samplename)
    p = p + theme_bw(base_size=35) + ggtitle(samplename) + ylab("Total reads") + xlab("Fraction of mutant reads")
    createPng(p, paste(outpath, samplename, "_depth.vs.frac.mutCount.png", sep=""), width=1500, height=500*length(subsamplenames))
  }

  # Create summary of AF space
  z = as.data.frame(res$mutCount/(res$mutCount+res$WTCount))
  if (ncol(z) == 1) {
    af.summary = data.frame("Min."=min(z[,1]), "1st Qu."=quantile(z[,1],0.25), "Median"=quantile(z[,1],0.5), "Mean"=mean(z[,1],0.5), "3rd Qu."=quantile(z[,1],0.75), "Max."=max(z[,1]))
    colnames(af.summary) = c("Min.","1st.Qu.","Median","Mean","3rd.Qu.","Max.")
  } else {
    af.summary = as.data.frame(t(apply(z, 1, summary)))
  }

  af.summary$samplename = rep(samplename, ncol(res$mutCount))
  af.summary$subsample = subsamplenames
  write.table(af.summary, file=paste(outpath, samplename, "_af_summary.txt", sep=""), quote=F)

  # Plot AF divided by the cellularity
  dat = res$mutCount/(res$mutCount+res$WTCount)/cellularity
  p = createHistFacetPlot(meltFacetPlotData(dat, subsamplenames), paste(samplename, "Cellularity Corr. AF"), "Cellularity Corrected Allele Frequency", "Count", binwidth=0.01)
  createPng(p, paste(outpath, samplename, "_cellularityCorrectedAF.png", sep=""), width=1500, height=500*length(subsamplenames))


  
  
#   
#   
#   
#   d = res$mutCount/(res$mutCount+res$WTCount)
#   for (i in 1:length(subsamplenames)) {
#     for (j in 1:length(subsamplenames)) {
#       if (!(i==j) && !(i>j)) {
#         d.loc = as.data.frame(d[,c(i,j)])
#         colnames(d.loc) = c(subsamplenames[i], subsamplenames[j])
#         p = ggplot(d.loc) + aes_string(x=subsamplenames[i], y=subsamplenames[j]) + geom_point()
#         print(p)
#       }
#     }
#   }
#   dev.off()
}


# if (infile == "simulated") {
# 	samplename = "sim_multisample_noCN_2"
# 	subsamples = list("first", "second", "third")
# 	load(file='/lustre/scratch110/sanger/sd11/dirichlet/simulated/sim_multisample_noCN_2.RData')
# 	createQCDocument(combined, samplename, subsamples, outpath, cellularity)
# } else {
#     sample2purity = read.table(infile, header=T, stringsAsFactors=F)
#     for (samplename in unique(sample2purity$sample)) {
#       datafiles = sample2purity[sample2purity$sample==samplename,]$datafile
#       subsamples = sample2purity[sample2purity$sample==samplename,]$subsample
#       cellularity = sample2purity[sample2purity$sample==samplename,]$cellularity
#       dataset = load.data(datpath,"",datafiles, cellularity=cellularity, Chromosome="chr", WT.count="WTCount", mut.count="mutCount", subclonal.CN="subclonal.CN", no.chrs.bearing.mut="no.chrs.bearing.mut", mutation.copy.number="mutation.copy.number", data_file_suffix="")
#       createQCDocument(dataset, samplename, subsamples, outpath, cellularity)
#     }
# }

# c("PD7404a_allDirichletProcessInfo.txt","PD7405a_allDirichletProcessInfo.txt","PD7407a_allDirichletProcessInfo.txt")

# dat = load.data(datpath, "", list(infile), c(0.35224), "chr", "WT.count", "mut.count", "subclonal.CN", "no.chrs.bearing.mut", "mutation.copy.number", "subclonal.fraction", "")
# createQCDocument(dat, "PD7404", c("a"), outpath, c(0.35224))

sample2purity = read.table(infile, header=T, stringsAsFactors=F)
for (samplename in unique(sample2purity$sample)) {
  print(samplename)
  datafiles = sample2purity[sample2purity$sample==samplename,]$datafile
  subsamples = sample2purity[sample2purity$sample==samplename,]$subsample
  cellularity = sample2purity[sample2purity$sample==samplename,]$cellularity
  if (file.exists(paste(datpath,datafiles, sep=""))) {
    dataset = load.data(c(paste(datpath, datafiles, sep="")), cellularity=cellularity, Chromosome="chr", position="end", WT.count="WT.count", mut.count="mut.count", subclonal.CN="subclonal.CN", no.chrs.bearing.mut="no.chrs.bearing.mut", mutation.copy.number="mutation.copy.number", subclonal.fraction="subclonal.fraction")
    createQCDocument(dataset, samplename, subsamples, outpath, cellularity)
  }
}
