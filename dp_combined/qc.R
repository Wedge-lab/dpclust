library(ggplot2)
library(reshape2)

args=commandArgs(TRUE)
infile = toString(args[1])
datpath = toString(args[2])
outpath = toString(args[3])

wd = getwd()
setwd("~/repo/dirichlet/dp_combined")
source('interconvertMutationBurdens.R')
source('LoadData.R')
setwd(wd)

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
  p = ggplot(data) + aes(x=value, y=..count..) + geom_histogram(binwidth=binwidth) + facet_grid(variable ~ .)
  p = p + theme_bw(base_size=35) + ggtitle(title) + xlab(xlab) + ylab(ylab)
  return(p)
}

createBoxFacetPlot <- function(data, title, xlab, ylab) {
  p = ggplot(data) + aes(x=variable, y=value) + geom_boxplot() + facet_grid(subsample ~ .)
  p = p + theme_bw() + ggtitle(title) + xlab(xlab) + ylab(ylab)
  return(p)
}

createQCDocument <- function(res, samplename, subsamplenames, outpath, cellularity) {
  p = createHistFacetPlot(meltFacetPlotData(res$totalCopyNumber, subsamplenames), paste(samplename, "totalCopyNumber"), "Copynumber", "Count", binwidth=1)
  createPng(p, paste(outpath, samplename, "_totalCopyNumber.png", sep=""), width=1500, height=500*length(subsamplenames))
  
  p = createHistFacetPlot(meltFacetPlotData(res$mutation.copy.number, subsamplenames), paste(samplename, "mutation.copy.number"), "mutation.copy.number", "Count", binwidth=0.1)
  createPng(p, paste(outpath, samplename, "_mutation.copy.number.png", sep=""), width=1500, height=500*length(subsamplenames))

  p = createHistFacetPlot(meltFacetPlotData(res$copyNumberAdjustment, subsamplenames), paste(samplename, "copyNumberAdjustment"), "copyNumberAdjustment", "Count (log10)", binwidth=1)
  p = p + scale_y_log10()
  createPng(p, paste(outpath, samplename, "_copyNumberAdjustment.png", sep=""), width=1500, height=500*length(subsamplenames))
  
  p = createHistFacetPlot(meltFacetPlotData(res$mutCount/(res$mutCount+res$WTCount), subsamplenames), paste(samplename, "alleleFrequency"), "Allele Frequency", "Count", binwidth=0.01)
  createPng(p, paste(outpath, samplename, "_alleleFrequency.png", sep=""), width=1500, height=500*length(subsamplenames))
  
  p = createHistFacetPlot(meltFacetPlotData(res$kappa, subsamplenames), paste(samplename, "kappa"), "Kappa", "Count", binwidth=0.01)
  createPng(p, paste(outpath, samplename, "_kappa.png", sep=""), width=1500, height=500*length(subsamplenames))
  
#   p = createHistFacetPlot(meltFacetPlotData(res$subclonal.fraction, subsamplenames), paste(samplename, "subclonal.fraction"), "Subclonal fraction", "Count", binwidth=0.01)
#   createPng(p, paste(outpath, samplename, "_subclonal.fraction.png", sep=""), width=1500, height=500*length(subsamplenames))

  p = createHistFacetPlot(meltFacetPlotData((res$mutCount+res$WTCount), subsamplenames), paste(samplename, "depth"), "Depth", "Count", binwidth=5)
  createPng(p, paste(outpath, samplename, "_depth.png", sep=""), width=1500, height=500*length(subsamplenames))

#   fractionOfCells = res$mutation.copy.number / res$copyNumberAdjustment
#   meltFacetPlotData(fractionOfCells, subsamplenames)
  p = createHistFacetPlot(meltFacetPlotData(res$subclonal.fraction, subsamplenames), paste(samplename, "Fraction Of Cells"), "Fraction of Cells", "Count (log10)", binwidth=0.03)
  p = p + scale_y_log10() + geom_vline(xintercept=0.5, colour="red")+ geom_vline(xintercept=1.5, colour="red")
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
    d.loc = data.frame(subsample=subsamplenames[i],variable=factor(dataset$chromosome[,i]), value=dataset$subclonal.fraction[,i])
    d.m = rbind(d.m, d.loc)
  }
  p = createBoxFacetPlot(d.m, paste(samplename, "subclonal fraction per chrom"), "Chromosome", "Subclonal fraction")
  createPng(p, paste(outpath, samplename, "_subclonalFractionPerChromosome.png", sep=""), width=1500, height=500*length(subsamplenames))
  
  # Manually melt the data
  d = as.data.frame(dataset$chromosome)
  colnames(d) = subsamplenames
  d.m = data.frame()
  for (i in 1:ncol(d)) {
    d.loc = d[dataset$subclonal.fraction > 1.5,i]
    d.m = rbind(d.m, data.frame(variable=colnames(d)[i], value=factor(d.loc)))
  }
  
  p = createHistFacetPlot(d.m, paste(samplename, "subclonal fraction > 1.5"), "Chromosome", "Count", binwidth=1)
  createPng(p, paste(outpath, samplename, "_large.subclonal.fraction.by.chrom.png", sep=""), width=1500, height=500*length(subsamplenames))
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
  dataset = load.data(datpath,"",datafiles, cellularity=cellularity, Chromosome="chr", WT.count="WT.count", mut.count="mut.count", subclonal.CN="subclonal.CN", no.chrs.bearing.mut="no.chrs.bearing.mut", mutation.copy.number="mutation.copy.number", subclonal.fraction="subclonal.fraction", data_file_suffix="")
  createQCDocument(dataset, samplename, subsamples, outpath, cellularity)
}
