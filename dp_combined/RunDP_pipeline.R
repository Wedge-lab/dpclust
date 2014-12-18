outdir = getwd()
args=commandArgs(TRUE)
libdir = toString(args[1]) # Directory where the pipeline is installed
run = as.integer(args[2]) # The sample to be run. Integer that selects from a list of unique samplenames
no.iters = as.integer(args[3]) # Number of iters
no.iters.burn.in = as.double(args[4]) # Number of iters used for burnin
datpath = toString(args[5]) # Where are input files stored
purity_file = toString(args[6]) # A file containing samplenames and purity
analysis_type = toString(args[7]) # String that represents the analysis that is to be run
parallel = as.logical(args[8]) # Supply true or false whether to run parts of the method in parallel
no.of.threads = as.integer(args[9]) # Integer that determines how many threads to use when running parts in parallel

# Optional arguments
if (length(args) >= 10) {
  bin.size = as.double(args[10])
  if (length(args) >= 11) {
    blockid = as.integer(args[11])
    if (length(args) >= 12) {
      no.of.blocks = as.integer(args[12])
    } else {
      no.of.blocks = 1
    }
  } else { 
    blockid = 1
    no.of.blocks = 1
  }
} else {
  bin.size = NA
  blockid = 1
  no.of.blocks = 1
}

# Check whether a supported analysis_type was supplied
supported_commands = c('nd_dp', "tree_dp", 'tree', 'cons', 'replot_1d', 'replot_nd')
if (!(analysis_type %in% supported_commands)) {
  print(paste("Type of analysis", analysis_type, "unknown."))
  print(paste(c("Specify either ", supported_commands)), sep=" ")
  q(save="no", status=1)
}

# Source the required files
setwd(libdir)
source("RunDP.R")
source("LoadData.R")
setwd(outdir)

# Parse the input file and obtain the required data for this run
sample2purity = read.table(purity_file, header=T, stringsAsFactors=F)
samplename = unique(sample2purity$sample)[run]
datafiles = sample2purity[sample2purity$sample==samplename,]$datafile
subsamples = sample2purity[sample2purity$sample==samplename,]$subsample
cellularity = sample2purity[sample2purity$sample==samplename,]$cellularity

print("")
print(paste("Running:", samplename, sep=" "))
print(paste("Working dir:", outdir, sep=" "))
print(paste("Analysis type:", analysis_type, sep=" "))
print("Datafiles:")
print(datafiles)
print("")


if (analysis_type == "tree_dp" | analysis_type == 'tree' | analysis_type == 'cons') {
  outdir = paste(outdir, "/", samplename, "_DPoutput_treeBased_", no.iters,"iters_",no.iters.burn.in,"burnin_", sep="")
  if (!is.na(bin.size)) {
    outdir = paste(outdir, "_",bin.size, "binsize")
  }
} else if (analysis_type == 'nd_dp') {
  outdir = paste(outdir, "/", samplename, "_DPoutput_", no.iters,"iters_",no.iters.burn.in,"burnin", sep="")
}

dataset = load.data(datpath,
                    "",
                    datafiles, 
                    cellularity=cellularity, 
                    Chromosome="chr", 
                    position="end",
                    WT.count="WT.count", 
                    mut.count="mut.count", 
                    subclonal.CN="subclonal.CN", 
                    no.chrs.bearing.mut="no.chrs.bearing.mut", 
                    mutation.copy.number="mutation.copy.number", 
                    subclonal.fraction="subclonal.fraction", 
                    data_file_suffix="")

RunDP(analysis_type=analysis_type, 
      dataset=dataset, 
      samplename=samplename, 
      subsamples=subsamples, 
      no.iters=no.iters, 
      no.iters.burn.in=no.iters.burn.in, 
      outdir=outdir, 
      conc_param=0.01, 
      cluster_conc=5, 
      resort.mutations=T, 
      parallel=parallel, 
      blockid=blockid, 
      no.of.blocks=no.of.blocks, 
      remove.node.frequency=12, #12
      remove.branch.frequency=51, #51
      annotation=vector(mode="character",length=nrow(dataset$mutCount)),
      init.alpha=0.01, 
      shrinkage.threshold=0.1,
      bin.size=bin.size)
