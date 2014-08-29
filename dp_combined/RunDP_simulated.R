outdir = getwd()
args=commandArgs(TRUE)
samplename = toString(args[1]) # Sample name
no.iters = as.integer(args[2]) # Number of iters
no.iters.burn.in = as.double(args[3]) # Number of iters used for burnin
datpath = toString(args[4]) # Where are input files stored
purity_file = toString(args[5]) # A file containing samplenames and purity
analysis_type = toString(args[6]) # String that represents the analysis that is to be run
parallel = as.logical(args[7]) # Supply true or false whether to run parts of the method in parallel
no.of.threads = as.integer(args[8]) # Integer that determines how many threads to use when running parts in parallel

# Optional arguments
if (length(args) >= 9) {
  bin.size = as.double(args[9])
  if (length(args) >= 10) {
    blockid = as.integer(args[10])
    if (length(args) >= 11) {
      no.of.blocks = as.integer(args[11])
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

# Source the required files
setwd("~/repo/dirichlet/dp_combined")
source("RunDP.R")
source("LoadData.R")
source("SimulateData.R")
setwd(outdir)

# first = getSimpleSubclones(c(100,25,25), c(100,100,100), c(0.9,0.5,0.3))
# second = getSimpleSubclones(c(100,25,25), c(100,100,0), c(0.9,0.5,0.0))
# third = getSimpleSubclones(c(100,25,25), c(100,0,1000), c(0.9,0.0,0.3))
# combined = mergeData(list(first, second, third))
# cellularity = c(1,1,1)

# first = getSimpleSubclones(c(100,25,25), c(100,100,100), c(0.5,0.2,0.3), 1)
# second = getSimpleSubclones(c(100,25,25), c(100,100,50), c(0.5,0.5,0.0), 1)
# third = getSimpleSubclones(c(100,25,25), c(100,50,100), c(0.5,0.0,0.5), 1)
# combined = mergeData(list(first, second, third))
# cellularity = c(1,1,1)
# save(file='/lustre/scratch110/sanger/sd11/dirichlet/simulated/sim_multisample_noCN_2.RData', combined, cellularity)

samplename = "sim_multisample_noCN_2"
subsamples = list("first", "second", "third")
load(file='/lustre/scratch110/sanger/sd11/dirichlet/simulated/sim_multisample_noCN_2.RData')

RunDP(analysis_type=analysis_type, dataset=combined, samplename=samplename, subsamples=subsamples, no.iters=no.iters, no.iters.burn.in=no.iters.burn.in, outdir=outdir, conc_param=0.01, cluster_conc=5, resort.mutations=T, parallel=T, blockid=blockid, no.of.blocks=no.of.blocks)

