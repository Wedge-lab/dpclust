suppressMessages(library(DPClust))
library(optparse)

#####################################################################################
# Command line options
#####################################################################################
option_list = list(
  make_option(c("-s", "--samplename"), type="character", default=NULL, help="Name of sample to run", metavar="character"),
  make_option(c("-d", "--data_path"), type="character", default=NULL, help="Path to where dpinput data files are stored", metavar="character"),
  make_option(c("-o", "--outputdir"), type="character", default=NULL, help="Directory where the output is saved", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL, help="Datafile with sample information", metavar="character"),
  make_option(c("--choose.from"), type="integer", default=NULL, help="The total number of subsamples to choose triplets from", metavar="character"),
  make_option(c("--choose.number"), type="integer", default=NULL, help="The number to be picked for each triplet (always 3)"),
  make_option(c("--choose.index"), type="integer", default=NULL, help="Of all possible combinations of subsamples, run this combination", metavar="character"),
  make_option(c("--no.iters"), type="integer", default=1000, help="Number of MCMC iterations. Must match the value used in step 1", metavar="character"),
  make_option(c("--burn.in"), type="integer", default=NULL, help="Number of iterations discarded as burn in. Defaults to ceiling(no.iters/5) when not given. Must match the value used in step 1", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

samplename = opt$samplename
output_folder = opt$outputdir
datpath = opt$data_path
purity_file = opt$input
choose.from = opt$choose.from
choose.number = opt$choose.number
choose.index = opt$choose.index
no.iters = opt$no.iters
no.iters.burn.in = ifelse(is.null(opt$burn.in), ceiling(no.iters/5), opt$burn.in)

#####################################################################################
# Process input
#####################################################################################
sample2purity = read.table(purity_file, header=T, stringsAsFactors=F,
                           colClasses=c(sample="character", subsample="character",
                                        datafile="character", cellularity="numeric"))

if (!samplename %in% sample2purity$sample) {
  stop("Given samplename not found in sample column of samplesheet")
}

datafiles = sample2purity[sample2purity$sample==samplename,]$datafile
# we assume the samplesheet contains just the filename, here we prepend the path to it
datafiles = file.path(datpath, datafiles)
subsamples = sample2purity[sample2purity$sample==samplename,]$subsample
subsamples = paste0(samplename, "_", subsamples)
cellularity = sample2purity[sample2purity$sample==samplename,]$cellularity

subsample.indices = get.subsample.indices(choose.index, choose.number, choose.from)

#####################################################################################
# Print status message
#####################################################################################
print("")
print(paste("Running:", samplename, sep=" "))
print(paste("Output dir:", output_folder, sep=" "))
print(paste("Analysis type:", "triplet - density", sep=" "))
print(paste("Triplet:", paste(subsample.indices, collapse=","), sep=" "))
print("")

#####################################################################################
# Run density estimation
#####################################################################################
sample_params = make_sample_params(datafiles=datafiles, cellularity=cellularity, is.male=T, samplename=samplename, subsamples=subsamples)
triplet_params = make_triplet_params(no.iters=no.iters, no.iters.burn.in=no.iters.burn.in)

RunTripletDP(step="density", sample_params=sample_params, triplet_params=triplet_params, outdir=output_folder, subsample.indices=subsample.indices)
