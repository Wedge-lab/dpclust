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
  make_option(c("-k", "--keep_temp_files"), type="logical", action="store_true", default=FALSE, help="Keep the md_out working directory (per-triplet intermediate files) after copying the final result files to outputdir", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

samplename = opt$samplename
output_folder = opt$outputdir
datpath = opt$data_path
purity_file = opt$input
keep_md_out = opt$keep_temp_files

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

#####################################################################################
# Print status message
#####################################################################################
print("")
print(paste("Running:", samplename, sep=" "))
print(paste("Output dir:", output_folder, sep=" "))
print(paste("Analysis type:", "triplet - combine", sep=" "))
print("")

#####################################################################################
# Run combining
#####################################################################################
sample_params = make_sample_params(datafiles=datafiles, cellularity=cellularity, is.male=T, samplename=samplename, subsamples=subsamples)
triplet_params = make_triplet_params(keep_md_out=keep_md_out)

RunTripletDP(step="combine", sample_params=sample_params, triplet_params=triplet_params, outdir=output_folder)
