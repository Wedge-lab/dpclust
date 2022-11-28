#begin test
library(DPClust)
source("/gpfs3/users/woodcock/dus183/Data/test/DPClust/AssignMutations.R")
source("/gpfs3/users/woodcock/dus183/Data/test/DPClust/DensityEstimator.R")
source("/gpfs3/users/woodcock/dus183/Data/test/DPClust/DirichletProcessClustering.R")
source("/gpfs3/users/woodcock/dus183/Data/test/DPClust/InformationCriterions.R")
source("/gpfs3/users/woodcock/dus183/Data/test/DPClust/interconvertMutationBurdens.R")
source("/gpfs3/users/woodcock/dus183/Data/test/DPClust/LoadData.R")
source("/gpfs3/users/woodcock/dus183/Data/test/DPClust/pipeline_functions.R")
source("/gpfs3/users/woodcock/dus183/Data/test/DPClust/PlotDensities.R")
source("/gpfs3/users/woodcock/dus183/Data/test/DPClust/SampleMutations.R")
source("/gpfs3/users/woodcock/dus183/Data/test/DPClust/subclone_Dirichlet_Gibbs_sampler_nD_binomial.R")
source("/gpfs3/users/woodcock/dus183/Data/test/DPClust/util.R")

source("/users/woodcock/dus183/Data/test/core_phasing_assign.R")
source("/users/woodcock/dus183/Data/test/RunDP_phasingversion.R")
source("/users/woodcock/dus183/Data/test/function_phasingversion_sample.R")


library(optparse)
option_list = list(
  make_option(c("-r", "--run_sample"), type="integer", default=NULL, help="Sample to run", metavar="character"),
  make_option(c("-d", "--data_path"), type="character", default=NULL, help="Path to where dpinput data files are stored", metavar="character"),
  make_option(c("-o", "--outputdir"), type="character", default=NULL, help="Directory where the output is saved", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL, help="Datafile with sample information", metavar="character"),
  make_option(c("-k", "--keep_temp_files"), type="logical", action="store_true", default=FALSE, help="Keep intermediate output files", metavar="character"),
  make_option(c("-a", "--analysis_type"), type="character", default="nd_dp", help="Analysis type to run", metavar="character"),
  make_option(c("--num_clus_sampler"), type="integer", default=10, help="number of clusters in sampler", metavar="character"),
  make_option(c("--iterations"), type="integer", default=2000, help="Number of iterations to run the MCMC chain", metavar="character"),
  make_option(c("--burnin"), type="integer", default=1000, help="Number of iterations to discard as burnin", metavar="character"),
  make_option(c("--mut_assignment_type"), type="integer", default=1, help="Mutation assignment method", metavar="character"),
  make_option(c("--min_muts_cluster"), type="integer", default=-1, help="Minimum number of mutations per cluster required for it to be kept in the final output, set to -1 to disable (default), see also --min_frac_muts_cluster", metavar="character"),
  make_option(c("--min_frac_muts_cluster"), type="numeric", default=0.01, help="Minimum fraction of mutations per cluster required for it to be kept in the final output, set to -1 to disable, see also --min_muts_cluster", metavar="character"),
  make_option(c("--num_muts_sample"), type="integer", default=50000, help="Number of mutations from which downsampling starts", metavar="character"),
  make_option(c("--bin_size"), type="double", default=NULL, help="Binsize to use when constructing multi-dimensional density - only used when number of samples > 1", metavar="character"),
  make_option(c("--seed"), type="integer", default=123, help="Provide a seed", metavar="character"),
  make_option(c("--assign_sampled_muts"), type="integer", default=TRUE, help="Whether to assign mutations that have been removed during sampling", metavar="character"),
  make_option(c("--error_rate"), type="double", default=0.005, help="Error rate in sequencing for phasing information", metavar="character"),
  make_option(c("--kill_clu_frac"), type="double", default=0.005, help="The threshold for dropping clusters in phasing version", metavar="character")
)



opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

print("test")

run = opt$run_sample
datpath = opt$data_path
outdir = opt$outputdir
purity_file = opt$input
analysis_type = opt$analysis_type #phasing_assign
no.iters = opt$iterations
no.iters.burn.in = opt$burnin
mut.assignment.type = opt$mut_assignment_type
num_muts_sample = opt$num_muts_sample
bin.size = opt$bin_size
seed = opt$seed
assign_sampled_muts = opt$assign_sampled_muts
keep_temp_files = opt$keep_temp_files
min_muts_cluster = opt$min_muts_cluster # set absolute minimum number of mutations per cluster
min_frac_muts_cluster = opt$min_frac_muts_cluster # set proportional minimum number of mutations per cluster

num_clu_samplers = opt$num_clus_sampler #C in gibbs samplers, number of clusters
error.rate = opt$error_rate 
kill.clu.frac =  opt$kill_clu_frac 

if (is.null(outdir)) { outdir = getwd() }

#####################################################################################
# Fixed parameters
#####################################################################################
is.male = T
#not used: min_sampling_factor = 1.1
#not used: is.vcf = F
# assign_sampled_muts = T
sample.snvs.only = T # Perform sampling on just the SNVs and not on CNAs
remove.snvs = F # Clear all SNVs, to perform clustering on CNAs only - This needs a better solution
generate_cluster_ordering = F
species = "human" # mouse also supported, just changes the chromosomes on which mutations are kept, has not effect on functionality

# Cocluster CNA parameters
co_cluster_cna = F
add.conflicts = F # Make the conflicts matrix in a dataset - Flag pertains to both copy number and mut2mut phasing
cna.conflicting.events.only = F # Add only those CNAs that are conflicting
num.clonal.events.to.add = 1 # Add this many clonal CNA events to the clustering
min.cna.size = 100 # Minim size in 10kb for a CNA event to be included

# Remove RGL X11 error
options(rgl.useNULL=TRUE) 
suppressMessages(library(DPClust))

#####################################################################################
# Process input
#####################################################################################
# Parse the input file and obtain the required data for this run
sample2purity = read.table(purity_file, header=T, stringsAsFactors=F)
samplename = unique(sample2purity$sample)[run]
datafiles = sample2purity[sample2purity$sample==samplename,]$datafile
subsamples = sample2purity[sample2purity$sample==samplename,]$subsample
cellularity = sample2purity[sample2purity$sample==samplename,]$cellularity


if ("sex" %in% colnames(sample2purity)) {
  is.male = (sample2purity[sample2purity$sample==samplename,]$sex=="male")[1]
  cndatafiles = sample2purity[sample2purity$sample==samplename,]$cndatafile
} else {
  is.male = T
  cndatafiles = NA
}

if ("mutphasing" %in% colnames(sample2purity)) {
  mutphasingfiles = sample2purity[sample2purity$sample==samplename,]$mutphasing #provide phasing files 
} else {
  mutphasingfiles = NA
}

#####################################################################################
# Print status message
#####################################################################################
print("")
print(paste("Running:", samplename, sep=" "))
print(paste("Working dir:", outdir, sep=" "))
print(paste("Analysis type:", analysis_type, sep=" "))
print("Datafiles:")
print(datafiles)
print("")

#####################################################################################
# Create output path
#####################################################################################
outdir = file.path(outdir, paste(samplename, "_DPoutput_", no.iters,"iters_", no.iters.burn.in, "burnin_seed", seed, "/", sep=""))

#####################################################################################
# Setup parameters
#####################################################################################
run_params = make_run_params(no.iters, no.iters.burn.in, mut.assignment.type, num_muts_sample, is.male=is.male, min_muts_cluster=min_muts_cluster, min_frac_muts_cluster=min_frac_muts_cluster, species=species, assign_sampled_muts=assign_sampled_muts, keep_temp_files=keep_temp_files, generate_cluster_ordering=generate_cluster_ordering)
sample_params = make_sample_params(datafiles, cellularity, is.male, samplename, subsamples, mutphasingfiles)
advanced_params = make_advanced_params(seed)


#####################################################################################
# Run clustering
#####################################################################################
RunDP(analysis_type=analysis_type,
      datpath = datpath,
      run_params=run_params,
      sample_params=sample_params, 
      advanced_params=advanced_params,
      outdir=outdir, 
      cna_params=NULL, 
      mutphasingfiles=mutphasingfiles)
