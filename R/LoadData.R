#' This function loads a series of data tables or VCFs into a dataset object
#' which is in essense a list of tables, one for each type of information
#' 
#' Takes:
#' @param list_of_data_files A list of data file names to be read in
#' @param cellularity A vector containing the cellularity estimates for each sample
#' @param Chromosome String that contains the name of the column with chromosome denomination for each mutation
#' @param position String that contains the name of the column with position denomination for each mutation
#' @param WT.count String the colname of the column in the data files that contains the number of wild type reads
#' @param mut.count String, the colname of the column in the data files that contains the number of mutant reads
#' @param subclonal.CN String, the name of the total copynumber column
#' @param no.chrs.bearing.mut String, the column that contains the number of chromosomes not bearing any mutations
#' @param mutation.copy.number String, colname with mutation copy number estimates
#' @param subclonal.fraction String, colname of the column that contains the subclonal fraction (or CCF) estimate for each mutation
#' @param is.male Optional boolean that represents the sex, set to TRUE if male, FALSE if female. This information is used to decide whether to include X chromosome mutations
#' @param is.vcf Optional boolean parameter whether the files to be read in are in VCF format
#' @param ref.genome.version Optional string that represents the reference genome, required when reading in VCF files
#' @author sdentro
#' @return A list of tables, one for each type of information
# REMOVED datpath, samplename, data_file_suffix num_muts_sample=NA, 
load.data <- function(list_of_data_files, cellularity, Chromosome, position, WT.count, mut.count, subclonal.CN, no.chrs.bearing.mut, mutation.copy.number, subclonal.fraction, is.male=T, is.vcf=F, ref.genome.version="hg19") {
  data=list()
  
  if (!is.vcf) {
    for(s in 1:length(list_of_data_files)) {
      data[[s]] = read.table(list_of_data_files[s], header=T, stringsAsFactors=F, sep="\t")
    }
  } else {
    for(s in 1:length(list_of_data_files)) {
      v = readVcf(list_of_data_files[s], genome=ref.genome.version)
      # Transform the VCF into the format that the original load.data function understands
      data[[s]] = data.frame(chr=as.vector(seqnames(v)), start=as.vector(start(v))-1, end=as.vector(end(v)),
                             WT.count=as.vector(info(v)$WC), mut.count=as.vector(info(v)$MC), subclonal.CN=as.vector(info(v))$TSC,
                             nMaj1=as.vector(info(v)$NMA1), nMin1=as.vector(info(v)$NMI1), frac1=as.vector(info(v)$FR1),
                             nMaj2=as.vector(info(v)$NMA2), nMin2=as.vector(info(v)$NMI2), frac2=as.vector(info(v)$FR2),
                             phase=as.vector(info(v)$PHS), mutation.copy.number=as.vector(info(v)$MCN), subclonal.fraction=as.vector(info(v)$CCF),
                             no.chrs.bearing.mut=as.vector(info(v)$NCBM))
    }
  }

  # Ofload combining of the tables per sample into a series of tables per data type
  return(load.data.inner(data, cellularity, is.male))
}
  
#' This inner function takes a list of loaded data tables and transforms them into
#' a dataset, which is a list that contains a table per data type
#' @noRD
load.data.inner = function(list_of_tables, cellularity, is.male) {
  no.subsamples = length(list_of_tables)
  no.muts = nrow(list_of_tables[[1]])
  
  # One matrix for each data type and propagate it
  chromosome = matrix(0,no.muts,no.subsamples)
  mut.position = matrix(0,no.muts,no.subsamples)
  WTCount = matrix(0,no.muts,no.subsamples)
  mutCount = matrix(0,no.muts,no.subsamples)
  totalCopyNumber = matrix(0,no.muts,no.subsamples)
  copyNumberAdjustment = matrix(0,no.muts,no.subsamples)
  non.deleted.muts = vector(mode="logical",length=nrow(data[[1]]))
  mutationCopyNumber = matrix(NA,no.muts,no.subsamples)
  subclonalFraction = matrix(NA,no.muts,no.subsamples)
  for(s in 1:length(list_of_data_files)){
    chromosome[,s] = data[[s]][,Chromosome]
    mut.position[,s] = as.numeric(data[[s]][,position])
    WTCount[,s] = as.numeric(data[[s]][,WT.count])
    mutCount[,s] = as.numeric(data[[s]][,mut.count])
    totalCopyNumber[,s] = as.numeric(data[[s]][,subclonal.CN])
    copyNumberAdjustment[,s] = as.numeric(data[[s]][,no.chrs.bearing.mut])
    non.deleted.muts[data[[s]][,no.chrs.bearing.mut]>0]=T
    mutationCopyNumber[,s] = as.numeric(data[[s]][,mutation.copy.number])
    subclonalFraction[,s] = as.numeric(data[[s]][,subclonal.fraction])
  }
  
  # Calculate the kappa, in essense the correction component for the allele frequency of each mutation
  kappa = matrix(1,no.muts,no.subsamples)
  for(i in 1:length(list_of_data_files)){
    #multiply by no.chrs.bearing.mut, so that kappa is the fraction of reads required for fully clonal mutations, rather than mutation at MCN = 1
    kappa[,i] = mutationCopyNumberToMutationBurden(1,data[[i]][,subclonal.CN],cellularity[i]) * data[[i]][,no.chrs.bearing.mut]
  }

  # Remove those mutations that have missing values
  not.there.wt = apply(WTCount, 1, function(x) { sum(is.na(x))>0 })
  not.there.mut = apply(mutCount, 1, function(x) { sum(is.na(x))>0 })
  not.there.cn = apply(totalCopyNumber, 1, function(x) { sum(is.na(x))>0 })
  not.there.cna = apply(copyNumberAdjustment, 1, function(x) { sum(is.na(x))>0 })
  not.there.kappa = apply(kappa, 1, function(x) { sum(is.na(x))>0 })
  # Remove those mutations that have no coverage. These cause for trouble lateron.
  not.coverage = apply(WTCount+mutCount, 1, function(x) { any(x==0) })
  not.cna = apply(copyNumberAdjustment, 1, function(x) { any(x==0) })
  if (is.male) {
    not.on.supported.chrom = apply(chromosome, 1, function(x) { ! any(x %in% as.character(1:22)) })
  } else {
    not.on.supported.chrom = apply(chromosome, 1, function(x) { ! any(x %in% c(1:22, "X")) })
  }

  cov.mean = mean(colMeans(WTCount+mutCount))
  cov.std = mean(apply((WTCount+mutCount), 2, sd))

  too.high.coverage = apply(WTCount+mutCount, 1, function(x) { any(x > cov.mean+6*cov.std)})

  print(paste("Removed", sum(not.there.wt),"with missing WTCount", sep=" "))
  print(paste("Removed", sum(not.there.mut),"with missing mutCount", sep=" "))
  print(paste("Removed", sum(not.there.cn),"with missing totalCopyNumber", sep=" "))
  print(paste("Removed", sum(not.there.cna),"with missing copyNumberAdjustment", sep=" "))
  print(paste("Removed", sum(not.there.kappa),"with missing kappa", sep=" "))
  print(paste("Removed", sum(not.coverage),"with no coverage", sep=" "))
  print(paste("Removed", sum(not.cna),"with zero copyNumberAdjustment", sep=" "))
  print(paste("Removed", sum(not.on.supported.chrom), "on not supported genomic regions", sep=" "))
  print(paste("Removed", sum(too.high.coverage), "mutations with coverage over",cov.mean+6*cov.std, sep=" "))

  select = !(not.there.wt | not.there.mut | not.there.cn | not.there.cna | not.there.kappa | not.coverage | not.cna | not.on.supported.chrom | too.high.coverage)
  
  # Keep indices of removed mutations to 'spike in' lateron when constructing the output
  removed_indices = which(!select)
  chromosome.not.filtered = chromosome
  mut.position.not.filtered = mut.position

  # Remove mutations that have been flagged for various reasons
  chromosome = as.matrix(chromosome[select,])
  mut.position = as.matrix(mut.position[select,])
  WTCount = as.matrix(WTCount[select,])
  mutCount = as.matrix(mutCount[select,])
  totalCopyNumber = as.matrix(totalCopyNumber[select,])
  copyNumberAdjustment = as.matrix(copyNumberAdjustment[select,])
  non.deleted.muts = as.matrix(non.deleted.muts[select])
  kappa = as.matrix(kappa[select,])
  mutationCopyNumber = as.matrix(mutationCopyNumber[select,])
  subclonalFraction = as.matrix(subclonalFraction[select,])
  print(paste("Removed",no.muts-nrow(WTCount), "mutations with missing data"))

#   # Sample mutations
#   if (!is.na(num_muts_sample) & (nrow(chromosome) > 2*num_muts_sample)) {
#   	print(paste("Sampling mutations:", num_muts_sample))
#   	# Store the original mutations
#   	full_data = list(chromosome=chromosome, position=mut.position, WTCount=WTCount, mutCount=mutCount,
#   	                  totalCopyNumber=totalCopyNumber, copyNumberAdjustment=copyNumberAdjustment,
#                       non.deleted.muts=non.deleted.muts, kappa=kappa, mutation.copy.number=mutationCopyNumber,
#   	                  subclonal.fraction=subclonalFraction, removed_indices=removed_indices,
#   		                chromosome.not.filtered=chromosome.not.filtered, mut.position.not.filtered=mut.position.not.filtered,
#   		                sampling.selection=NA, full.data=NA, most.similar.mut=NA)
# 
#   	# Do the sampling
#   	selection = sample(1:nrow(chromosome))[1:num_muts_sample]
#   	selection = sort(selection)
# 
#   	# Select all the data from the various matrices
#   	chromosome = as.matrix(chromosome[selection,])
#   	mut.position = as.matrix(mut.position[selection,])
#   	WTCount = as.matrix(WTCount[selection,])
#   	mutCount = as.matrix(mutCount[selection,])
#   	totalCopyNumber = as.matrix(totalCopyNumber[selection,])
#   	copyNumberAdjustment = as.matrix(copyNumberAdjustment[selection,])
#   	non.deleted.muts = as.matrix(non.deleted.muts[selection,])
#   	kappa = as.matrix(kappa[selection,])
#   	mutationCopyNumber = as.matrix(mutationCopyNumber[selection,])
#   	subclonalFraction = as.matrix(subclonalFraction[selection,])
# 
#   	# for each muation not sampled, find the most similar mutation that was sampled
#   	most.similar.mut = rep(1, nrow(full_data$chromosome))
#   	for (i in 1:nrow(full_data$chromosome)) {
#   		if (i %in% selection) {
# 			# Save index of this mutation within selection - i.e. this row of the eventual mutation assignments must be selected
#   			most.similar.mut[i] = which(selection==i)
#   		} else {
#   			# Find mutation with closest kappa
#   			kappa.diff = matrix(full_data$kappa[selection,]-full_data$kappa[i,], ncol=ncol(full_data$mutCount))
#   			curr = selection[which.min(abs(rowSums(kappa.diff)))]
#   			# Select all mutations with this kappa - a bit of trickery needed to make this work properly with a single column matrix
#   			curr = selection[which(rowSums(matrix(full_data$kappa[selection,], ncol=ncol(full_data$kappa)))==sum(full_data$kappa[curr,]))]
#   			# Pick the mutation with the most similar AF as the the most similar mutation for i
#   			af.i = full_data$mutCount[i,] / (full_data$mutCount[i,] + full_data$WTCount[i,])
#   			af = full_data$mutCount[curr,] / (full_data$mutCount[curr,] + full_data$WTCount[curr,])
#   			af.diff = matrix(af-af.i, ncol=ncol(full_data$mutCount))
#   			curr = curr[which.min(abs(rowSums(af.diff)))]
#   			most.similar.mut[i] = which(selection==curr) # Saving index of most similar mut in the sampled data here for expansion at the end
#   		}
#   	}
  	#print(cbind(1:nrow(full_data$chromosome), most.similar.mut))

#   } else {
	  selection = NA
	  full_data = NA
	  most.similar.mut = NA
#   }
  
  return(list(chromosome=chromosome, position=mut.position, WTCount=WTCount, mutCount=mutCount, 
              totalCopyNumber=totalCopyNumber, copyNumberAdjustment=copyNumberAdjustment, 
              non.deleted.muts=non.deleted.muts, kappa=kappa, mutation.copy.number=mutationCopyNumber, 
              subclonal.fraction=subclonalFraction, removed_indices=removed_indices,
              chromosome.not.filtered=chromosome.not.filtered, mut.position.not.filtered=mut.position.not.filtered,
	            sampling.selection=selection, full.data=full_data, most.similar.mut=most.similar.mut))
}
