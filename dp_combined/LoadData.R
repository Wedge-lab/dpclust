load.data <- function(datpath, samplename, list_of_data_files, cellularity, Chromosome, position, WT.count, mut.count, subclonal.CN, no.chrs.bearing.mut, mutation.copy.number, subclonal.fraction, data_file_suffix, num_muts_sample=NA) {
  # Takes:
  # - An array of data file names to be read in
  # - WT.count, the colname of the column in the data files that contains the number of WT reads
  # - mut.count, the name of the column in the data files that contains the number of mutation bearing reads
  # - subclonal.CN, the name of the total copynumber column
  # - no.chrs.bearing.mut, the column that contains the number of chromosomes not bearing any mutations
  data=list()
  for(s in 1:length(list_of_data_files)){
    data[[s]] = read.table(paste(datpath,samplename,list_of_data_files[s],data_file_suffix,sep=""),header=T,sep="\t",stringsAsFactors=F)
    #data[[s]] = data[[s]][data[[s]][,Chromosome] %in% 1:22,]
  }

  no.subsamples = length(list_of_data_files)
  no.muts = nrow(data[[1]])
  
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
  #not.mut = apply(mutCount, 1, function(x) { any(x<5) })
  not.cna = apply(copyNumberAdjustment, 1, function(x) { any(x==0) })
  not.on.supported.chrom = apply(chromosome, 1, function(x) { ! any(x %in% as.character(1:22)) })

  cov.mean = mean(colMeans(WTCount+mutCount))
  cov.std = mean(apply((WTCount+mutCount), 2, sd))

  too.high.coverage = apply(WTCount+mutCount, 1, function(x) { any(x > cov.mean+6*cov.std)})

  print(paste("Removed", sum(not.there.wt),"with missing WTCount", sep=" "))
  print(paste("Removed", sum(not.there.mut),"with missing mutCount", sep=" "))
  print(paste("Removed", sum(not.there.cn),"with missing totalCopyNumber", sep=" "))
  print(paste("Removed", sum(not.there.cna),"with missing copyNumberAdjustment", sep=" "))
  print(paste("Removed", sum(not.there.kappa),"with missing kappa", sep=" "))
  print(paste("Removed", sum(not.coverage),"with no coverage", sep=" "))
  #print(paste("Removed", sum(not.mut),"with not enough quality coverage of mut", sep=" "))
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

  # Sample mutations
  if (!is.na(num_muts_sample) & (nrow(chromosome) > 2*num_muts_sample)) {
  	print(paste("Sampling mutations:", num_muts_sample))
  	# Store the original mutations
  	full_data = list(chromosome=chromosome, position=mut.position, WTCount=WTCount, mutCount=mutCount,
  	                  totalCopyNumber=totalCopyNumber, copyNumberAdjustment=copyNumberAdjustment,
                      non.deleted.muts=non.deleted.muts, kappa=kappa, mutation.copy.number=mutationCopyNumber,
  	                  subclonal.fraction=subclonalFraction, removed_indices=removed_indices,
  		                chromosome.not.filtered=chromosome.not.filtered, mut.position.not.filtered=mut.position.not.filtered,
  		                sampling.selection=NA, full.data=NA, most.similar.mut=NA)

  	# Do the sampling
  	selection = sample(1:nrow(chromosome))[1:num_muts_sample]
  	selection = sort(selection)

  	# Select all the data from the various matrices
  	chromosome = as.matrix(chromosome[selection,])
  	mut.position = as.matrix(mut.position[selection,])
  	WTCount = as.matrix(WTCount[selection,])
  	mutCount = as.matrix(mutCount[selection,])
  	totalCopyNumber = as.matrix(totalCopyNumber[selection,])
  	copyNumberAdjustment = as.matrix(copyNumberAdjustment[selection,])
  	non.deleted.muts = as.matrix(non.deleted.muts[selection,])
  	kappa = as.matrix(kappa[selection,])
  	mutationCopyNumber = as.matrix(mutationCopyNumber[selection,])
  	subclonalFraction = as.matrix(subclonalFraction[selection,])

  	# for each muation not sampled, find the most similar mutation that was sampled
  	most.similar.mut = rep(1, nrow(full_data$chromosome))
  	for (i in 1:nrow(full_data$chromosome)) {
  		if (i %in% selection) {
			# Save index of this mutation within selection - i.e. this row of the eventual mutation assignments must be selected
  			most.similar.mut[i] = which(selection==i)
  		} else {
  			# Find mutation with closest kappa
  			kappa.diff = matrix(full_data$kappa[selection,]-full_data$kappa[i,], ncol=ncol(full_data$mutCount))
  			curr = selection[which.min(abs(rowSums(kappa.diff)))]
  			# Select all mutations with this kappa - a bit of trickery needed to make this work properly with a single column matrix
  			curr = selection[which(rowSums(matrix(full_data$kappa[selection,], ncol=ncol(full_data$kappa)))==sum(full_data$kappa[curr,]))]
  			# Pick the mutation with the most similar AF as the the most similar mutation for i
  			af.i = full_data$mutCount[i,] / (full_data$mutCount[i,] + full_data$WTCount[i,])
  			af = full_data$mutCount[curr,] / (full_data$mutCount[curr,] + full_data$WTCount[curr,])
  			af.diff = matrix(af-af.i, ncol=ncol(full_data$mutCount))
  			curr = curr[which.min(abs(rowSums(af.diff)))]
  			most.similar.mut[i] = which(selection==curr) # Saving index of most similar mut in the sampled data here for expansion at the end
  		}
  	}
  	#print(cbind(1:nrow(full_data$chromosome), most.similar.mut))

  } else {
	  selection = NA
	  full_data = NA
	  most.similar.mut = NA
  }
  
  return(list(chromosome=chromosome, position=mut.position, WTCount=WTCount, mutCount=mutCount, 
              totalCopyNumber=totalCopyNumber, copyNumberAdjustment=copyNumberAdjustment, 
              non.deleted.muts=non.deleted.muts, kappa=kappa, mutation.copy.number=mutationCopyNumber, 
              subclonal.fraction=subclonalFraction, removed_indices=removed_indices,
              chromosome.not.filtered=chromosome.not.filtered, mut.position.not.filtered=mut.position.not.filtered,
	            sampling.selection=selection, full.data=full_data, most.similar.mut=most.similar.mut))
}
