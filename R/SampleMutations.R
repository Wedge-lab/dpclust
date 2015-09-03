#' Sample a number of mutations from the dataset to reduce its size. The original
#' dataset object is kept within the returned dataset with label full.data
#' 
#' Note: Resampling an already sampled dataset will not work and returns the original
#' Note2: A conflict array will not be updated.
#' @return A dataset object with only the sampled mutations and a full.data field that contains the original dataset
sample_mutations = function(dataset, num_muts_sample) {
  # Check if sampling already was done
  if (!is.na(dataset$sampling.selection)) {
    return(dataset)
  }
  
  print(paste("Sampling mutations:", num_muts_sample))
  # Store the original mutations
  full_data = list(chromosome=dataset$chromosome, position=dataset$position, WTCount=dataset$WTCount, mutCount=dataset$mutCount,
                   totalCopyNumber=dataset$totalCopyNumber, copyNumberAdjustment=dataset$copyNumberAdjustment,
                   non.deleted.muts=dataset$non.deleted.muts, kappa=dataset$kappa, mutation.copy.number=dataset$mutation.copy.number,
                   subclonal.fraction=dataset$subclonal.fraction, removed_indices=dataset$removed_indices,
                   chromosome.not.filtered=dataset$chromosome.not.filtered, mut.position.not.filtered=dataset$mut.position.not.filtered,
                   sampling.selection=NA, full.data=NA, most.similar.mut=NA, mutationType=dataset$mutationType, cellularity=dataset$cellularity,
                   conflict.array=dataset$conflict.array, phase=dataset$phase)
  
  # Do the sampling - only sample SNVs, leave the CNAs in there (if available)
  selection = sample(1:length(dataset$chromosome[dataset$mutationType=="SNV"]))[1:num_muts_sample]
  selection = c(selection, which(dataset$mutationType=="CNA"))
  selection = sort(selection)
  print(length(selection))
  print(head(dataset$kappa))
  
  # Select all the data from the various matrices
  chromosome = as.matrix(dataset$chromosome[selection,])
  position = as.matrix(dataset$position[selection,])
  WTCount = as.matrix(dataset$WTCount[selection,])
  mutCount = as.matrix(dataset$mutCount[selection,])
  totalCopyNumber = as.matrix(dataset$totalCopyNumber[selection,])
  copyNumberAdjustment = as.matrix(dataset$copyNumberAdjustment[selection,])
  non.deleted.muts = dataset$non.deleted.muts[selection]
  kappa = as.matrix(dataset$kappa[selection,])
  mutation.copy.number = as.matrix(dataset$mutation.copy.number[selection,])
  subclonal.fraction = as.matrix(dataset$subclonal.fraction[selection,])
  mutationType = dataset$mutationType[selection]
  phase = dataset$phase[selection,]
  if (!is.na(dataset$conflict.array)) {
    conflict.array = dataset$conflict.array[selection, selection]
  } else {
    conflict.array = NA
  }

  # Don't update these - maybe this should be done, but not like this as the removed_indices matrix remains the same size as CNAs are added
  #removed_indices = as.matrix(dataset$removed_indices[selection])
  removed_indices = dataset$removed_indices
  
  # for each muation not sampled, find the most similar mutation that was sampled
  # TODO: add in code that at least selects the conflicting SNVs to help find branching trees
  most.similar.mut = rep(1, nrow(full_data$chromosome))
  for (i in 1:length(full_data$chromosome[full_data$mutationType=="SNV"])) {
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
  # Map CNAs back onto themselves
  most.similar.mut = c(most.similar.mut, which(full_data$mutationType=="CNA"))
  
  return(list(chromosome=chromosome, position=position, WTCount=WTCount, mutCount=mutCount, 
              totalCopyNumber=totalCopyNumber, copyNumberAdjustment=copyNumberAdjustment, 
              non.deleted.muts=non.deleted.muts, kappa=kappa, mutation.copy.number=mutation.copy.number, 
              subclonal.fraction=subclonal.fraction, removed_indices=removed_indices,
              chromosome.not.filtered=dataset$chromosome.not.filtered, mut.position.not.filtered=dataset$mut.position.not.filtered,
              sampling.selection=selection, full.data=full_data, most.similar.mut=most.similar.mut,
              mutationType=mutationType, cellularity=dataset$cellularity, conflict.array=conflict.array,
              phase=phase))
}

#' Unsample a sampled dataset and expand clustering results with the mutations that were not used during clustering.
#' Mutations are assigned to the same cluster as the most similar mutation is assigned to
#' @return A list containing the unsampled dataset and the clustering object with the unused mutations included
unsample_mutations = function(dataset, clustering_result) {
  best.node.assignments = rep(1, nrow(dataset$full.data$chromosome))
  best.assignment.likelihoods = rep(1, nrow(dataset$full.data$chromosome))
  best.node.assignments = clustering_result$best.node.assignments[dataset$most.similar.mut]
  best.assignment.likelihoods = clustering_result$best.assignment.likelihoods[dataset$most.similar.mut]
  clustering = list(best.node.assignments=best.node.assignments, best.assignment.likelihoods=best.assignment.likelihoods)
  return(list(dataset=dataset$full.data, clustering=clustering))
}

