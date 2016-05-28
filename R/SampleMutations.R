#' Sample a number of mutations from the dataset to reduce its size. The original
#' dataset object is kept within the returned dataset with label full.data
#' 
#' Note: Resampling an already sampled dataset will not work and returns the original
#' Note2: A conflict array will not be updated.
#' @param min_sampling_factor num_muts_sample*min_sampling_factor is the minimum number of mutations to have before sampling is applied. Use this multiplier to make sure we're not just sampling out a very low fraction of mutations (Default: 1.5)
#' @param sampling_method Integer selecting a sampling method. 1 is uniform sampling, 2 is a hybrid of uniform and bin sampling
#' @return A dataset object with only the sampled mutations and a full.data field that contains the original dataset
sample_mutations = function(dataset, num_muts_sample, min_sampling_factor=1.5, sampling_method=1, sample.snvs.only=T, remove.snvs=F) {

  # Check if sampling already was done
  if (!is.na(dataset$sampling.selection)) {
    return(dataset)
  }
  
  print(paste("Sampling mutations:", num_muts_sample))
  
  # Get the data that is available for sampling
  avail_for_sampling = get_data_avail_for_sampling(dataset, sample.snvs.only, remove.snvs)
  
  # Check that the amount of data available for sampling is sufficient, if not, return the original dataset
  if (length(avail_for_sampling) < floor(min_sampling_factor*num_muts_sample)) {
    print(paste("Num muts smaller than", min_sampling_factor, "*threshold, not performing sampling", sep=""))
    return(dataset)
  }
  
  # Store the original mutations
  full_data = list(chromosome=dataset$chromosome, position=dataset$position, WTCount=dataset$WTCount, mutCount=dataset$mutCount,
                   totalCopyNumber=dataset$totalCopyNumber, copyNumberAdjustment=dataset$copyNumberAdjustment,
                   non.deleted.muts=dataset$non.deleted.muts, kappa=dataset$kappa, mutation.copy.number=dataset$mutation.copy.number,
                   subclonal.fraction=dataset$subclonal.fraction, removed_indices=dataset$removed_indices,
                   chromosome.not.filtered=dataset$chromosome.not.filtered, mut.position.not.filtered=dataset$mut.position.not.filtered,
                   sampling.selection=NA, full.data=NA, most.similar.mut=NA, mutationType=dataset$mutationType, cellularity=dataset$cellularity,
                   conflict.array=dataset$conflict.array, phase=dataset$phase)
  
  
  if (sampling_method==1) {
    # Perform regular uniform sampling
    selection = do_uniform_sampling(avail_for_sampling, num_muts_sample)
    
  } else if (sampling_method==2) {
    # Perform sampling of mutations per bin - 1/5th of the total
    selection_bins = do_ccf_bin_sampling(dataset$subclonal.fraction[avail_for_sampling], floor(num_muts_sample/5), 250, 0.05, max_ccf_bin=0.9)
    if (length(selection_bins) < num_muts_sample) {
      # Sample remaining mutations randomly
      num_to_sample_extra = num_muts_sample-length(selection_bins)
      selection_random = do_uniform_sampling(avail_for_sampling, num_muts_sample)
      # Remove mutations already selected through bin sampling
      selection_random = selection_random[!(selection_random %in% selection_bins)]
      selection = sort(c(selection_bins, selection_random[1:num_to_sample_extra]))
    } else {
      selection = selection_bins
    }
    
  } else if (sampling_method==3) {
    # Take only subclonal mutations
    if (ncol(dataset$mutCount)==1) {
      print("Taking only subclonal data. This only works in single sample cases")
      selection = avail_for_sampling[which(dataset$subclonal.fraction[avail_for_sampling, 1] < 0.9)]
    } else {
      print("Taking only subclonal data does not work for multi-sample cases, returning all data available for sampling")
      selection = avail_for_sampling
    }
    
  } else {
    print("Unsupported sampling method supplied. No sampling performed.")
    return(dataset)
  }
  print(paste("Subsampled number of mutations: ", length(selection)))
  
  # Add CNA pseudo-SNVs if these were not to be sampled
  if (sample.snvs.only & !remove.snvs) {
    selection = c(selection, which(dataset$mutationType=="CNA"))
  }
  
  
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
  
  for (i in 1:nrow(full_data$chromosome)) { #[full_data$mutationType=="SNV"]
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
  # # Map CNAs back onto themselves
  # most.similar.mut = c(most.similar.mut, which(full_data$mutationType=="CNA"))
  
  return(list(chromosome=chromosome, position=position, WTCount=WTCount, mutCount=mutCount, 
              totalCopyNumber=totalCopyNumber, copyNumberAdjustment=copyNumberAdjustment, 
              non.deleted.muts=non.deleted.muts, kappa=kappa, mutation.copy.number=mutation.copy.number, 
              subclonal.fraction=subclonal.fraction, removed_indices=dataset$removed_indices,
              chromosome.not.filtered=dataset$chromosome.not.filtered, mut.position.not.filtered=dataset$mut.position.not.filtered,
              sampling.selection=selection, full.data=full_data, most.similar.mut=most.similar.mut,
              mutationType=mutationType, cellularity=dataset$cellularity, conflict.array=conflict.array,
              phase=phase))
}

#' Unsample a sampled dataset and expand clustering results with the mutations that were not used during clustering.
#' Mutations are assigned to the same cluster as the most similar mutation is assigned to
#' @dataset A dataset object in which mutations have been sampled
#' @clustering_result A clustering result object based on the downsampled mutations
#' @return A list containing the unsampled dataset and the clustering object with the unused mutations included
unsample_mutations = function(dataset, clustering_result) {
  # Update the cluster summary table with the new assignment counts
  best.node.assignments = clustering_result$best.node.assignments[dataset$most.similar.mut]
  cluster.locations = clustering_result$cluster.locations
  new_assignment_counts = table(best.node.assignments)
  for (cluster in names(new_assignment_counts)) {
    cluster.locations[cluster.locations[,1]==as.numeric(cluster), 3] = new_assignment_counts[cluster]
  }
  
  clustering = list(all.assignment.likelihoods=clustering_result$all.assignment.likelihoods[dataset$most.similar.mut,],
                    best.node.assignments=best.node.assignments, 
                    best.assignment.likelihoods=clustering_result$best.assignment.likelihoods[dataset$most.similar.mut],
                    cluster.locations=cluster.locations)
  # Save the cndata, if available
  new_dataset = dataset$full.data
  new_dataset$cndata = dataset$cndata
  return(list(dataset=new_dataset, clustering=clustering))
}

#' Function that returns the indices of data that is available for sampling given
#' that one can choose to only sample SNVs, sample both SNVs and CNAs or remove SNVs
#' alltogether
#' @param dataset A dataset object
#' @param sample.snvs.only Boolean that if set to TRUE samples SNVs only. CNA indices should then lateron be added manually if they need to stay in the dataset
#' @param remove.snvs A boolean that if set to TRUE removes the SNVs from the dataset and returns just the CNA indices
#' @return A list of indices representing the data that can be sampled under the given restrictions
#' @author sd11
get_data_avail_for_sampling = function(dataset, sample.snvs.only, remove.snvs) {
  # Make inventory of what can be sampled
  if (sample.snvs.only & !remove.snvs) {
    print("Sampling only SNVs")
    avail_for_sampling = which(dataset$mutationType=="SNV")
    
  } else if (remove.snvs) {
    print("Sampling only CNAs, removing all SNVs")
    # Remove SNVs and sample CNAs
    avail_for_sampling = which(dataset$mutationType=="CNA")
    
  } else {
    print("Sampling all data")
    avail_for_sampling = 1:nrow(dataset$chromosome)
  }
  return(avail_for_sampling)
}

#' Uniformly sample mutations
#' @param num_muts The total number of mutations in the dataset
#' @param num_muts_sample The number of mutations to sample
#' @return A vector with indices of mutations to select
#' @author sd11
do_uniform_sampling = function(avail_for_sampling, num_muts_sample) {
  selection = sample(avail_for_sampling)[1:num_muts_sample]
  selection = sort(selection)
  return(selection)
}

#' Sample bins in CCF space, then sample mutations from each bin - This does not work very well
#' as it often finds incorrect clusters due to the uneven sampling from both sides of the clonal
#' cluster. Seems to become worse with larger clonal peaks
do_ccf_bin_sampling = function(dat, num_muts_sample, max_bin_selection, bin_size, max_ccf_bin=NA) {
  selection = c()
  # dat contains CCF
  min_ccf = min(dat)
  if (is.na(max_ccf_bin)) {
    max_ccf = max(dat)
  } else {
    max_ccf = max_ccf_bin
  }
  bins = cut(dat[,1], seq(min_ccf, max_ccf, bin_size))
  bin_counts = table(bins)
  
  # If the number of mutations to sample is larger than what is available in the bins, return all mutations
  if (sum(bin_counts) <= num_muts_sample) {
    print("do_ccf_bin_sampling: Number of mutations available smaller than number to sample, returning all")
    return(1:nrow(dat))
  }
  
  bin_selection = sample(1:length(bin_counts), num_muts_sample, replace=T)
  bin_selection = table(bin_selection)
  
  # If a bin has more selections than is allowed, redistribute some selections to other bins
  bins_over_max = c()
  if (any(bin_selection > max_bin_selection)) {
    # Resample the number of mutations each of the bins is over the threshold
    bins_over_max = which(bin_selection > max_bin_selection)
    if (length(bins_over_max) > 0) {
      resamples = c()
      for (i in bins_over_max) {
        num_to_resample = bin_selection[i] - max_bin_selection
        # Remove the bins that are over the max so we don't accidentally select them again
        bin_ids = 1:length(bin_counts)[-bins_over_max]
        resamples = c(resamples, sample(bin_ids, num_to_resample))
        # Set the selection to the allowed max
        bin_selection[i] = max_bin_selection
      }
    }
    
    # Add the resamples to the table
    for (selection in resamples) {
      bin_selection[selection] = bin_selection[selection] + 1
    }
  }
  
  # Check if any bin has been selected more often then there are mutations in it
  # redistribute selections if this is the case
  num_to_resample = 0
  bin_oversampled = T
  max.iters = 10
  iter = 1
  print("Redistributing selections from oversampled bins")
  while (bin_oversampled && iter<=max.iters) {
    print(iter)
    for (i in names(bin_selection)) {
      if (bin_selection[i] > bin_counts[as.numeric(i)]) {
        num_to_resample = num_to_resample + (bin_selection[i] - bin_counts[as.numeric(i)])
        bins_over_max = c(bins_over_max, as.numeric(i))
        # Set selection to total number of mutations in bin
        bin_selection[i] = bin_counts[as.numeric(i)]
      }
    }
    if (num_to_resample > 0) {
      bin_ids = (1:length(bin_counts))[-bins_over_max]
      resamples = sample(bin_ids, num_to_resample, replace=T, prob=bin_counts[bin_ids]/sum(bin_counts[bin_ids]))
      # Add the resamples to the table
      for (selection in resamples) {
        bin_selection[selection] = bin_selection[selection] + 1
      }
    } else {
      bin_oversampled = F
    }
    num_to_resample = 0
    iter = iter + 1
  }
  
  print("Bin counts")
  print(rbind(bin_counts, bin_selection))
  
  selection = c()
  # Select the determined number of mutations from each bin
  for (index in as.numeric(names(bin_selection))) {
    bin_name = names(bin_counts)[index]
    num_sample_bin = bin_selection[index]
    selection = c(selection, sample(which(bins==bin_name), num_sample_bin))
  }
  
  return(sort(selection))
}
