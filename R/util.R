##############################################################################################################################
# Additional parsers
##############################################################################################################################

#' Helper function to read in the DP algorithm files
#' @param indir The directory in which the files were stored
#' @return A GS.data object with fields S.i, V.h, pi.h and alpha
#' @author sd11
read_gsdata_object = function(indir, samplename, no.iters, conc_param, cluster_conc) {
  GS.data = list()
  filename_prefix = file.path(indir, paste(samplename, "_2D_iters", no.iters, "_concParam", conc_param, "_clusterWidth", 1/cluster_conc, sep=""))
  GS.data$S.i = as.matrix(read.csv(paste(filename_prefix, "_states.csv", sep=""), row.names=1))
  GS.data$V.h = as.matrix(read.csv(paste(filename_prefix, "_stickbreaking_weights.csv", sep=""), row.names=1))
  GS.data$pi.h = as.matrix(read.csv(paste(filename_prefix, "_discreteMutationCopyNumbers.csv", sep=""), row.names=1))
  GS.data$alpha = read.csv(paste(filename_prefix, "_alpha.csv", sep=""), row.names=1)
  
  # Fix the dimensions of this array
  pi.h = array(NA, c(nrow(GS.data$pi.h), ncol(GS.data$pi.h), 1))
  pi.h[1:nrow(GS.data$pi.h), 1:ncol(GS.data$pi.h), 1] = GS.data$pi.h
  GS.data$pi.h = pi.h
  return(GS.data)
}

##############################################################################################################################
# Functions to build mutation coassignment probability matrices
##############################################################################################################################

#' Helper function that builds a mutation to mutation coassignment probability matrix out of mutation to cluster assignments
#' @param GS.data Main clustering output directly from the chain
#' @param density Posterior density through the MCMC output
#' @param no.muts The number of mutations in the data set
#' @param sampledIters Iterations to consider
#' @return A no.muts x no.muts matrix where each cell contains the fraction of iterations the pair of mutations was assigned to the same cluster
build_coassignment_prob_matrix_preferences = function(GS.data, density, no.muts, sampledIters) {
  res = getLocalOptima(density, hypercube.size=5)
  localOptima = res$localOptima
  peak.indices = res$peak.indices
  no.optima = length(localOptima)
  
  boundary = array(NA,no.optima-1)
  mutation.preferences = array(0,c(no.muts,no.optima))
  for(i in 1:(no.optima-1)){
    min.density = min(density$median.density[(peak.indices[i]+1):(peak.indices[i+1]-1)])
    min.indices = intersect(which(density$median.density == min.density),(peak.indices[i]+1):(peak.indices[i+1]-1))
    
    #what distance along the line between a pair of optima do we have to go to reach the minimum density
    boundary[i] = (density$fraction.of.tumour.cells[max(min.indices)] + density$fraction.of.tumour.cells[min(min.indices)])/2
  }
  
  # Adapt this to make a table with the preferred cluster for each mutation in each iteration
  S.i = data.matrix(GS.data$S.i)
  pi.h = GS.data$pi.h[,,1]
  mutation.preferences = array(0, c(no.muts, length(sampledIters)))
  coassignments = array(0, c(no.muts, no.muts))
  for(s in sampledIters) {
    #temp.preferences = array(0,c(no.muts,no.optima))
    for(c in unique(S.i[s,])){
      
      bestOptimum = sum(pi.h[s,c]>boundary)+1
      #temp.preferences[S.i[s,]==c,bestOptimum] = temp.preferences[S.i[s,]==c,bestOptimum] + 1
      assigned.muts = which(S.i[s,]==c)
      coassignments[assigned.muts, assigned.muts] = coassignments[assigned.muts, assigned.muts] + 1
    }
    #iter.preferences = t(apply(temp.preferences, 1, function(p) { as.integer(p==max(p)) / sum(p==max(p)) }))
    #mutation.preferences = mutation.preferences + iter.preferences
  }
  coassignments = coassignments / length(sampledIters)
  return(coassignments)
}

#' Fetch the CCFs of the clusters to which SNVs have been assigned during MCMC into a big table
#' @param pi.h The cluster CCFs
#' @param S.i The mutation assignments to cluster ids
#' @param no.muts Number of total SNVs
#' @param no.iters The total number of iterations
#' @param no.iters.burn.in The total number of iterations to use as burnin
#' @return A matrix with a column for each SNV and a row for each to consider iteration with the cell containing the CCF
#' @author sd11
get_snv_assignment_ccfs = function(pi.h, S.i, no.muts, no.timepoints, no.iters, no.iters.burn.in) {
  #' Get assignment ccfs for each snv
  no.iters.post.burnin = no.iters-no.iters.burn.in
  snv_ccfs = array(NA, c(no.iters.post.burnin, no.muts, no.timepoints))
  x = (no.iters.burn.in+1):no.iters
  for (t in 1:no.timepoints) {
    for (i in 1:no.muts) {
      snv_ccfs[, i, t] = pi.h[cbind(x, S.i[-c(1:no.iters.burn.in), i], t)]
    }
  }
  return(snv_ccfs)
}

#' Helper function that builds the a density over assignment CCFs for each mutation
get_snv_ccf_assignmnent_density = function(S.i, pi.h, no.iters.burn.in, ccf_max_value=10) {
  snv_ccfs = get_snv_assignment_ccfs(pi.h, S.i, ncol(S.i), dim(pi.h)[3], nrow(S.i), no.iters.burn.in)
  snv_ccfs = snv_ccfs[,,1]
  snv_densities = lapply(1:ncol(snv_ccfs), function(i) { 
    if (all(snv_ccfs[,i] < ccf_max_value)) {
      ggplot_build(ggplot(data.frame(ccf=snv_ccfs[,i])) + aes(x=ccf, y=..density..) + geom_density() + xlim(0, ccf_max_value))$data[[1]]$y 
    } else {
      NA
    }
  })
  snv_densities = lapply(snv_densities, function(dat) { dat/sum(dat) })
  return(snv_densities)
}

#' Build a coassignment probability matrix using SNV specific densities over
#' the assigned CCFs during MCMC
#' @param S.i The mutation assignments to cluster ids
#' @param pi.h The cluster CCFs
#' @param no.iters.burn.in The number of iterations to discard as burn in
#' @return Square matrix with a row and column for every mutation, each cell contains the probability of the pair of mutations belonging to the same cluster
#' @author dw9, sd11
build_coassignment_prob_matrix_densities = function(S.i, pi.h, no.iters.burn.in) {
  snv_densities = get_snv_ccf_assignmnent_density(S.i, pi.h, no.iters.burn.in)
  no.muts = length(snv_densities)
  identity.strengths = array(0,c(no.muts,no.muts))
  for (i in 1:(no.muts-1)) {
    identity.strengths[i,i] = 1
    for(j in (i+1):no.muts) {
      if (is.na(snv_densities[[i]]) | is.na(snv_densities[[j]])) {
        identity.strengths[i,j] = identity.strengths[j,i] = 0
      } else {
        identity.strengths[i,j] = identity.strengths[j,i] = 1-(sum(abs(snv_densities[[i]]-snv_densities[[j]])) / 2)
      }
    }
  }
  identity.strengths[no.muts,no.muts] = 1
  return(identity.strengths)
}