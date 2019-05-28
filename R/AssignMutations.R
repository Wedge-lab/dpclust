#
# This file contains various functions to assign mutations to clusters
#

#' Identify clusters and assign mutations for 1-dimensional clustering
#' 
#' @param samplename Sample identifier, used for writing out data to disk
#' @param subclonal.fraction Cancer cell fraction estimates for each mutation
#' @param GS.data Output from the MCMC chain
#' @param density Posterior density estimate of where clusters are located
#' @param no.iters Number of total iterations the MCMC chain was run for
#' @param no.iters.burn.in Number of iterations to discard as burn in
#' @return Standardised mutation clustering output, including clusters, mutation assignments and likelihoods
#' @author dw9, sd11
oneDimensionalClustering <- function(samplename, subclonal.fraction, GS.data, density, no.iters, no.iters.burn.in) {
  no.muts = length(subclonal.fraction)
  normal.copy.number = rep(2,no.muts)
  post.burn.in.start = no.iters.burn.in
  
  S.i = GS.data$S.i
  V.h = GS.data$V.h
  pi.h = GS.data$pi.h[,,1]
  
  # Obtain local optima and peak indices
  res = getLocalOptima(density, hypercube.size=5)
  localOptima = res$localOptima
  peak.indices = res$peak.indices
  write.table(localOptima,paste(samplename,"_localOptima.txt",sep=""), quote=F, sep="\t")
  
  # Assign mutations to clusters
  no.optima = length(localOptima)
  if (no.optima>1) {  
    boundary = array(NA,no.optima-1)
    mutation.preferences = array(0,c(no.muts,no.optima))
    for(i in 1:(no.optima-1)){
      min.density = min(density$median.density[(peak.indices[i]+1):(peak.indices[i+1]-1)])
      min.indices = intersect(which(density$median.density == min.density),(peak.indices[i]+1):(peak.indices[i+1]-1))
      
      #what distance along the line between a pair of optima do we have to go to reach the minimum density
      boundary[i] = (density$fraction.of.tumour.cells[max(min.indices)] + density$fraction.of.tumour.cells[min(min.indices)])/2
    }
    
    sampledIters = (no.iters.burn.in + 1) : no.iters
    #don't use the intitial state
    sampledIters = sampledIters[sampledIters!=1]	
    if(length(sampledIters) > 1000){
      sampledIters=floor(post.burn.in.start + (1:1000) * (no.iters - no.iters.burn.in)/1000)  
    }
    
    S.i = data.matrix(S.i)
    for(s in sampledIters){
      temp.preferences = array(0,c(no.muts,no.optima))
      for(c in unique(S.i[s,])){
        bestOptimum = sum(pi.h[s,c]>boundary)+1
        temp.preferences[S.i[s,]==c,bestOptimum] = temp.preferences[S.i[s,]==c,bestOptimum] + 1		
      }
      iter.preferences = t(apply(temp.preferences, 1, function(p) { as.integer(p==max(p)) / sum(p==max(p)) }))
      mutation.preferences = mutation.preferences + iter.preferences
    }
    mutation.preferences = mutation.preferences / length(sampledIters)
    
    # Drop clusters with all probs for mutations zero
    not.is.empty = !apply(mutation.preferences, 2, function(x) { all(x==0) })
    mutation.preferences = mutation.preferences[, not.is.empty, drop=F]
    localOptima = localOptima[not.is.empty]
    no.optima = length(localOptima)
    
    # Save the cluster assignment probabilities table
    most.likely.cluster = max.col(mutation.preferences)
    out = cbind(mutation.preferences, most.likely.cluster)
    colnames(out)[(ncol(out)-no.optima):ncol(out)] = c(paste("prob.cluster", 1:ncol(mutation.preferences),sep=""),"most.likely.cluster")
    write.table(out, paste(samplename,"_DP_and_cluster_info.txt",sep=""), sep="\t", row.names=F, quote=F)

    # Assemble a table with mutation assignments to each cluster
    cluster_assignment_counts = sapply(1:ncol(mutation.preferences), function(x, m) { sum(m==x) }, m=most.likely.cluster) #table(most.likely.cluster)
    cluster_locations = array(NA, c(length(cluster_assignment_counts), 3))
    cluster_locations[,1] = 1:length(cluster_assignment_counts)
    cluster_locations[,2] = localOptima
    cluster_locations[,3] = cluster_assignment_counts
    
    # Keep a record of all clusters, in case it is of interest
    write.table(cluster_locations, paste(samplename,"_optimaInfo.txt",sep=""), col.names=c("cluster.no","location","no.of.mutations"), row.names=F, sep="\t", quote=F)		
    
    # Clear clusters with no mutations assigned
    non_empty_clusters = which(cluster_locations[,3] > 0)
    cluster_locations = cluster_locations[non_empty_clusters,,drop=F]
    mutation.preferences = mutation.preferences[,non_empty_clusters, drop=F]
    mutation.preferences = mutation.preferences / rowSums(mutation.preferences)
    
    # Sort clusters by CCF
    clust_order = order(cluster_locations[,2], decreasing=T)
    cluster_locations = cluster_locations[clust_order,,drop=F]
    cluster_locations[,1] = 1:nrow(cluster_locations)
    mutation.preferences = mutation.preferences[,clust_order, drop=F]
    
    # get most likely cluster
    most.likely.cluster = max.col(mutation.preferences)
    
    # Obtain likelyhood of most likely cluster assignments
    most.likely.cluster.likelihood = apply(mutation.preferences, 1, max)
    
  }else{
    warning("No local optima found when assigning mutations to clusters")
    most.likely.cluster = rep(1, no.muts)
    most.likely.cluster.likelihood = rep(1, no.muts)
    best.assignment.likelihoods = rep(1, no.muts)
    cluster_locations = NA
  }
  
  return(list(best.node.assignments=most.likely.cluster, best.assignment.likelihoods=most.likely.cluster.likelihood, cluster.locations=cluster_locations, all.assignment.likelihoods=mutation.preferences))
}

#' Assign mutations for multi-dimensional clustering by adding clusters until the likelihood no longer improves
#' 
#' @param GS.data Output from the MCMC chain
#' @param mutCount Matrix containing counts of the mutated allele
#' @param WTCount Matrix containing counts of the wild type allele
#' @param subclonal.fraction Matrix with cancer cell fraction estimates for each mutation in each sample
#' @param node.assignments Matrix with assignments of mutations to clusters during MCMC (i.e. GS.data$S.i)
#' @param opts List with various parameters, including samplename and number of iterations
#' @return Standardised clustering output, including mutation assignments, likelihoods and cluster positions
#' @author dw9, sd11
mutation_assignment_em = function(GS.data, mutCount, WTCount, subclonal.fraction, node.assignments, opts) {
  
  # Unpack analysis options required
  no.iters = opts$no.iters
  no.iters.burn.in = opts$no.iters.burn.in
  no.iters.post.burn.in = opts$no.iters.post.burn.in
  outdir = opts$outdir
  subsamplenames = opts$subsamplenames
  samplename = opts$samplename
  
  # Disabled because it rerutns different values than the original code below
  #identity.strengths = build_coassignment_prob_matrix_densities(GS.data$S.i, GS.data$pi.h, no.iters.burn.in)
  #identity.strengths = identity.strengths*no.iters.post.burn.in
  
  print("Setting up the data")
  no.muts = nrow(mutCount)
  no.subsamples = ncol(mutCount)
  
  if(no.muts<=1){
    stop(paste(samplename," has only 0 or 1 mutations",sep=""))
  }
  
  # Determine mutation strengths across all iterations, discarding burnin
  identity.strengths = array(0,c(no.muts,no.muts))
  for(m in 1:(no.muts-1)){
    identity.strengths[m,m] = no.iters.post.burn.in
    for(n in (m+1):no.muts){
      identity.strengths[m,n] = identity.strengths[n,m] = sum(node.assignments[(1+no.iters-no.iters.post.burn.in):no.iters,m] == node.assignments[(1+no.iters-no.iters.post.burn.in):no.iters,n])
    }
  }
  identity.strengths[no.muts,no.muts] = no.iters-no.iters.post.burn.in
  
  #initialise: assume all mutations are assigned to a single node, with mean subclonal fractions
  likelihoods = 0 
  #subclonal.fraction = mutCount/(mutCount+WTCount)
  #subclonal.fraction[is.nan(subclonal.fraction)]=0
  mean.subclonal.fractions = colMeans(subclonal.fraction)
  
  for(i in 1:no.muts){
    lfoy = log.f.of.y(mutCount[i,], mutCount[i,] + WTCount[i,], rep(1,no.subsamples), mean.subclonal.fractions)
    if(!is.nan(lfoy)){
      likelihoods <- likelihoods + lfoy
    }
  }
  
  print("Opening devices for plotting")
  pdf(file.path(outdir, paste(samplename, "_", no.iters, "iters_", no.iters.burn.in, "burnin_histograms.pdf",sep="")),height=4,width=4*no.subsamples)
  hist.device=dev.cur()
  par(mfrow=c(2,no.subsamples))   
  pdf(file.path(outdir, paste(samplename, "_", no.iters, "iters_", no.iters.burn.in, "burnin_densities.pdf",sep="")),height=4,width=4*no.subsamples)
  density.device=dev.cur()        
  par(mfrow=c(2,no.subsamples))
  if(no.subsamples>1){
    pdf(file.path(outdir, paste(samplename, "_", no.iters, "iters_", no.iters.burn.in, "burnin_consensus_scatter.pdf",sep="")),height=4,width=no.subsamples*(no.subsamples-1)*2)
    scatter.device=dev.cur()
    par(mfrow=c(1,no.subsamples*(no.subsamples-1)/2))
  }
  
  print("Initialising storage")
  consensus.assignments = rep(1,no.muts)
  no.nodes=1
  current.agreement = sum(identity.strengths)
  
  fractional.current.agreement = current.agreement/(no.muts*no.muts*no.iters.post.burn.in)
  print(paste("1 node:",current.agreement,fractional.current.agreement))
  
  #initialise all mutations in one node, so the number of pairwise agreements is just the number of times a pair of mutations appear in the same node
  pairwise.agreements = identity.strengths
  
  all.consensus.assignments = list()
  all.consensus.assignments[[no.nodes]] = consensus.assignments
  all.node.positions = list()
  all.node.positions[[no.nodes]] = array(mean.subclonal.fractions,c(1,no.subsamples))
  all.likelihoods = list()
  all.likelihoods[[no.nodes]] = matrix(rep(1, no.muts*no.subsamples), ncol=no.subsamples)
  
  print("Start adding nodes")
  node.added=T
  while(node.added){      
    unique.nodes = unique(consensus.assignments)
    no.nodes = length(unique.nodes)
    
    new.pairwise.agreements = pairwise.agreements
    new.consensus.assignments = consensus.assignments
    new.node = max(unique.nodes)+1
    new.unique.nodes = c(unique.nodes, new.node)
    
    # iteratively move muts to the new node or back again
    mut.moved=T
    count=1
    saved.consensus.assignments = new.consensus.assignments
    #we may need to avoid infinite cycling by checking whether a set of node assignments has been repeated
    while(mut.moved){
      count=count+1
      mut.moved=F
      rand.inds = sample(no.muts)
      for(r in rand.inds){
        old.agreement = sum(new.pairwise.agreements[r,])
        if(new.consensus.assignments[r]==new.node){
          new.ass = saved.consensus.assignments[r]
        }else{
          new.ass = new.node
        }
        temp.ass = new.consensus.assignments
        temp.ass[r] = new.ass
        new.agreement = sum(identity.strengths[r,temp.ass==new.ass]) + sum(no.iters.post.burn.in - identity.strengths[r,temp.ass!=new.ass])             
        
        if(new.agreement > old.agreement){
          mut.moved=T
          new.consensus.assignments[r] = new.ass
          new.pairwise.agreements[r,]=NA
          new.pairwise.agreements[,r]=NA
          new.pairwise.agreements[r,new.consensus.assignments==new.ass] = identity.strengths[r,new.consensus.assignments==new.ass]
          new.pairwise.agreements[new.consensus.assignments==new.ass,r] = identity.strengths[new.consensus.assignments==new.ass,r]
          new.pairwise.agreements[r,new.consensus.assignments!=new.ass] = no.iters.post.burn.in - identity.strengths[r,new.consensus.assignments!=new.ass]
          new.pairwise.agreements[new.consensus.assignments!=new.ass,r] = no.iters.post.burn.in - identity.strengths[new.consensus.assignments!=new.ass,r]                                        
          old.agreement = new.agreement
        }
      }
    }
    new.agreements = sum(new.pairwise.agreements)
    #don't move a whole node of mutations
    if(length(unique(new.consensus.assignments))<=no.nodes){
      new.agreements = 0
    }
    muts.to.move = which(new.consensus.assignments == new.node)
    if(new.agreements > current.agreement){
      consensus.assignments[muts.to.move] = new.node
      current.agreement = new.agreements
      pairwise.agreements = new.pairwise.agreements
      fractional.current.agreement = current.agreement/(no.muts*no.muts*no.iters.post.burn.in)
      
      #fairly crude - use mean position
      node.position = array(NA,c(no.nodes+1,no.subsamples))
      for(n in 1:(no.nodes+1)){
        for(s in 1:no.subsamples){
          node.position[n,s] = mean(subclonal.fraction[consensus.assignments==n,s])
        }
      }
      
      all.likelihoods[[no.nodes+1]] = calc.new.likelihood(mutCount, mutCount+WTCount, matrix(1, nrow=nrow(mutCount), ncol=ncol(mutCount)), node.position[consensus.assignments,])
      likelihoods = c(likelihoods, sum(all.likelihoods[[no.nodes+1]]))     
      
      all.consensus.assignments[[no.nodes+1]] = consensus.assignments
      all.node.positions[[no.nodes+1]] = node.position
      
      if(no.subsamples>1){
        #its hard to distinguish more than 8 different colours
        max.cols = 8
        cols = rainbow(min(max.cols,no.nodes))
        dev.set(which = scatter.device)
        plot.data = subclonal.fraction
        plot.data[is.na(plot.data)]=0
        for(i in 1:(no.subsamples-1)){
          for(j in (i+1):no.subsamples){
            plot(plot.data[,i],plot.data[,j],type = "n",xlab = paste(samplename,subsamplenames[i]," allele fraction",sep=""), ylab = paste(samplename,subsamplenames[j]," allele fraction",sep=""),xlim = c(0,max(plot.data[,i])*1.2))
            for(n in 1:(no.nodes+1)){
              points(plot.data[,i][consensus.assignments==n],plot.data[,j][consensus.assignments==n],col=cols[(n-1) %% max.cols + 1],pch=20 + floor((n-1)/max.cols))
            }
            #legend(max(plot.data[,i]),max(plot.data[,j]),legend = unique.nodes[n],col=cols[(n-1) %% max.cols + 1],pch=20 + floor((n-1)/max.cols))
            legend(max(plot.data[,i])*1.05,max(plot.data[,j]),legend = 1:(no.nodes+1),col=cols[(0:(no.nodes-1)) %% max.cols + 1],pch=20 + floor((0:(no.nodes-1))/max.cols),cex=1)
          }
        }
      }                                               
    }else{
      node.added = F
    }
  }
  print("Done adding nodes, cleaning up and writing output/last figures")
  dev.off(which = hist.device)
  dev.off(which = density.device)
  if(no.subsamples>1){
    dev.off(which = scatter.device)
  }
  
  tree.sizes = 1:no.nodes
  BIC = bic(likelihoods, no.subsamples, tree.sizes, no.muts)

  best.BIC.index = which.min(BIC)
  # print("likelihoods and BIC")
  # print(cbind(likelihoods,BIC))
  # print(paste("best BIC index=",best.BIC.index,sep=""))
  
  if(no.subsamples>1){
    consensus.assignments = all.consensus.assignments[[best.BIC.index]]
    pdf(file.path(outdir, paste(samplename, "_", no.iters, "iters_", no.iters.burn.in, "burnin_bestScatter.pdf",sep="")),height=4,width=4)
    
    #its hard to distinguish more than 8 different colours
    max.cols = 8
    cols = rainbow(min(max.cols,no.nodes))
    plot.data = subclonal.fraction
    plot.data[is.na(plot.data)]=0
    for(i in 1:(no.subsamples-1)){
      for(j in (i+1):no.subsamples){
        plot(plot.data[,i],plot.data[,j],type = "n",xlab = paste(samplename,subsamplenames[i]," subclonal fraction",sep=""), ylab = paste(samplename,subsamplenames[j]," subclonal fraction",sep=""),xlim = c(0,max(plot.data[,i])*1.25))
        for(n in 1:no.nodes){
          pch=20 + floor((n-1)/max.cols)
          #pch is not implmeneted above 25
          if(pch>25){
            pch=pch-20
          }
          points(plot.data[,i][consensus.assignments==n],plot.data[,j][consensus.assignments==n],col=cols[(n-1) %% max.cols + 1],pch=pch)
        }
        #legend(max(plot.data[,i]),max(plot.data[,j]),legend = unique.nodes[n],col=cols[(n-1) %% max.cols + 1],pch=20 + floor((n-1)/max.cols))
        pch=20 + floor((0:(no.nodes-1))/max.cols)
        pch[pch>25] = pch[pch>25]-20
        legend(max(plot.data[,i])*1.05,max(plot.data[,j]),legend = 1:no.nodes,col=cols[(0:(no.nodes-1)) %% max.cols + 1],pch=pch,cex=1)
      }
    }       
    dev.off()
    
    #
    # Following code commented out because the input for it is not (yet) available. For example of input, see:
    # /lustre/scratch109/sanger/dw9/2D_Dirichlet_Process/Lucy_heterogeneity/AllSampleSubsv0.4Jan23rd2014forDW9.txt
    #
    #     #png(paste("/nfs/team78pc11/dw9/Lucy_heterogeneity_23Jan2014_1000iters/",samplename,"_heterogeneity_linePlot.png",sep=""),width=2000,height=1500)
    #     png(paste(outdir, "/", samplename,"_heterogeneity_linePlot.png",sep=""),width=2000,height=1500)
    #     par(mar=c(10,6,2,2),cex=2)
    #     plot(rep(1:ncol(subclonal.fraction),nrow(subclonal.fraction)),c(subclonal.fraction),type="n",xlab = "sample",xaxt="n",ann = F,xlim=c(0.5,ncol(subclonal.fraction)+2))
    #     axis(1,at=1:ncol(subclonal.fraction),labels=paste(samplename,subsamplenames,sep=""),las=2)
    #     mtext(side = 1, text = "sample", line = 6,cex=3)
    #     mtext(side = 2, text = "allele fraction", line = 4,cex=3)
    #     for(i in 1:nrow(subclonal.fraction)){
    #       linetype = 3
    #       if(data$DRIVER_CATEGORY[i]=="ONCOGENIC"){
    #         linetype=1
    #       }else if(data$DRIVER_CATEGORY[i]=="POSSIBLE_ONCOGENIC"){
    #         linetype=2
    #       }
    #       lines(1:ncol(subclonal.fraction),subclonal.fraction[i,],col=cols[consensus.assignments[i]],lwd=3,lty=linetype)
    #       #points(1:ncol(subclonal.fraction),subclonal.fraction[i,],col=cols[consensus.assignments[i]],pch=20,cex=3)
    #       points(1:ncol(subclonal.fraction),subclonal.fraction[i,],col=cols[consensus.assignments[i]],pch=(i+14)%%25,cex=3)
    #     }
    #     #legend(ncol(subclonal.fraction)+0.2,max(subclonal.fraction),legend = data$geneRef,col=cols[consensus.assignments],pch=20 + floor((0:(no.nodes-1))/max.cols),cex=1)
    #     legend(ncol(subclonal.fraction)+0.2,max(subclonal.fraction),legend = data$geneRef,col=cols[consensus.assignments],pch=(14+(1:nrow(data)))%%25,cex=1)
    #     
    #     dev.off()               
    
  }
  most.likely.cluster = all.consensus.assignments[[best.BIC.index]]
  write.table(cbind(1:length(table(most.likely.cluster)), table(most.likely.cluster), all.node.positions[[best.BIC.index]]), paste(outdir, "/", samplename, "_optimaInfo.txt", sep=""), col.names=c("cluster.no", "no.muts.in.cluster", paste(samplename,subsamplenames,sep="")), sep="\t", quote=F, row.names=F)
  assignment_counts = table(most.likely.cluster)
  cluster.locations = data.frame(cbind(as.numeric(names(assignment_counts)), all.node.positions[[best.BIC.index]], assignment_counts), stringsAsFactors=F)
  return(list(best.node.assignments=most.likely.cluster, best.assignment.likelihoods=all.likelihoods[[best.BIC.index]], cluster.locations=cluster.locations, all.assignment.likelihoods=NA))
}

#' Establish clusters and assign mutations for multi-dimensional clustering using a multi-dimensional density
#' 
#' @param mutation.copy.number Mutation copy number values for all samples
#' @param copyNumberAdjustment Multiplicity values for all samples
#' @param GS.data Output from the MCMC chain containing mutation assignments, cluster positions and cluster weights
#' @param density.smooth Parameter that determines the amount of smoothing applied when establishing multi-dimensional density
#' @param opts List with parameters, including donorname (samplename), individual samplenames (subsamples) and iterations and burnin
#' @return List with various standardised clustering results, including cluster positions, assignments and probabilities
#' @author dw9, sd11
multiDimensionalClustering = function(mutation.copy.number, copyNumberAdjustment, GS.data, density.smooth, opts) {
  #
  # Uses clustering in multi dimensions to obtain a likelihood across all iterations for each mutation
  # The cluster where a mutation is assigned most often is deemed the most likeli destination.
  #
  
  # Unpack the opts
  samplename = opts$samplename
  subsamples = opts$subsamplenames
  no.iters = opts$no.iters
  burn.in = opts$no.iters.burn.in
  new_output_folder = opts$outdir
  
  post.burn.in.start = burn.in
  no.subsamples = length(subsamples)
  no.muts = nrow(mutation.copy.number)
  
  # Get multi-D density
  density.out = Gibbs.subclone.density.est.nd(mutation.copy.number/copyNumberAdjustment,GS.data,density.smooth,burn.in+1,no.iters,max.burden = 1.5)
  
  range = density.out$range
  gridsize = density.out$gridsize
  median.density = density.out$median.density
  lower.CI = density.out$lower.CI
  
  getHypercubeIndices<-function(gridsize,lastMin,hypercube.size){
    indices = array(0,(2*hypercube.size+1)^length(lastMin))
    pos.within.hypercube = array(0,c((2*hypercube.size+1)^length(lastMin),length(lastMin)))
    for(i in 1:length(lastMin)){
      pos.within.hypercube[,i]=rep(0:(2*hypercube.size),each=(2*hypercube.size+1)^(i-1), times = (2*hypercube.size+1)^(length(lastMin)-i))
    }
    indices = pos.within.hypercube[,1] + lastMin[1]
    for(i in 2:length(lastMin)){
      indices = indices + (lastMin[i]+pos.within.hypercube[,i]-1) * prod(gridsize[1:(i-1)])
    }
    return(indices)
  }
  
  hypercube.size = 2
  localMins = array(NA,c(0,no.subsamples))
  lastMin = rep(1,no.subsamples)
  if(median.density[rbind(lastMin+hypercube.size)]>0 & median.density[rbind(lastMin+hypercube.size)] == max(median.density[getHypercubeIndices(gridsize,lastMin,hypercube.size)])){
    localMins = rbind(localMins,lastMin)
  }
  
  getNextHyperCube<-function(gridsize,lastMin,hypercube.size){
    current.dimension = length(gridsize)
    lastMin[current.dimension] = lastMin[current.dimension] + 1
    while(T){
      if(lastMin[current.dimension]==gridsize[current.dimension] - 2 * hypercube.size + 1){
        if(current.dimension==1){
          return(NULL)
        }
        lastMin[current.dimension] = 1
        current.dimension = current.dimension-1
        lastMin[current.dimension] = lastMin[current.dimension] + 1
      }else{
        return(lastMin)
      }
    }
  }
  
  localMins = array(NA,c(0,no.subsamples))
  above95confidence = NULL
  lastMin = rep(1,no.subsamples)
  if(median.density[rbind(lastMin+hypercube.size)]>0 & median.density[rbind(lastMin+hypercube.size)] == max(median.density[getHypercubeIndices(gridsize,lastMin,hypercube.size)])){
    localMins = rbind(localMins,lastMin)
    above95confidence = c(above95confidence,lower.CI[rbind(lastMin+hypercube.size)]>0)
  }
  
  while(!is.null(lastMin)){
    if(median.density[rbind(lastMin+hypercube.size)]>0 & median.density[rbind(lastMin+hypercube.size)] == max(median.density[getHypercubeIndices(gridsize,lastMin,hypercube.size)])){
      localMins = rbind(localMins,lastMin)
      above95confidence = c(above95confidence,lower.CI[rbind(lastMin+hypercube.size)]>0)
    }
    lastMin = getNextHyperCube(gridsize,lastMin,hypercube.size)
  }
  
  localMins = localMins + hypercube.size
  localOptima = array(rep(range[,1],each=no.subsamples),dim(localMins))
  localOptima = localOptima + array(rep((range[,2] - range[,1])/(gridsize-1),each=no.subsamples),dim(localMins)) * localMins
  
  write.table(cbind(localOptima,above95confidence),paste(new_output_folder,"/",samplename,"_localMultidimensionalOptima_",density.smooth,".txt",sep=""),quote=F,sep="\t",row.names=F,col.names = c(paste(samplename,subsamples,sep=""),"above95percentConfidence"))
  write.table(localOptima[above95confidence,],paste(new_output_folder,"/",samplename,"_localHighConfidenceMultidimensionalOptima_",density.smooth,".txt",sep=""),quote=F,sep="\t",row.names=F,col.names = paste(samplename,subsamples,sep=""))
  
  no.optima = nrow(localOptima)
  
  if(no.optima>1){        
    stepsize = (range[,2]-range[,1])/(gridsize-1)
    peak.indices = round((localOptima - rep(range[,1],each = nrow(localOptima))) / rep(stepsize,each = nrow(localOptima)))
    
    #peak.heights are not currently used, but they could be used to construct a mixture model
    #and then estimate posterior probs of belonging to each cluster
    peak.heights = median.density[peak.indices]
    
    #if T, j is 'above' i, relative to the plane through the origin
    vector.direction = array(NA,c(no.optima,no.optima))
    
    boundary = array(NA,c(no.optima,no.optima))
    vector.length = array(NA,c(no.optima,no.optima))
    plane.vector = array(NA,c(no.optima,no.optima,no.subsamples+1))
    for(i in 1:(no.optima-1)){
      for(j in (i+1):no.optima){
        #coefficients describe a plane ax + by + cz + d = 0,
        #perpendicular to the line between optimum i and optimum j, passing through optimum i
        plane.vector[i,j,1:no.subsamples] = localOptima[j,] - localOptima[i,]
        plane.vector[i,j,no.subsamples+1] = -sum(plane.vector[i,j,1:no.subsamples] * localOptima[i,])
        
        #not sure this is needed - it may always be true
        vector.direction[i,j] = sum(plane.vector[i,j,1:no.subsamples] * localOptima[j,]) > sum(plane.vector[i,j,1:no.subsamples] * localOptima[i,])
        #this is needed, in order to normalise the distance
        vector.length[i,j] = sqrt(sum(plane.vector[i,j,1:no.subsamples]^2))
        
        longest.dimension = which.max(abs(plane.vector[i,j,1:no.subsamples])/stepsize)                  
        no.steps = max(abs(plane.vector[i,j,1:no.subsamples])/stepsize) - 1
        #how long is each step?
        Euclidean.stepsize = sqrt(sum((plane.vector[i,j,1:no.subsamples])^2)) / (no.steps + 1)
        
        step.coords = array(NA,c(no.steps,no.subsamples))
        step.coords = t(sapply(1:no.steps,function(o,x){o[i,] + x * (o[j,] - o[i,]) / (no.steps + 1)},o=localOptima))
        step.indices = round((step.coords - rep(range[,1],each = nrow(step.coords))) / rep(stepsize,each = nrow(step.coords)))
        densities.on.line = median.density[step.indices]
        min.indices = which(densities.on.line == min(densities.on.line))
        #what distance along the line between a pair of optima do we have to go to reach the minimum density
        boundary[i,j] = (max(min.indices) + min(min.indices))/2 * Euclidean.stepsize
        #boundary[j,i] = Euclidean.stepsize * (no.steps+1) - boundary[i,j]
      }
    }
    #boundary = boundary - plane.vector[,,no.subsamples + 1] # make distance relative to a plane through the origin
    boundary = boundary - plane.vector[,,no.subsamples + 1] / vector.length # 020714 - normalise adjustment
    
    mutation.preferences = array(0,c(no.muts,no.optima))
    
    sampledIters = (burn.in + 1) : no.iters
    #don't use the intitial state
    sampledIters = sampledIters[sampledIters!=1]            
    if(length(sampledIters) > 1000){
      sampledIters=floor(post.burn.in.start + (1:1000) * (no.iters - burn.in)/1000)                   
    }
    
    S.i = data.matrix(GS.data$S.i)
    for(s in sampledIters){
      temp.preferences = array(0,c(no.muts,no.optima))
      for(c in unique(S.i[s,])){
        for(i in 1:(no.optima-1)){
          for(j in (i+1):no.optima){
            distance.from.plane = sum(GS.data$pi.h[s,c,] * plane.vector[i,j,1:no.subsamples]) / vector.length[i,j]
            if(distance.from.plane<=boundary[i,j]){
              if(vector.direction[i,j])
              {
                temp.preferences[S.i[s,]==c,i] = temp.preferences[S.i[s,]==c,i] + 1
              }else{
                temp.preferences[S.i[s,]==c,j] = temp.preferences[S.i[s,]==c,j] + 1
              }
            }else{
              if(vector.direction[i,j])
              {
                temp.preferences[S.i[s,]==c,j] = temp.preferences[S.i[s,]==c,j] + 1
              }else{
                temp.preferences[S.i[s,]==c,i] = temp.preferences[S.i[s,]==c,i] + 1
              }
            }
          }
        }                       
      }
      iter.preferences = t(sapply(1:no.muts,function(p,i){as.integer(p[i,]==max(p[i,])) / sum(p[i,]==max(p[i,]))},p=temp.preferences))
      mutation.preferences = mutation.preferences + iter.preferences
    }
    mutation.preferences = mutation.preferences / length(sampledIters)
    most.likely.cluster = sapply(1:no.muts,function(m,i){which.max(m[i,])},m=mutation.preferences)
    assignment.likelihood = sapply(1:no.muts,function(m,c,i){ m[i,c[i]] },m=mutation.preferences, c=most.likely.cluster)
    
    #get confidence intervals and median
    #subclonal.fraction = data.matrix(mutation.copy.number/copyNumberAdjustment)
    #no.perms = 10000
    #sampled.vals = array(0,c(no.perms,no.optima,no.subsamples))
    #no.muts.per.cluster = array(0,c(no.perms,no.optima))
    #for(p in 1:no.muts){
    #       print(p)
    #       sampled.cluster = sample(1:no.optima,no.perms,mutation.preferences[p,],replace=T)
    #       for(c in unique(sampled.cluster)){
    #               sampled.vals[sampled.cluster == c,c,] = sampled.vals[sampled.cluster == c,c,] + rep(subclonal.fraction[p,],each = sum(sampled.cluster == c))
    #               no.muts.per.cluster[sampled.cluster == c,c] = no.muts.per.cluster[sampled.cluster == c,c] + 1
    #       }
    #}
    #sampled.vals = sampled.vals/rep(no.muts.per.cluster,times = no.subsamples)
    #quantiles = array(NA, c(no.optima,no.subsamples,3))
    #for(c in 1:no.optima){
    #       for(s in 1:no.subsamples){
    #               quantiles[c,s,] = quantile(sampled.vals[,c,s],probs=c(0.025,0.5,0.975),na.rm=T)
    #       }
    #}
    #new method - 180714
    #quantiles = array(NA, c(no.optima,no.subsamples,3))
    #sampled.thetas = list()
    #totals = table(factor(most.likely.cluster,levels = 1:no.optima))
    #for(i in 1:no.optima){
    #       sampled.thetas[[i]] = array(NA,c(length(sampledIters)*totals[i],no.subsamples))
    #       for(s in 1:length(sampledIters)){
    #               sampled.thetas[[i]][((s-1)*totals[i]+1):(s*totals[i]),] = pi.h[sampledIters[s],S.i[sampledIters[s],most.likely.cluster==i],]
    #       }
    #       for(s in 1:no.subsamples){
    #               quantiles[i,s,] = quantile(sampled.thetas[[i]][,s],probs=c(0.025,0.5,0.975),na.rm=T)
    #       }
    #}
    #new method - 210714 - should be intermediate between previous methods
    quantiles = array(NA, c(no.optima,no.subsamples,3))
    sampled.thetas = list()
    totals = table(factor(most.likely.cluster,levels = 1:no.optima))
    for(i in 1:no.optima){
      sampled.thetas[[i]] = array(NA,c(length(sampledIters),totals[i],no.subsamples))
      for(s in 1:length(sampledIters)){
        sampled.thetas[[i]][s,,] = GS.data$pi.h[sampledIters[s],S.i[sampledIters[s],most.likely.cluster==i],]
      }
      for(s in 1:no.subsamples){
        median.sampled.vals = sapply(1:length(sampledIters),function(x){median(sampled.thetas[[i]][x,,s])})
        quantiles[i,s,] = quantile(median.sampled.vals,probs=c(0.025,0.5,0.975),na.rm=T)
      }
    }
    
    CIs = array(quantiles[,,c(1,3)],c(no.optima,no.subsamples*2))
    CIs = CIs[,rep(1:no.subsamples,each=2) + rep(c(0,no.subsamples),no.subsamples)]
    #out = cbind(data[[1]][,1:2],mutation.preferences,most.likely.cluster)
    #names(out)[(ncol(out)-no.optima):ncol(out)] = c(paste("prob.cluster",1:no.optima,sep=""),"most.likely.cluster")
    out = cbind(mutation.preferences,most.likely.cluster)
    names(out) = c(paste("prob.cluster",1:no.optima,sep=""),"most.likely.cluster")
    
    cluster.locations = cbind(1:ncol(mutation.preferences), quantiles[,,2], colSums(mutation.preferences), table(factor(most.likely.cluster,levels = 1:no.optima)))
    write.table(cluster.locations,
                paste(new_output_folder,"/",samplename,"_optimaInfo_",density.smooth,".txt",sep=""),
                col.names = c("cluster.no", paste(samplename, subsamples, sep=""), "estimated.no.of.mutations", "no.of.mutations.assigned"),
                row.names=F,
                sep="\t",
                quote=F)
    
    write.table(out,paste(new_output_folder,"/",samplename,"_DP_and cluster_info_",density.smooth,".txt",sep=""),sep="\t",row.names=F,quote=F)
    write.table(CIs,paste(new_output_folder,"/",samplename,"_confInts_",density.smooth,".txt",sep=""),col.names = paste(rep(paste(samplename,subsamples,sep=""),each=2),rep(c(".lower.CI",".upper.CI"),no.subsamples),sep=""),row.names=F,sep="\t",quote=F)
  }else{
    most.likely.cluster = rep(1,no.muts)
    warning("No local optima found")
  }
  
  # Remove the estimated number of SNVs per cluster from the cluster locations table  
  cluster.locations = as.data.frame(cluster.locations[,c(1:(length(subsamples)+1),ncol(cluster.locations))])
  # Report only the clusters that have mutations assigned
  non_empty_clusters = cluster.locations[,ncol(cluster.locations)] > 0
  cluster.locations = cluster.locations[non_empty_clusters,,drop=F]
  mutation.preferences = mutation.preferences[,non_empty_clusters, drop=F]
  mutation.preferences = mutation.preferences / rowSums(mutation.preferences)
  
  # Sort clusters by summed CCF (as a proxy for most clonal cluster first)
  clust_order = order(rowSums(cluster.locations[,2:(ncol(cluster.locations)-1)]), decreasing=T)
  cluster.locations = cluster.locations[clust_order,,drop=F]
  cluster.locations[,1] = 1:nrow(cluster.locations)
  mutation.preferences = mutation.preferences[,clust_order]

  # get most likely cluster
  most.likely.cluster = max.col(mutation.preferences)

  # Obtain likelyhood of most likely cluster assignments
  assignment.likelihood = sapply(1:no.muts,function(m,c,i){ m[i,c[i]] },m=mutation.preferences, c=most.likely.cluster)
  
  no.optima = nrow(cluster.locations)
  pdf(paste(new_output_folder,"/",samplename,"_most_likely_cluster_assignment_",density.smooth,".pdf",sep=""),height=4,width=4)
  #its hard to distinguish more than 8 different colours
  max.cols = 8
  cols = rainbow(min(max.cols,no.optima))
  plot.data = mutation.copy.number/copyNumberAdjustment
  plot.data[is.na(plot.data)]=0
  for(i in 1:(no.subsamples-1)){
    for(j in (i+1):no.subsamples){
      plot(plot.data[,i],plot.data[,j],type = "n",xlab = paste(samplename,subsamples[i]," subclonal fraction",sep=""), ylab = paste(samplename,subsamples[j]," subclonal fraction",sep=""),xlim = c(0,max(plot.data[,i])*1.25))
      for(n in 1:no.optima){
        pch=20 + floor((n-1)/max.cols)
        #pch is not implmeneted above 25
        if(pch>25){
          pch=pch-20
        }
        points(plot.data[,i][most.likely.cluster==n],plot.data[,j][most.likely.cluster==n],col=cols[(n-1) %% max.cols + 1],pch=pch)
      }
      #legend(max(plot.data[,i]),max(plot.data[,j]),legend = unique.nodes[n],col=cols[(n-1) %% max.cols + 1],pch=20 + floor((n-1)/max.cols))
      #legend(max(plot.data[,i])*1.05,max(plot.data[,j]),legend = 1:no.optima,col=cols[(0:(no.optima-1)) %% max.cols + 1],pch=20 + floor((0:(no.optima-1))/max.cols),cex=1)
      pch=20 + floor((0:(no.optima-1))/max.cols)
      pch[pch>25] = pch[pch>25]-20
      legend(max(plot.data[,i])*1.05,max(plot.data[,j]),legend = 1:no.optima,col=cols[(0:(no.optima-1)) %% max.cols + 1],pch=pch,cex=1)
      
    }
  }       
  dev.off()

  return(list(best.node.assignments=most.likely.cluster, best.assignment.likelihoods=assignment.likelihood, all.assignment.likelihoods=mutation.preferences, cluster.locations=cluster.locations))
}

#######################################################################################################################
# Binomial based mutation assignment
#######################################################################################################################

#' Note: This doesn't work better than the density based approach
#'
#' Assign mutations to clusters by looking at the binomial probability of each cluster for generating a mutation
#' This for now only works with a single timepoint
#' @noRd
mutation_assignment_binom = function(clustering_density, mutCount, WTCount, copyNumberAdjustment, tumourCopyNumber, normalCopyNumber, cellularity, samplename) {
  # Define convenience variables
  num.timepoints = ncol(mutCount)
  num.muts = nrow(mutCount)

  if (num.timepoints > 1) {
    warning("Assigment of mutations through binomial only implemented for a single timepoint")
    q(save="no")
  }

  # Obtain peak locations whtin the given clustering density
  res = getLocalOptima(clustering_density, hypercube.size=5)
  cluster_locations = res$localOptima

  # Strip out clusters with a small density
  cluster_density = getClusterDensity(clustering_density, cluster_locations, min.window.density=1)
  # Take all clusters with at least 1% of the density
  cluster_locations = cluster_locations[cluster_density > 0.01]
  num.clusters = length(cluster_locations)

  # Calculate log likelihoods for each mutation to be part of each cluster location
  assignment_ll = array(NA, c(num.muts, num.clusters))
  for (t in 1:num.timepoints) {
    for (c in 1:num.clusters) {
      mutBurdens = mutationCopyNumberToMutationBurden(cluster_locations[c] * copyNumberAdjustment[,t], tumourCopyNumber[,t], cellularity[t], normalCopyNumber[,t])
      assignment_ll[,c] = sapply(1:num.muts, function(k, mc, wt, mb) {  mc[k]*log(mb[k]) + wt[k]*log(1-mb[k]) }, mc=mutCount[,t], wt=WTCount[,t], mb=mutBurdens)
    }
  }

  # Convert ll to prob
  assignment_probs = assignment_ll
  assignment_probs[is.na(assignment_probs)] = 0
  assignment_probs = t(apply(assignment_probs, 1, function(assignment_probs_k) { assignment_probs_k - max(assignment_probs_k) }))
  assignment_probs = exp(assignment_probs)
  assignment_probs = matrix(assignment_probs / rowSums(assignment_probs), ncol=num.clusters)

  # Hard assign mutations
  most.likely.cluster = sapply(1:num.muts, function(k, assignment_probs) { which.max(assignment_probs[k,]) }, assignment_probs=assignment_probs)
  assignment.likelihood = sapply(1:num.muts, function(k, assignment_probs, most.likely.cluster) { assignment_probs[k, most.likely.cluster[k]] }, assignment_probs=assignment_probs, most.likely.cluster=most.likely.cluster)

  # Save a table with the output as a summary
  cluster_assignments = table(most.likely.cluster)
  output = array(NA, c(length(cluster_locations), 3))
  
  for (c in 1:num.clusters) {
    cluster_id = names(cluster_assignments)[c]
    cluster_id = as.character(c)

    output[c,1] = as.numeric(cluster_id)
    output[c,2] = cluster_locations[c]
    # Check if there are mutations assigned to the cluster, i.e. it's in the cluster_assignments table
    if (cluster_id %in% names(cluster_assignments)) {
      output[c,3] = cluster_assignments[names(cluster_assignments)==cluster_id]
    } else {
      output[c,3] = 0
    }
  }
  write.table(output, paste(samplename,"_optimaInfo.txt",sep=""), col.names=c("cluster.no","location","no.of.mutations"), row.names=F, sep="\t", quote=F)

  return(list(best.node.assignments=most.likely.cluster, best.assignment.likelihoods=assignment.likelihood, all.assignment.likelihoods=assignment_probs, cluster.locations=output))
}

#' Function that fetches the local optima from a density function call output
#' 
#' @param cluster_density Density estimate of where clusters are located
#' @param hypercube.size Stepsize to use when stepping through the density looking for optima
#' @return A list containing two fields: localOptima, the location of a peak and peak.indices, the index of the peak within a hypercube
#' @author sd11, dw9
getLocalOptima = function(cluster_density, hypercube.size=20) {
  localOptima = NULL
  peak.indices = NULL
  for (i in (1+hypercube.size):(nrow(cluster_density)-hypercube.size)) {
    if (cluster_density$median.density[i] == max(cluster_density$median.density[(i-hypercube.size):(i+hypercube.size)])) {
      localOptima = c(localOptima, cluster_density$fraction.of.tumour.cells[i])
      peak.indices = c(peak.indices, i)
    }
  }
  return(list(localOptima=localOptima, peak.indices=peak.indices))
}

#' Obtain the mean density of each cluster. 
#' 
#' This function takes the cluster_locations and for
#' each cluster it walks from the cluster peak along CCF space until the median density
#' drops below the supplied minimum in both directions. After obtaining the CCF space that a
#' cluster takes up we calculate the mean density in that space. Finally across all cluster
#' densities we normalise to obtain the fraction of total density that each cluster represents
#' @param clustering_density Posterior density estimate of where clusters are
#' @param cluster_locations Estimated cluster locations within the density
#' @param min.window.density Minimum stepsize used
#' @return A vector with for each cluster the fraction of density
#' @author sd11, dw9
getClusterDensity = function(clustering_density, cluster_locations, min.window.density) {
  cluster_density = array(NA, length(cluster_locations))
  for (c in 1:length(cluster_locations)) {
    cluster_location = cluster_locations[c]
    x.cluster = which.min(abs(clustering_density$fraction.of.tumour.cells-cluster_location))
    # Obtain the left most x.axis point of the cluster
    run = T
    i = x.cluster
    while (run) {
      if (i!=0 && clustering_density[i,]$median.density > min.window.density) {
        i = i-1
      } else {
        run = F
      }
    }
    x.cluster.min = i
    
    # Obtain the right most x.axis point of the cluster
    run = T
    i = x.cluster
    while (run) {
      if (i!=(nrow(clustering_density)+1) && clustering_density[i,]$median.density > min.window.density) {
        i = i+1
      } else {
        run = F
      }
    }
    x.cluster.max = i
    
    # Take the average density
    cluster_density[c] = mean(clustering_density[x.cluster.min:x.cluster.max,]$median.density)
  }
  
  # Normalise the densities across the clusters
  cluster_density = cluster_density/sum(cluster_density)
  return(cluster_density)
}

##################################################################
# Confidence intervals and cluster order probabilities
##################################################################

#' Calculate confidence intervals on the cluster location
#' @param GS.data MCMC output with assignments and cluster locations
#' @param mut_assignments Final mutation to cluster assignments
#' @param clusterids Clusterids to run through
#' @param no.muts The total number of mutations
#' @param no.iters Total number of iterations
#' @param no.timepoints Total number of samples in this dataset
#' @param no.iters.burn.in Number of iterations to use as burn-in
#' @return A data.frame with the confidence intervals for each cluster
#' @author sd11
calc_cluster_conf_intervals = function(GS.data, mut_assignments, clusterids, no.muts, no.timepoints, no.iters, no.iters.burn.in) {
  assign_ccfs = get_snv_assignment_ccfs(GS.data$pi.h, GS.data$S.i, no.muts, no.timepoints, no.iters, no.iters.burn.in)
  cluster_intervals = data.frame()
  for (t in 1:no.timepoints) {
    for (i in 1:length(clusterids)) {
      # get all SNVs assigned to this cluster
      clusterid = clusterids[i]
      assigned = which(mut_assignments==clusterid)
      ccfs = assign_ccfs[, assigned, t]
      # Flatten the matrix
      dim(ccfs) = NULL
      quants = t(data.frame(quantile(ccfs, c(.05, .5, .95))))
      cluster_intervals = rbind(cluster_intervals, data.frame(clusterid=clusterid, timepoint=t, quants))
    }
  }
  return(cluster_intervals)
}

#' Get mutation preferences table given a density and cluster locations. The output table
#' contains the CCF of the preferred given cluster locations, i.e. the cluster to which
#' the SNV would've been assigned during clustering if those were the cluster locations
get_mutation_preferences = function(GS.data, density, mut_assignments, clusterids, cluster_ccfs, no.muts, no.timepoints, no.iters, no.iters.burn.in) {
  sampledIters = (no.iters.burn.in+1):no.iters
  
  res = getLocalOptima(density, hypercube.size=5)
  localOptima = res$localOptima
  peak.indices = res$peak.indices
  
  # Check if any corresponding peaks to found clusters
  peak_is_cluster = unlist(lapply(localOptima, function(x) any(abs(x-cluster_ccfs) < .Machine$double.eps ^ 0.5)))
  if (sum(peak_is_cluster)==0) {
    print("No corresponding cluster locations when calculating cluster order probs, this function only works with the density mutation assignment")
    dummy_matrix = array(0, c(length(sampledIters), no.muts, no.timepoints))
    return(dummy_matrix)
  }
  
  # Take only the already found clusters
  peak.indices = peak.indices[peak_is_cluster]
  localOptima = localOptima[peak_is_cluster]
  no.optima = length(localOptima)
  
  S.i = data.matrix(GS.data$S.i)
  pi.h = GS.data$pi.h #[,,1]
  
  if (no.optima == 1) {
    # If only a single optimum was found we assign all muts to that one cluster
    assign_ccfs = array(0, c(length(sampledIters), no.muts, no.timepoints))
    for (t in 1:no.timepoints) {
      for (s in sampledIters) {
        assign_ccfs[s-no.iters.burn.in, 1:no.muts, t] = localOptima[1]
      }
    }
      
  } else {
    
    boundary = array(NA,no.optima-1)
    for(i in 1:(no.optima-1)){
      min.density = min(density$median.density[(peak.indices[i]+1):(peak.indices[i+1]-1)])
      min.indices = intersect(which(density$median.density == min.density),(peak.indices[i]+1):(peak.indices[i+1]-1))
      
      #what distance along the line between a pair of optima do we have to go to reach the minimum density
      boundary[i] = (density$fraction.of.tumour.cells[max(min.indices)] + density$fraction.of.tumour.cells[min(min.indices)])/2
    }
    
    # Get a table with the preferred cluster CCF
    assign_ccfs = array(0, c(length(sampledIters), no.muts, no.timepoints))
    for (t in 1:no.timepoints) {
      for (s in sampledIters) {
        for (c in unique(S.i[s,])) {
          bestOptimum = sum(pi.h[s, c, t]>boundary)+1
          assigned.muts = which(S.i[s,]==c)
          assign_ccfs[s-no.iters.burn.in, assigned.muts, t] = localOptima[bestOptimum] #pi.h[s, c, t]
        }
      }
    }
  }
  return(assign_ccfs)
}


#' Calculate for a pair of clusters whether a has a higher CCF than b.
#' 
#' This function calculates cluster order probabilities by obtaining mutation preferences throughout the MCMC iterations for each provided cluster
#' and then classifies a pair of clusters in groups: Greater than / equal (GT-EQ), less than / equal (LT-EQ), greater than (GT), less than (LT),
#' equal (EQ) or undertain (uncertain). These classifications are obtained by sampling pairs of SNVs from either cluster and account how often
#' SNV 1 is assigned a higher CCF in the preferences than SNV 2.
#' @param GS.data MCMC output with assignments and cluster locations
#' @param clusterids Clusterids to run through
#' @param no.muts The total number of mutations
#' @param no.timepoints Total number of samples in this dataset
#' @param no.iters Total number of iterations
#' @param no.iters.burn.in Number of iterations to use as burn-in
#' @return A array multi-dimensional array with in each cell whether the column cluster has a higher CCF than the row cluster across the samples in the third dimension
#' @author sd11
#' Note: This approach only works with the density based mutation assignment strategy
calc_cluster_order_probs = function(GS.data, density, mut_assignments, clusterids, cluster_ccfs, no.muts, no.timepoints, no.iters, no.iters.burn.in, no.samples=1000) {
  assign_ccfs = get_mutation_preferences(GS.data, density, mut_assignments, clusterids, cluster_ccfs, no.muts, no.timepoints, no.iters, no.iters.burn.in)
  
  num_clusters = length(clusterids)
  sampledIters = no.iters-no.iters.burn.in
  if (num_clusters > 1) {
    
    probs_gt = array(NA, c(length(clusterids), length(clusterids), no.timepoints))
    probs_lt = array(NA, c(length(clusterids), length(clusterids), no.timepoints))
    probs_eq = array(NA, c(length(clusterids), length(clusterids), no.timepoints))
    for (t in 1:no.timepoints) {
      for (c in 1:(num_clusters)) {
        # for (k in (c+1):num_clusters) {
        for (k in 1:(num_clusters)) {
          
          snvs_a = which(mut_assignments==clusterids[c])
          snvs_b = which(mut_assignments==clusterids[k])
          
          if (length(snvs_a)==1) {
            sampled_a = rep(snvs_a, no.samples)
          } else {
            sampled_a = sample(snvs_a, no.samples, replace=T)
          }
          if (length(snvs_b)==1) {
            sampled_b = rep(snvs_b, no.samples)
          } else {
            sampled_b = sample(snvs_b, no.samples, replace=T)
          }
          
          frac_gt = sum(sapply(1:no.samples, function(i) { (sum(assign_ccfs[, sampled_a[i], t] > assign_ccfs[, sampled_b[i], t])) })) / (no.samples*sampledIters)
          frac_lt = sum(sapply(1:no.samples, function(i) { (sum(assign_ccfs[, sampled_a[i], t] < assign_ccfs[, sampled_b[i], t])) })) / (no.samples*sampledIters)
          frac_eq = sum(sapply(1:no.samples, function(i) { (sum(assign_ccfs[, sampled_a[i], t] == assign_ccfs[, sampled_b[i], t])) })) / (no.samples*sampledIters)
          
          # Filling the probability matrices as row-vs-column
          probs_gt[c, k, t] = frac_gt
          probs_lt[c, k, t] = frac_lt
          probs_eq[c, k, t] = frac_eq
        }
      }
    }
    
    #' Should return classification of each cluster/cluster pair (row vs column):
    #'  * GT / LT / EQ / GT-EQ / EQ-LT / uncertain / NA
    classification = array(NA, c(length(clusterids), length(clusterids), no.timepoints))
    classification[probs_gt+probs_eq > 0.95] = "GT-EQ"
    classification[probs_lt+probs_eq > 0.95] = "LT-EQ"
    classification[probs_eq > 0.95] = "EQ"
    classification[probs_lt > 0.95] = "LT"
    classification[probs_gt > 0.95] = "GT"
    classification[probs_lt > 0.95 & probs_gt > 0.95 & probs_eq > 0.95] = "EQ-special"
    classification[is.na(classification)] = "uncertain"
    
    return(list(classification=classification, probs_gt=probs_gt, probs_lt=probs_lt, probs_eq=probs_eq))
  } else {
    dummy_matrix = array(NA, c(length(clusterids), length(clusterids), no.timepoints))
    return(list(classification=dummy_matrix, probs_gt=dummy_matrix, probs_lt=dummy_matrix, probs_eq=dummy_matrix))
  }
}
