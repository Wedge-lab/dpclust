
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
  
  print("localOptima")
  print(localOptima)
  
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
    # Clear clusters with no mutations assigned
    cluster_locations = cluster_locations[cluster_locations[,3] > 0,,drop=F]
    print(cluster_locations)
    write.table(cluster_locations, paste(samplename,"_optimaInfo.txt",sep=""), col.names=c("cluster.no","location","no.of.mutations"), row.names=F, sep="\t", quote=F)		
    
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


# twoDimensionalClustering <- function(samplename, subclonal.fraction, GS.data, density, no.iters, no.iters.burn.in) {
#   
#   hypercube.size = 20
#   localOptima.x = NULL
#   localOptima.y = NULL
#   peak.indices.x = NULL
#   peak.indices.y = NULL
#   # Obtain density for each point in the data.frame
#   for(i in (1+hypercube.size):(nrow(density$median.density)-hypercube.size)){
#     for(j in (1+hypercube.size):(ncol(density$median.density)-hypercube.size)){
#       if(density$median.density[i,j] == max(density$median.density[(i-hypercube.size):(i+hypercube.size), (j-hypercube.size):(j+hypercube.size)])){
#         # Only take those clusters that have a non-trivial density
#         if(density$median.density[i,j] > 0.1) { 
#           localOptima.x = c(localOptima.x,density$xvals[i])
#           localOptima.y = c(localOptima.y,density$yvals[j])
#           peak.indices.x = c(peak.indices.x,i)
#           peak.indices.y = c(peak.indices.y,j)
#         }
#       }
#     }
#   }
# 
#   localOptima = data.frame(x=localOptima.x, y=localOptima.y)
#   print("localOptima")
#   print(localOptima)
#   
#   write.table(localOptima,paste(samplename,"_localOptima.txt",sep=""),quote=F,sep="\t")
#   
#   # Below is not adapted to 2D case
# #   no.optima = nrow(localOptima)
# #   if(no.optima>1){  
# #     boundary = array(NA,no.optima-1)
# #     mutation.preferences = array(0,c(no.muts,no.optima))
# #     for(i in 1:(no.optima-1)){
# #       min.density = min(density$median.density[(peak.indices.x[i]+1):(peak.indices.x[i+1]-1), (peak.indices.y[i]+1):(peak.indices.y[i+1]-1)])
# #       min.indices = intersect(which(density$median.density == min.density),(peak.indices[i]+1):(peak.indices[i+1]-1))
# #       
# #       #what distance along the line between a pair of optima do we have to go to reach the minimum density
# #       boundary[i] = (density$fraction.of.tumour.cells[max(min.indices)] + density$fraction.of.tumour.cells[min(min.indices)])/2
# #     }
# #     
# #     sampledIters = (no.iters.burn.in + 1) : no.iters
# #     #don't use the intitial state
# #     sampledIters = sampledIters[sampledIters!=1]	
# #     if(length(sampledIters) > 1000){
# #       sampledIters=floor(post.burn.in.start + (1:1000) * (no.iters - no.iters.burn.in)/1000)  
# #     }
# #     
# #     S.i = data.matrix(S.i)
# #     for(s in sampledIters){
# #       temp.preferences = array(0,c(no.muts,no.optima))
# #       for(c in unique(S.i[s,])){
# #         bestOptimum = sum(pi.h[s,c]>boundary)+1
# #         temp.preferences[S.i[s,]==c,bestOptimum] = temp.preferences[S.i[s,]==c,bestOptimum] + 1		
# #       }
# #       iter.preferences = t(apply(temp.preferences, 1, function(p) { as.integer(p==max(p)) / sum(p==max(p)) }))
# #       mutation.preferences = mutation.preferences + iter.preferences
# #     }
# #     mutation.preferences = mutation.preferences / length(sampledIters)
# #     most.likely.cluster = max.col(mutation.preferences)
# #     
# #     out = cbind(mutation.preferences,most.likely.cluster)
# #     colnames(out)[(ncol(out)-no.optima):ncol(out)] = c(paste("prob.cluster",1:no.optima,sep=""),"most.likely.cluster")
# #     write.table(cbind(1:no.optima,localOptima,colSums(mutation.preferences)),paste(samplename,"_optimaInfo.txt",sep=""),col.names = c("cluster.no","MCN","no.of.mutations"),row.names=F,sep="\t",quote=F)		
# #     write.table(out,paste(samplename,"_DP_and_cluster_info.txt",sep=""),sep="\t",row.names=F,quote=F)
# #     
# #     most.likely.cluster.likelihood = apply(mutation.preferences, 1, max)
# #     
# #   }else{
# #     most.likely.cluster = rep(1,no.muts)
# #     most.likely.cluster.likelihood = rep(1,no.muts)
# #   }
# }


# q(save="no")

mutation_assignment_em = function(GS.data, mutCount, WTCount, subclonal.fraction, node.assignments, opts) {
  
  # Unpack analysis options required
  no.iters = opts$no.iters
  no.iters.burn.in = opts$no.iters.burn.in
  no.iters.post.burn.in = opts$no.iters.post.burn.in
  outdir = opts$outdir
  subsamplenames = opts$subsamplenames
  samplename = opts$samplename
  
  identity.strengths = build_coassignment_prob_matrix_densities(GS.data$S.i, GS.data$pi.h, no.iters.burn.in)
  identity.strengths = identity.strengths*no.iters.post.burn.in
  
  
  print("Setting up the data")
  no.muts = nrow(mutCount)
  no.subsamples = ncol(mutCount)
  
  if(no.muts<=1){
    print(paste(samplename," has only 0 or 1 mutations",sep=""))
    next()
  }
  
  # # Determine mutation strengths across all iterations, discarding burnin
  # identity.strengths = array(0,c(no.muts,no.muts))
  # for(m in 1:(no.muts-1)){
  #   identity.strengths[m,m] = no.iters.post.burn.in
  #   for(n in (m+1):no.muts){
  #     identity.strengths[m,n] = identity.strengths[n,m] = sum(node.assignments[(1+no.iters-no.iters.post.burn.in):no.iters,m] == node.assignments[(1+no.iters-no.iters.post.burn.in):no.iters,n])
  #   }
  # }
  # identity.strengths[no.muts,no.muts] = no.iters-no.iters.post.burn.in
  
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
  pdf(paste(outdir, "/", samplename, "_", no.iters, "iters_", no.iters.burn.in, "burnin_histograms.pdf",sep=""),height=4,width=4*no.subsamples)
  hist.device=dev.cur()
  par(mfrow=c(2,no.subsamples))   
  pdf(paste(outdir, "/", samplename, "_", no.iters, "iters_", no.iters.burn.in, "burnin_densities.pdf",sep=""),height=4,width=4*no.subsamples)
  density.device=dev.cur()        
  par(mfrow=c(2,no.subsamples))
  if(no.subsamples>1){
    pdf(paste(outdir, "/", samplename, "_", no.iters, "iters_", no.iters.burn.in, "burnin_consensus_scatter.pdf",sep=""),height=4,width=no.subsamples*(no.subsamples-1)*2)
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
  
  print("Starting EM")
  node.added=T
  while(node.added){      
    unique.nodes = unique(consensus.assignments)
    no.nodes = length(unique.nodes)
    print(no.nodes)
    new.pairwise.agreements = pairwise.agreements
    new.consensus.assignments = consensus.assignments
    new.node = max(unique.nodes)+1
    new.unique.nodes = c(unique.nodes, new.node)
    
    #EM algorithm - iteratively move muts to the new node or back again
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
      print(paste(new.node,current.agreement,fractional.current.agreement))
      
      #fairly crude - use mean position
      node.position = array(NA,c(no.nodes+1,no.subsamples))
      for(n in 1:(no.nodes+1)){
        for(s in 1:no.subsamples){
          node.position[n,s] = mean(subclonal.fraction[consensus.assignments==n,s])
        }
      }
      
      all.likelihoods[[no.nodes+1]] = calc.new.likelihood2(mutCount, mutCount+WTCount, matrix(1, nrow=nrow(mutCount), ncol=ncol(mutCount)), node.position[consensus.assignments,])
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
  print("Done EM, cleaning up and writing output/last figures")
  dev.off(which = hist.device)
  dev.off(which = density.device)
  if(no.subsamples>1){
    dev.off(which = scatter.device)
  }
  
  tree.sizes = 1:no.nodes
  BIC = bic(likelihoods, no.subsamples, tree.sizes, no.muts)

  best.BIC.index = which.min(BIC)
  print("likelihoods and BIC")
  print(cbind(likelihoods,BIC))
  print(paste("best BIC index=",best.BIC.index,sep=""))
  
  #
  # Following code commented out because the input for it is not (yet) available. For example of input, see:
  # /lustre/scratch109/sanger/dw9/2D_Dirichlet_Process/Lucy_heterogeneity/AllSampleSubsv0.4Jan23rd2014forDW9.txt
  #
  #   write.table(cbind(data[,c(2,3,4,5,6,9,10,11)],subclonal.fraction,all.consensus.assignments[[best.BIC.index]]),paste(outdir, "/", samplename,"_muts_withDPAssignments.txt",sep=""),col.names = c(names(data)[c(2,3,4,5,6,9,10,11)],colnames(subclonal.fraction),"DirichletProcessCluster"),sep="\t",quote=F,row.names=F)
  
  if(no.subsamples>1){
    consensus.assignments = all.consensus.assignments[[best.BIC.index]]
    pdf(paste(outdir, "/", samplename, "_", no.iters, "iters_", no.iters.burn.in, "burnin_bestScatter.pdf",sep=""),height=4,width=4)
    
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
  cluster.locations = data.frame(cbind(names(assignment_counts), all.node.positions[[best.BIC.index]], assignment_counts))

  return(list(best.node.assignments=most.likely.cluster, best.assignment.likelihoods=all.likelihoods[[best.BIC.index]], cluster.locations=cluster.locations, all.assignment.likelihoods=NA))
}

multiDimensionalClustering = function(mutation.copy.number, copyNumberAdjustment, GS.data, density.smooth, opts) {
  #
  # Uses clustering in multi dimensions to obtain a likelihood across all iterations for each mutation
  # The cluster where a mutation is assigned most often is deemed the most likeli destination.
  #
  
  # Unpack the opts
  samplename = opts$samplename
  subsamples = opts$subsamplenames
  noiters = opts$no.iters
  burn.in = opts$no.iters.burn.in
  new_output_folder = opts$outdir
  
  no.subsamples = length(subsamples)
  
  # Get multi-D density
  density.out = Gibbs.subclone.density.est.nd(mutation.copy.number/copyNumberAdjustment,GS.data,density.smooth,burn.in+1,noiters,max.burden = 1.5)
  
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
      print("local maximum indices:")
      print(lastMin)
      above95confidence = c(above95confidence,lower.CI[rbind(lastMin+hypercube.size)]>0)
    }
    lastMin = getNextHyperCube(gridsize,lastMin,hypercube.size)
  }
  
  localMins = localMins + hypercube.size
  localOptima = array(rep(range[,1],each=no.subsamples),dim(localMins))
  localOptima = localOptima + array(rep((range[,2] - range[,1])/(gridsize-1),each=no.subsamples),dim(localMins)) * localMins
  
  print("localMins")
  print(localMins)
  print("localOptima")
  print(localOptima)
  
  write.table(cbind(localOptima,above95confidence),paste(new_output_folder,"/",samplename,"_localMultidimensionalOptima_",density.smooth,".txt",sep=""),quote=F,sep="\t",row.names=F,col.names = c(paste(samplename,subsamples,sep=""),"above95percentConfidence"))
  print(localOptima)
  print(localOptima[above95confidence,])
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
    
    no.muts = nrow(mutation.copy.number)
    
    mutation.preferences = array(0,c(no.muts,no.optima))
    
    sampledIters = (burn.in + 1) : noiters
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

  # Remove the estimated number of SNVs per cluster from the cluster locations table  
  cluster.locations = as.data.frame(cluster.locations[,c(1:(length(subsamples)+1),ncol(cluster.locations))])
  # Report only the clusters that have mutations assigned
  cluster.locations = cluster.locations[cluster.locations[,ncol(cluster.locations)] > 0,]
  
  return(list(best.node.assignments=most.likely.cluster, best.assignment.likelihoods=assignment.likelihood, all.assignment.likelihoods=mutation.preferences, cluster.locations=cluster.locations))
}

# #' Function that creates various final output files that contain fixed information about mutations.
# #' This includes adding removed mutations (either due to sampling or data characteristics) back in.
# produceMutAssignmentOutput = function(dataset, clustering, outfiles.prefix, most.similar.mut=NA, write_tree=F) {
#   # Check if mutation sampling has been done, if so, unpack and assign here
#   if (!is.na(most.similar.mut)) {
#     res = unsample_mutations(dataset, clustering)
#     dataset = res$dataset
#     clustering = res$clustering
#   }
#   
#   # Create final output data matrix
#   output = cbind(dataset$chromosome[,1], dataset$position[,1]-1, dataset$position[,1], clustering$best.node.assignments, clustering$best.assignment.likelihoods)
#   
#   # Add the removed mutations back in
#   for (i in dataset$removed_indices) {
#     if (i==1) {
#       output = rbind(c(dataset$chromosome.not.filtered[i], dataset$mut.position.not.filtered[i]-1, dataset$mut.position.not.filtered[i], NA, NA), output)
#     } else if (i >= nrow(output)) {
#       output = rbind(output, c(dataset$chromosome.not.filtered[i], dataset$mut.position.not.filtered[i]-1, dataset$mut.position.not.filtered[i], NA, NA))
#     } else {
#       output = rbind(output[1:(i-1),], c(dataset$chromosome.not.filtered[i], dataset$mut.position.not.filtered[i]-1, dataset$mut.position.not.filtered[i], NA, NA), output[i:nrow(output),])
#     }
#   }
#   
#   # Save the indices of the mutations that were not used during the analysis
#   write.table(data.frame(mut.index=dataset$removed_indices), file=paste(outfiles.prefix,"_removedMutationsIndex.txt", sep=""), row.names=F, quote=F)
#   
#   # Save the consensus mutation assignments
#   save(file=paste(outfiles.prefix, "_bestConsensusResults.RData", sep=""), output, clustering)
#   colnames(output) = c("chr", "start", "end", "cluster", "likelihood")
#   write.table(output, file=paste(outfiles.prefix, "_bestConsensusAssignments.bed", sep=""), quote=F, row.names=F, sep="\t")
#   
#   # If tree based analysis, also save the tree
#   if (write_tree) {
#     write.table(clustering$best.tree, file=paste(outfiles.prefix, "_bestConsensusTree.txt", sep=""), quote=F, row.names=F, sep="\t")
#   }
# }

#######################################################################################################################
# Binomial based mutation assignment
#######################################################################################################################

#' Assign mutations to clusters by looking at the binomial probability of each cluster for generating a mutation
#' This for now only works with a single timepoint
mutation_assignment_binom = function(clustering_density, mutCount, WTCount, copyNumberAdjustment, tumourCopyNumber, normalCopyNumber, cellularity) {
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
  print(cluster_locations)
  
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
  print(cluster_assignments)
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
#' @return A list containing two fields: localOptima, the location of a peak and peak.indices, the index of the peak within a hypercube
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

#' Obtain the mean density of each cluster. This function takes the cluster_locations and for
#' each cluster it walks from the cluster peak along CCF space until the median density
#' drops below the supplied minimum in both directions. After obtaining the CCF space that a
#' cluster takes up we calculate the mean density in that space. Finally across all cluster
#' densities we normalise to obtain the fraction of total density that each cluster represents
#' @return A vector with for each cluster the fraction of density
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

#######################################################################################################################
# MPEAR based mutation assignment
#######################################################################################################################

#' A function that takes SNVs assigned to clusters and works out a PDF for each cluster, in 512 bins across the given space. These densities can then be used to obtain the probability of each cluster for not assigned SNVs.
#' @param best.node.assignments A vector with hard assignments of SNVs to clusters
#' @param subclonal.fraction CCF values of the SNVs
#' @param min_ccf The minimum CCF value to consider for constructing the density
#' @param max_ccf The minimum CCF value to consider for constructing the density
#' @return A combined density distribution, a column per cluster
get_cluster_densities = function(best.node.assignments, subclonal.fraction, min_ccf=0, max_ccf=1.5) {
  print("Building densities")
  print("Clusters:")
  print(unique(best.node.assignments[!is.na(best.node.assignments)]))
  print("Subclonal fraction:")
  print(head(subclonal.fraction))
  
  #' Get the per cluster densities
  cluster_densities = list()
  for (clusterid in unique(best.node.assignments[!is.na(best.node.assignments)])) {
    ccfs_snvs_cluster = subclonal.fraction[best.node.assignments==clusterid,]
    print("Clusterid:")
    print(clusterid)
    print(head(ccfs_snvs_cluster))
    p = ggplot_build(ggplot(data.frame(x=ccfs_snvs_cluster)) + aes(x=x, y=..density..) + geom_density() + xlim(min_ccf, max_ccf))
    cluster_densities[[length(cluster_densities)+1]] = p$data[[1]]
  }
  
  #' Combining the densities into a matrix
  combined_densities = data.frame(cbind(cluster_densities[[1]][,c("x")],
                                        do.call(cbind, lapply(cluster_densities, function(x) { x[,c("y")] }))))
  
  no.optima = length(unique(best.node.assignments[!is.na(best.node.assignments)]))
  # Rescale the height of the density by the total combined height to obtain fraction of the density per bin
  if (ncol(combined_densities)==2) {
    # Only 1 cluster, so always probability 1
    combined_densities[,(1:no.optima)+1] = 1
  } else {
    combined_densities[,(1:no.optima)+1] = combined_densities[,(1:no.optima)+1] / colSums(combined_densities[,(1:no.optima)+1])
  }
  colnames(combined_densities) = c("x", paste("cluster_", unique(best.node.assignments[!is.na(best.node.assignments)]), sep=""))
  return(combined_densities)
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

#' Function to obtain cluster preferences for individual SNVs via the CCF of assigned clusters, known cluster 
#' locations and a density that approximates the probability of belonging to a certain cluster
#' 
#' This function for now only works in 1D cases
#' @param snv_ccfs
#' @param combined_densities
#' @param no.optima
#' @param no.muts
#' @param no.iters
#' @param no.iters.burn.in
#' @return A list with two tables. One with the fraction of iterations that each SNV would've been assigned to each of the clusters and one where SNVs are not hard assigned, but the probabilities of the clusters are added up per cluster
#' @author sd11
get_preferences_matrix = function(snv_ccfs, combined_densities, no.optima, no.muts, no.iters, no.iters.burn.in) {
  no.iters.post.burnin = no.iters-no.iters.burn.in
  #' Determine the cluster the snv would be assigned to - given the densities
  # clusterids = cluster_locations_table$cluster.no
  # no.optima = length(clusterids)
  mutation.preferences = array(0, c(no.muts, no.optima))
  mutation.sum.probs = array(0, c(no.muts, no.optima))
  for (iter in 1:(no.iters.post.burnin)) {
    assignment_ccf = snv_ccfs[iter, 1:no.muts, 1]
    bins = sapply(1:no.muts, function(i) { which.min(abs(combined_densities$x-assignment_ccf[i])) })
    snv.most.likely.cluster = sapply(bins, function(bin) { which.max(combined_densities[bin,(1:no.optima)+1]) })
    # Get the max probability cluster for this assignment_ccf and count the assignment
    mutation.preferences[cbind(1:no.muts, snv.most.likely.cluster)] = mutation.preferences[cbind(1:no.muts, snv.most.likely.cluster)] + 1
    # Sum the probabilities of each cluster
    mutation.sum.probs = mutation.sum.probs + combined_densities[bins,(1:no.optima)+1]
  }
  mutation.preferences = mutation.preferences / no.iters.post.burnin
  mutation.sum.probs = mutation.sum.probs / no.iters.post.burnin
  return(list(mutation.preferences=mutation.preferences, mutation.sum.probs=mutation.sum.probs))
}

#' Function that works out where SNVs would've been assigned, given an SNV to cluster assignment
get_mutation_preferences_mpear = function(pi.h, S.i, best.node.assignments, subclonal.fraction, no.iters, no.iters.burn.in) {
  no.optima = length(unique(best.node.assignments[!is.na(best.node.assignments)]))
  no.muts = nrow(subclonal.fraction)
  # Obtain densities for each cluster so we know the probability of each cluster vs CCF  
  combined_densities = get_cluster_densities(best.node.assignments, subclonal.fraction, max_ccf=10)
  # Get the CCFs of the clusters to which SNVs have been assigned during MCMC
  snv_ccfs = get_snv_assignment_ccfs(pi.h, S.i, no.muts, ncol(subclonal.fraction), no.iters, no.iters.burn.in)
  # Combine the density with the assignment CCFs to work out where each SNV would've been assigned during MCMC
  prefs = get_preferences_matrix(snv_ccfs, combined_densities, no.optima, no.muts, no.iters, no.iters.burn.in)
  return(list(snv_ccfs=snv_ccfs, mutation.preferences=prefs$mutation.preferences, mutation.sum.probs=prefs$mutation.sum.probs, combined_densities=combined_densities))
}

#' Helper function that builds the a density over assignment CCFs for each mutation
get_snv_ccf_assignmnent_density = function(S.i, pi.h, no.iters.burn.in, ccf_max_value=10) {
  snv_ccfs = DPClust::get_snv_assignment_ccfs(pi.h, S.i, ncol(S.i), dim(pi.h)[3], nrow(S.i), no.iters.burn.in)
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
#' @param S.i
#' @param pi.h
#' @param no.iters.burn.in
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

#' Use the max-PEAR metric to obtain mutation clusters
#' @param GS.data
#' @param no.iters
#' @param no.iters.burn.in
#' @param min.frac.snvs.cluster
#' @return 
#' @author sd11
mutation_assignment_mpear = function(GS.data, no.iters, no.iters.burn.in, min.frac.snvs.cluster, dataset, samplename, outdir, density) {
  
  #' Take a list of assigned labels and a dataset and create the cluster locations table by 
  #' using the mean CCF of the assigned mutations per cluster
  make_clusters = function(filtered_label_assignments, dataset) {
    # Build a DPClust clusters table
    cluster_locs = data.frame()
    for (i in unique(filtered_label_assignments)) {
      if (!is.na(i)) {
        assigned_muts = filtered_label_assignments==i
        assigned_muts[is.na(assigned_muts)] = F
        cluster_locs = rbind(cluster_locs, data.frame(cluster.no=i, 
                                                      location=mean(dataset$subclonal.fraction[assigned_muts,]), 
                                                      no.mutations=sum(assigned_muts)))
      }
    }
    return(cluster_locs)
  }
  
  #' Set assigned labels to small clusters to NA
  remove_small_clusters = function(label_assignments, min.frac.snvs.cluster) {
    num.muts = length(label_assignments)
    cluster_sizes = table(label_assignments)
    filtered_label_assignments = label_assignments
    for (clusterid in unique(label_assignments[!is.na(label_assignments)])) {
      if (cluster_sizes[as.character(clusterid)] < floor(min.frac.snvs.cluster*num.muts) & cluster_sizes[as.character(clusterid)] < 30) {
        filtered_label_assignments[label_assignments==clusterid] = NA
      }
    }
    return(filtered_label_assignments)
  }
  
  #' Convenience function to make the plot
  make_figure = function(dataset, mutation_assignments, filename, max_ccf=1.5) {
    dat_ccf_mpear = cbind(data.frame(ccf=dataset$subclonal.fraction[,1]), data.frame(cluster=factor(mutation_assignments)))
    p = ggplot(dat_ccf_mpear) + aes(x=ccf, y=..count.., fill=cluster) + geom_bar(binwidth=0.05, colour="black") + xlim(0, max_ccf)
    png(filename, width=700, height=400)
    print(p)
    dev.off()
  }
  
  
  num.muts = ncol(GS.data$S.i)
  # Obtain coassignment posterior estimates from the trace
  # coassignments = mcclust::comp.psm(GS.data$S.i[no.iters.burn.in:no.iters, ])
  # coassignments = build_coassignment_prob_matrix_densities(GS.data$S.i, GS.data$pi.h, no.iters.burn.in)
  coassignments = build_coassignment_prob_matrix_preferences(GS.data, density, num.muts, no.iters.burn.in:no.iters)
  # Get the max PEAR
  mpear = mcclust::maxpear(coassignments)
  label_assignments = mpear$cl
  
  # Save the computationally heavy results
  save(file=paste(outdir, samplename, "_coassignment_matrix.RData", sep=""), coassignments, mpear)
  
  # Remove too small clusters
  filtered_label_assignments = remove_small_clusters(label_assignments, min.frac.snvs.cluster)
  
  # Build a DPClust clusters table
  cluster_locs = make_clusters(filtered_label_assignments, dataset)
  make_figure(dataset=dataset, 
              mutation_assignments=filtered_label_assignments, 
              filename=paste(outdir, samplename, "_mpear_large_clusters.png", sep="")) 
  
  # Build a mutation preferences table
  mutation.preferences = get_mutation_preferences_mpear(GS.data$pi.h, GS.data$S.i, filtered_label_assignments, dataset$subclonal.fraction, no.iters, no.iters.burn.in)
  save(file=paste(outdir, samplename, "_mutation_preferences.RData", sep=""), mutation.preferences)
  preferences = mutation.preferences$mutation.preferences
  
  # Assign the unassigned mutations
  unassigned_index = which(is.na(filtered_label_assignments))
  if (length(unassigned_index) > 0) {
    most_likely_cluster = sapply(unassigned_index, function(i) { which.max(preferences[i,]) })
    # Save the id of the most likely cluster - sort because the preferences table is sorted by cluster id, which may not be the same ordering as the clusterin
    filtered_label_assignments_final = filtered_label_assignments
    filtered_label_assignments_final[is.na(filtered_label_assignments_final)] = sort(cluster_locs$cluster.no)[most_likely_cluster]
  } else {
    filtered_label_assignments_final = filtered_label_assignments
  }
  
  # TODO: maybe update the cluster locations?
  
  make_figure(dataset=dataset, 
              mutation_assignments=filtered_label_assignments_final, 
              filename=paste(outdir, samplename, "_mpear_large_clusters_allsnvs.png", sep=""))
  
  
  # # Get probabilities for each mutation and cluster
  # assignment_probs = calcBinomProbs(dataset, cluster_locs$location)
  # 
  # # Assign the remaining mutations to the most likely cluster
  # most.likely.cluster = cluster_locs$cluster.no[apply(assignment_probs, 1, which.max)]
  # filtered_label_assignments[is.na(filtered_label_assignments)] = most.likely.cluster[is.na(filtered_label_assignments)]
  # 
  # Rebuild the clusters with all mutations assigned
  cluster_locs = make_clusters(filtered_label_assignments_final, dataset)
  
  # Rename the clusters to have the ids line up
  filtered_label_assignments_final.temp = filtered_label_assignments_final
  for (i in nrow(cluster_locs):1) {
    filtered_label_assignments_final[filtered_label_assignments_final.temp==cluster_locs$cluster.no[i]] = i
  }
  cluster_locs$cluster.no = nrow(cluster_locs):1
  
  # Obtain the probability of the most likely cluster to return to upstream
  best.assignment.likelihoods = preferences[cbind(1:length(filtered_label_assignments_final), filtered_label_assignments_final)]
  return(list(best.node.assignments=filtered_label_assignments_final, best.assignment.likelihoods=best.assignment.likelihoods, cluster.locations=cluster_locs, all.assignment.likelihoods=preferences))
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

# #' Calculate for a pair of clusters whether a has a higher CCF than b.
# #' 
# #' This function uses mutation assignments during MCMC of the mutations that have been assigned
# #' to clusters a and b after completion of the run. It samples 1000 mutations from a and b and
# #' looks how often a_i has a higher assignment CCF than b_i and combines this information in a
# #' fraction. The final table contains for each cell a probability whether the column has a higher
# #' CCF than the row.
# #' @param GS.data MCMC output with assignments and cluster locations
# #' @param clusterids Clusterids to run through
# #' @param no.muts The total number of mutations
# #' @param no.timepoints Total number of samples in this dataset
# #' @param no.iters Total number of iterations
# #' @param no.iters.burn.in Number of iterations to use as burn-in
# #' @return A array multi-dimensional array with in each cell whether the column cluster has a higher CCF than the row cluster across the samples in the third dimension
# #' @author sd11
# calc_cluster_order_probs = function(GS.data, density, mut_assignments, clusterids, cluster_ccfs, no.muts, no.timepoints, no.iters, no.iters.burn.in, no.samples=1000) {
#   num_clusters = length(clusterids)
#   if (num_clusters > 1) {
#     preferences = get_mutation_preferences(GS.data, density, mut_assignments, clusterids, cluster_ccfs, no.muts, no.timepoints, no.iters, no.iters.burn.in)
# 
#     probs = array(NA, c(length(clusterids), length(clusterids), no.timepoints))
#     for (t in 1:no.timepoints) {
#       for (c in 1:(num_clusters-1)) {
#         for (k in (c+1):num_clusters) {
# 
#           snvs_a = which(mut_assignments==clusterids[c])
#           snvs_b = which(mut_assignments==clusterids[k])
# 
#           sampled_a = sample(snvs_a, no.samples, replace=T)
#           sampled_b = sample(snvs_b, no.samples, replace=T)
# 
#           gt = sum(sapply(1:no.samples, function(i) { (sum(preferences[, sampled_a[i], t] >= preferences[, sampled_b[i], t])) }))
#           lt = sum(sapply(1:no.samples, function(i) { (sum(preferences[, sampled_a[i], t] <= preferences[, sampled_b[i], t])) }))
# 
#           probs[c, k, t] = lt / (lt+gt)
#           probs[k, c, t] = gt / (lt+gt)
#         }
#       }
#     }
#     return(probs)
#   } else {
#     return(array(NA, c(length(clusterids), length(clusterids), no.timepoints)))
#   }
# }

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
