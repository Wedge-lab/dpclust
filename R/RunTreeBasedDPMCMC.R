RunTreeBasedDPMCMC <- function(mutCount, WTCount, kappa, no.muts.input, annotation, samplename, no.iters, no.iters.burn.in, shrinkage.threshold, init.alpha, outdir, parallel, clp, blockid=1, bin.indices=NULL, remove.node.frequency=NA, remove.branch.frequency=NA, conflict.array=NA) {

  setwd(outdir)
  checkpoint_file = paste("MCMC_output_block", blockid, ".RData", sep="")
  if (!file.exists(checkpoint_file)) {
    # Run the Gibbs sampler
    temp.list = tree.struct.dirichlet.gibbs(y=mutCount,
                                            n=WTCount+mutCount,
                                            kappa=kappa,
                                            iter=no.iters,
                                            shrinkage.threshold=shrinkage.threshold,
                                            init.alpha=init.alpha, 
                                            remove.node.frequency=remove.node.frequency, 
                                            remove.branch.frequency=remove.branch.frequency, 
                                            parallel=parallel, 
                                            cluster=clp,
                                            conflict.array=conflict.array)
  
    # Save the output
    save(checkpoint_file, temp.list)
  } else {
    # Load the checkpoint
    load(checkpoint_file)
  }
  
  trees = temp.list$trees
  if(is.null(bin.indices)) {
    node.assignments = temp.list$assignments
    binned.node.assignments = NA
    
    print("Calculating mutation strengths")
    strengths = get.mut.ass.strengths(nrow(node.assignments), no.iters, no.iters-no.iters.burn.in, node.assignments)

  }else{
    # Unpack the bins
    binned.node.assignments = temp.list$assignments #[[2]]
    node.assignments = array(NA,c(no.muts.input,ncol(binned.node.assignments)))
    for(i in 1:length(bin.indices)) {
      node.assignments[bin.indices[[i]],] = binned.node.assignments[i,]
    }
    write.table(binned.node.assignments,paste("aggregated_node_assignments_",samplename,"_",no.iters,"iters_block",blockid,".txt",sep=""),sep="\t",row.names=F,quote=F)

    print("Calculating mutation strengths")
    strengths = get.mut.ass.strengths(nrow(binned.node.assignments), no.iters, no.iters-no.iters.burn.in, binned.node.assignments)
  }

  ancestor.strengths = strengths$ancestor.strengths
  sibling.strengths = strengths$sibling.strengths
  identity.strengths = strengths$identity.strengths
  parent.child.strengths = strengths$parent.child.strengths
  child.parent.strengths = strengths$child.parent.strengths

  # Save the strengths to disk
  write.table(ancestor.strengths, paste("ancestor.strengths_",samplename,"_",no.iters,"iters_block",blockid,".txt",sep=""),sep="\t",row.names=F,quote=F,col.names=F)
  write.table(sibling.strengths, paste("sibling.strengths_",samplename,"_",no.iters,"iters_block",blockid,".txt",sep=""),sep="\t",row.names=F,quote=F,col.names=F)
  write.table(identity.strengths, paste("identity.strengths_",samplename,"_",no.iters,"iters_block",blockid,".txt",sep=""),sep="\t",row.names=F,quote=F,col.names=F)
  write.table(parent.child.strengths, paste("parent.child.strengths_",samplename,"_",no.iters,"iters_block",blockid,".txt",sep=""),sep="\t",row.names=F,quote=F,col.names=F)
  write.table(child.parent.strengths, paste("child.parent.strengths_",samplename,"_",no.iters,"iters_block",blockid,".txt",sep=""),sep="\t",row.names=F,quote=F,col.names=F)

  alphas = temp.list$alpha
  lambdas = temp.list$lambda
  gammas = temp.list$gamma
  likelihoods = temp.list$likelihood
  BIC = temp.list$BIC
  AIC = temp.list$AIC
  DIC = temp.list$DIC
  mant.anc = temp.list$mant.anc
  mant.sib = temp.list$mant.sib
  mant.ide = temp.list$mant.ide
  
  tree.sizes = unlist(lapply(trees, nrow))
  
  best.index = which.max(likelihoods)
  best.BIC.index = which.min(BIC)
  best.AIC.index = which.min(AIC)
  best.DIC.index = which.min(DIC)
  
  pdf(paste("all_trees_",samplename,"_",no.iters,"iters_block",blockid,".pdf",sep=""),height=4,width=3*ncol(mutCount))
  for(iter in 1:no.iters){
    tree = trees[[iter]]
    tree$annotation = NA
    tree = annotateTree(tree,node.assignments[,iter],annotation)
    title = paste(samplename,"iter",iter)
    if(iter==best.index){
      title = paste(title,"(most likely tree)")
    }
    if(iter==best.BIC.index){
      title = paste(title,"(best BIC tree)")
    }
    if(iter==best.AIC.index) {
      title = paste(title,"(best AIC tree)")
    }
    if(iter==best.DIC.index) {
      title = paste(title,"(best DIC tree)")
    }
    plotTree(tree,title)
  }
  dev.off()
  

  plotScores("log_likelihood", samplename, no,iters, blockid, likelihoods, "log-likelihood", xlab="MCMC iteration",consensus=FALSE)
  plotScores("BIC", samplename, no,iters, blockid, BIC, "BIC", xlab="MCMC iteration",consensus=FALSE)
  plotScores("AIC", samplename, no,iters, blockid, AIC, "AIC", xlab="MCMC iteration",consensus=FALSE)
  plotScores("DIC", samplename, no,iters, blockid, DIC, "DIC", xlab="MCMC iteration",consensus=FALSE)
  plotScores("mantel_ancestors", samplename, no,iters, blockid, mant.anc, "Mantel correlation", xlab="MCMC iteration / 10",consensus=FALSE)
  plotScores("mantel_siblings", samplename, no,iters, blockid, mant.sib, "Mantel correlation", xlab="MCMC iteration / 10",consensus=FALSE)
  plotScores("mantel_identity", samplename, no,iters, blockid, mant.ide, "Mantel correlation", xlab="MCMC iteration / 10",consensus=FALSE)
  
  write.table(node.assignments,paste("node_assignments_",samplename,"_",no.iters,"iters_block",blockid,".txt",sep=""),sep="\t",row.names=F,quote=F)
  write.table(alphas,paste("alphas_",samplename,"_",no.iters,"iters_block",blockid,".txt",sep=""),sep="\t",row.names=F,quote=F)
  write.table(lambdas,paste("lambdas_",samplename,"_",no.iters,"iters_block",blockid,".txt",sep=""),sep="\t",row.names=F,quote=F)
  write.table(gammas,paste("gammas_",samplename,"_",no.iters,"iters_block",blockid,".txt",sep=""),sep="\t",row.names=F,quote=F)
  write.table(likelihoods,paste("likelihoods_",samplename,"_",no.iters,"iters_block",blockid,".txt",sep=""),sep="\t",row.names=F,quote=F)
  save(trees,file=paste("",samplename,"_trees_iters",no.iters,"_block",blockid,".Rdata",sep=""))
  return(list(trees=trees, 
              node.assignments=node.assignments, 
              alphas=alphas, 
              lambdas=lambdas,
              gammas=gammas, 
              likelihoods=likelihoods, 
              binned.node.assignments=binned.node.assignments, 
              strengths=strengths))
}

plotScores <- function(file_prefix, samplename, no,iters, blockid, scores, ylab, xlab="posterior tree index",consensus=TRUE) {
  #   png(paste("log_likelihood_plot_",samplename,"_consensus_",no.iters,"iters.png",sep=""),width=1000)
  if (consensus) {
    filename = paste(file_prefix,"_plot_",samplename,"_consensus_",no.iters,"iters_block",blockid,".png",sep="")
  } else {
    filename = paste(file_prefix,"_plot_",samplename,"_",no.iters,"iters_block",blockid,".png",sep="")
  }
  png(filename,width=1000)
  par(cex.lab=3,cex.axis=3,lwd=3,mar=c(7,7,5,2))
  plot(1:length(scores),scores,type="l",col="red",xlab=xlab, ylab=ylab,main=samplename)
  dev.off()
}

get.mut.ass.strengths = function(no.muts, no.iters, no.iters.post.burn.in, node.assignments) {
	#
	# Calculates for each pair of mutations how often the pair is assigned to:
	#   - The same node (identity)
	#   - Nodes in parent-offspring relation (ancestor)
	#   - Nodes that are siblings (sibling)
	#
	calc.ancestor.strengths <- function(m,ancestor.strengths, node.assignments, no.iters.since.burnin, identity.strengths, parent.child.strengths, child.parent.strengths) {
		#
		# Calculates the strengths for iteration m of the MCMC algorithm
		#
		node.assignments.all = node.assignments[,no.iters.since.burnin]
		node.assignments.m = node.assignments[m,no.iters.since.burnin]
		#node.assignments.m = node.assignments.all[m]

		temp.ancestor.or.identity.relationship = younger.direct.descendants(node.assignments.m,node.assignments.all)
		temp.ancestor.strengths = ancestor.strengths[m,] + (temp.ancestor.or.identity.relationship & (node.assignments.m != node.assignments.all))
		temp.identity.strengths = identity.strengths[m,] + (node.assignments.m == node.assignments.all)
		temp.parent.child.strengths = parent.child.strengths[m,] + ancestors.only(node.assignments.m,node.assignments.all)
		temp.child.parent.strengths = child.parent.strengths[m,] + younger.descendants.only(node.assignments.m,node.assignments.all)

		return(list(temp.ancestor.strengths, temp.ancestor.or.identity.relationship, temp.identity.strengths, temp.parent.child.strengths, temp.child.parent.strengths))
	}

	ancestor.strengths = array(0,c(no.muts,no.muts))
	sibling.strengths = array(0,c(no.muts,no.muts))
	identity.strengths = array(0,c(no.muts,no.muts))
	parent.child.strengths = array(0,c(no.muts,no.muts))
	child.parent.strengths = array(0,c(no.muts,no.muts))

	for(i in 1:no.iters.post.burn.in){ # for each iteration past burnin
		if (i %% 100 == 0) { print(i) }

		ancestor.or.identity.relationship = array(NA,c(no.muts,no.muts))

		res = sapply(1:no.muts, FUN=calc.ancestor.strengths, ancestor.strengths, node.assignments, i+no.iters-no.iters.post.burn.in, identity.strengths,  parent.child.strengths, child.parent.strengths)
		# res contains three lists of vectors. Below we take a list and rbind all vectors in it together into a matrix
		ancestor.strengths = do.call(rbind,res[1,])
		ancestor.or.identity.relationship = do.call(rbind,res[2,])
		identity.strengths = do.call(rbind,res[3,])
		parent.child.strengths = do.call(rbind,res[4,])
		child.parent.strengths = do.call(rbind,res[5,])

		sibling.strengths = sibling.strengths + as.numeric(!ancestor.or.identity.relationship & !t(ancestor.or.identity.relationship))
	}

	return(list(ancestor.strengths=ancestor.strengths, sibling.strengths=sibling.strengths, identity.strengths=identity.strengths, parent.child.strengths=parent.child.strengths, child.parent.strengths=child.parent.strengths))
}

# get.mut.ass.strengths = function(no.muts, no.iters, no.iters.post.burn.in, node.assignments) {
#   ancestor.strengths = array(0,c(no.muts,no.muts))
#   sibling.strengths = array(0,c(no.muts,no.muts))
#   identity.strengths = array(0,c(no.muts,no.muts))
#   parent.child.strengths = array(0,c(no.muts,no.muts))
#   child.parent.strengths = array(0,c(no.muts,no.muts))
#   #it would be faster to use apply
#   for(i in 1:no.iters.post.burn.in){
#     ancestor.or.identity.relationship = array(NA,c(no.muts,no.muts))
#     for(m in 1:no.muts){
#       ancestor.strengths[m,] = ancestor.strengths[m,] + (younger.direct.descendants(node.assignments[m,i+no.iters-no.iters.post.burn.in], node.assignments[,i+no.iters-no.iters.post.burn.in]) & (node.assignments[m,i+no.iters-no.iters.post.burn.in] != node.assignments[,i+no.iters-no.iters.post.burn.in]))
#       ancestor.or.identity.relationship[m,] = younger.direct.descendants(node.assignments[m,i+no.iters-no.iters.post.burn.in], node.assignments[,i+no.iters-no.iters.post.burn.in])
#       identity.strengths[m,] = identity.strengths[m,] + (node.assignments[m,i+no.iters-no.iters.post.burn.in] == node.assignments[,i+no.iters-no.iters.post.burn.in])
#       parent.child.strengths[m,] = parent.child.strengths[m,] + ancestors.only(node.assignments[m,i+no.iters-no.iters.post.burn.in], node.assignments[,i+no.iters-no.iters.post.burn.in])
#       child.parent.strengths[m,] = child.parent.strengths[m,] + younger.descendants.only(node.assignments[m,i+no.iters-no.iters.post.burn.in], node.assignments[,i+no.iters-no.iters.post.burn.in])
#     }
#     sibling.strengths = sibling.strengths + as.numeric(!ancestor.or.identity.relationship & !t(ancestor.or.identity.relationship))
#   }
#   
#   return(list(ancestor.strengths=ancestor.strengths, sibling.strengths=sibling.strengths, identity.strengths=identity.strengths, parent.child.strengths=parent.child.strengths, child.parent.strengths=child.parent.strengths))
# }

# get.mut.ass.strengths = function(no.muts, no.iters, no.iters.post.burn.in, node.assignments) {
#   #
#   # Calculates for each pair of mutations how often the pair is assigned to:
#   #   - The same node (identity)
#   #   - Nodes in parent-offspring relation (ancestor)
#   #   - Nodes that are siblings (sibling)
#   #
#   calc.ancestor.strengths <- function(m, node.assignments, no.iters.since.burnin) {
#     #
#     # Calculates the strengths for iteration m of the MCMC algorithm
#     #
#     node.assignments.all = node.assignments[,no.iters.since.burnin]
#     node.assignments.m = node.assignments[m,no.iters.since.burnin]
#     
#     temp.ancestor.or.identity.relationship = younger.direct.descendants(node.assignments.m,node.assignments.all)
#     temp.ancestor.strengths = (temp.ancestor.or.identity.relationship & (node.assignments.m != node.assignments.all))
#     temp.identity.strengths = (node.assignments.m == node.assignments.all)
#     temp.parent.child.strengths = ancestors.only(node.assignments.m,node.assignments.all)
#     temp.child.parent.strengths = younger.descendants.only(node.assignments.m,node.assignments.all)
#     
#     return(list(temp.ancestor.strengths, temp.ancestor.or.identity.relationship, temp.identity.strengths, temp.parent.child.strengths, temp.child.parent.strengths))
#   }
#   
#   # Combine function that takes a list of lists and returns a new list of lists where the first list contains each first item, second each second item, etc
#   comb <- function(x, ...) {
#     lapply(seq_along(x), function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
#   }
#   
#   chunkSize = ceiling(no.iters.post.burn.in/4)
#   r = foreach(i=1:4, .packages="doRNG", .export=c("comb", "calc.ancestor.strengths","younger.direct.descendants","ancestors.only","younger.descendants.only"), .combine='comb', .multicombine=T, .init=list(list(), list(), list(), list(), list())) %dorng% { 
#     # Set start and end points of this chunk
#     s = ((i-1)*chunkSize)+1
#     e = i*chunkSize
#     # Iterate over the items in the range. res will be a list of lists
#     ri = foreach(j=s:e, .combine='comb', .multicombine=T, .init=list(list(), list(), list(), list(), list())) %do% {	     
#       if (j <= no.iters.post.burn.in) {
#         # Calculate the 
#         res = sapply(1:no.muts, FUN=calc.ancestor.strengths, node.assignments, j+no.iters-no.iters.post.burn.in)
#         ancestor.or.identity.relationship = do.call(rbind,res[2,])*1
#         list(do.call(rbind,res[1,])*1, do.call(rbind,res[3,])*1, do.call(rbind,res[4,])*1, do.call(rbind,res[5,])*1, (!ancestor.or.identity.relationship & !t(ancestor.or.identity.relationship))*1)
#       } else {
#         # This case must be caught to ensure all chunks return a list of the equal size
#         temp = array(0,c(no.muts,no.muts))
#         list(temp, temp, temp, temp, temp)
#       }
#     }
#     # Partially add the matrices for this chunk per category
#     lapply(ri, function(x) { Reduce('+', x) })
#   }
#   
#   # Add the summed matrices per category, per chunk and prepare the output
#   res = lapply(r, function(x) { Reduce('+', x) })
#   ancestor.strengths = res[[1]]
#   identity.strengths = res[[2]]
#   parent.child.strengths = res[[3]]
#   child.parent.strengths = res[[4]]
#   sibling.strengths = res[[5]]
#   
#   return(list(ancestor.strengths=ancestor.strengths, sibling.strengths=sibling.strengths, identity.strengths=identity.strengths, parent.child.strengths=parent.child.strengths, child.parent.strengths=child.parent.strengths))
# }
