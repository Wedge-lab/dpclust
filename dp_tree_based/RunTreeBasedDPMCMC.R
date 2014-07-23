RunTreeBasedDPMCMC <- function(mutCount, WTCount, kappa, no.muts.input, annotation, no.iters, shrinkage.threshold, init.alpha, outdir, parallel, clp, blockid=1, bin.indices=NULL) {

#   if(is.na(bin.size)){
  temp.list = tree.struct.dirichlet.gibbs(y=mutCount,n=WTCount+mutCount,kappa=kappa,iter=no.iters,shrinkage.threshold=shrinkage.threshold,init.alpha=init.alpha, remove.node.frequency = NA, remove.branch.frequency = NA, parallel=parallel, cluster=clp)
#   }else{
#     save(bin.indices,file=paste(outdir,"/",samplename,"_",no.iters,"iters_burnin",no.iters.burn.in,"_binnedIndices.RData",sep=""))
#     temp.list = tree.struct.dirichlet.gibbs(y=binned.mutCount,n=binned.WTCount+binned.mutCount,kappa=binned.kappa,iter=no.iters,shrinkage.threshold=shrinkage.threshold,init.alpha=init.alpha, remove.node.frequency = NA, remove.branch.frequency = NA, parallel=parallel, cluster=clp)	
#   }
  
  setwd(outdir)
  
  trees = temp.list$trees #[[1]]
#   if(is.na(bin.size)){
  if(is.null(bin.indices)) {
    node.assignments = temp.list$assignments #[[2]]
  }else{
    binned.node.assignments = temp.list$assignments #[[2]]
    node.assignments = array(NA,c(no.muts.input,ncol(binned.node.assignments)))
    for(i in 1:length(bin.indices)){
#       print("new")
#       print(i)
#       print(bin.indices[[i]])
#       print(dim(node.assignments))
#       print(node.assignments[bin.indices[[i]],])
#       print(binned.node.assignments[i,])
      node.assignments[bin.indices[[i]],] = binned.node.assignments[i,]
      
      
#       Error in print(node.assignments[bin.indices[[i]], ]) : 
#         error in evaluating the argument 'x' in selecting a method for function 'print': Error in node.assignments[bin.indices[[i]], ] : subscript out of bounds
    }
    write.table(binned.node.assignments,paste("aggregated_node_assignments_",samplename,"_",no.iters,"iters_block",blockid,".txt",sep=""),sep="\t",row.names=F,quote=F)
  }
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
  
  pdf(paste("all_trees_",samplename,"_",no.iters,"iters_block",blockid,".pdf",sep=""),height=4,width=3*no.subsamples)
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
  
  #   	png(paste("log_likelihood_plot_",samplename,"_",no.iters,"iters.png",sep=""),width=1000)
  #   	par(cex.lab=3,cex.axis=3,lwd=3,mar=c(7,7,5,2))
  #   	plot(1:no.iters,likelihoods,type="l",col="red",xlab="MCMC iteration", ylab="log-likelihood",main=samplename)
  #   	dev.off()
  #   	png(paste("BIC_plot_",samplename,"_",no.iters,"iters.png",sep=""),width=1000)
  #   	par(cex.lab=3,cex.axis=3,lwd=3,mar=c(7,7,5,2))
  #   	plot(1:no.iters,BIC,type="l",col="red",xlab="MCMC iteration", ylab="BIC",main=samplename)
  #   	dev.off()
  #   	png(paste("AIC_plot_",samplename,"_",no.iters,"iters.png",sep=""),width=1000)
  #   	par(cex.lab=3,cex.axis=3,lwd=3,mar=c(7,7,5,2))
  #   	plot(1:no.iters,AIC,type="l",col="red",xlab="MCMC iteration", ylab="AIC",main=samplename)
  #   	dev.off()
  #     
  #     
  #   	png(paste("DIC_plot_",samplename,"_",no.iters,"iters.png",sep=""),width=1000)
  #   	par(cex.lab=3,cex.axis=3,lwd=3,mar=c(7,7,5,2))
  #   	plot(1:no.iters,DIC,type="l",col="red",xlab="MCMC iteration", ylab="DIC",main=samplename)
  #   	dev.off()
  
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
  #write.table(best.tree,paste("best_tree_",samplename,"_",no.iters,"iters_block",blockid,".txt",sep=""),sep="\t",row.names=F,quote=F)
  save(trees,file=paste("",samplename,"_trees_iters",no.iters,"_block",blockid,".Rdata",sep=""))
  
}

plotScores <- function(file_prefix, samplename, no,iters, blockid, scores, ylab, xlab="posterior tree index",consensus=TRUE) {
  #   png(paste("log_likelihood_plot_",samplename,"_consensus_",no.iters,"iters.png",sep=""),width=1000)
  if (consensus) {
    filename = paste(file_prefix,"_plot_",samplename,"_consensus_",no.iters,"iters_block",blockid,".png",sep="")
  } else {
    filename = paste(file_prefix,"_plot_",samplename,"_",no.iters,"iters.png",sep="")
  }
  png(filename,width=1000)
  par(cex.lab=3,cex.axis=3,lwd=3,mar=c(7,7,5,2))
  plot(1:length(scores),scores,type="l",col="red",xlab=xlab, ylab=ylab,main=samplename)
  dev.off()
}