RunTreeBasedDPConsensus <- function(trees, node.assignments, mutCount, WTCount, kappa, samplename, subsamplenames, annotation, no.iters, no.iters.burn.in, resort.mutations, bin.indices=NULL) {

  no.subsamples = ncol(mutCount)
  temp.list = GetConsensusTrees(trees, node.assignments, mutCount, WTCount, subclonal.fraction=NULL, kappa = kappa, samplename = samplename, subsamplenames = subsamplenames, no.iters = no.iters, no.iters.burn.in = no.iters.burn.in, resort.mutations = resort.mutations)
  all.consensus.trees = temp.list$all.consensus.trees
  save(all.consensus.trees,file=paste(samplename,"_",no.iters,"iters_",no.iters.burn.in,"_allConsensusTrees.RData",sep=""))

  print(names(temp.list))
  
  if(!is.null(bin.indices)) {
    all.consensus.assignments = temp.list$all.consensus.assignments
    all.disaggregated.consensus.assignments = temp.list$all.disaggregated.consensus.assignments
    save(all.consensus.assignments,file=paste(samplename,"_",no.iters,"iters_",no.iters.burn.in,"_allBinnedConsensusAssignments.RData",sep=""))
    save(all.disaggregated.consensus.assignments,file=paste(samplename,"_",no.iters,"iters_",no.iters.burn.in,"_allConsensusAssignments.RData",sep=""))		
  }else{
    all.consensus.assignments = temp.list$all.consensus.assignments
    save(all.consensus.assignments,file=paste(samplename,"_",no.iters,"iters_",no.iters.burn.in,"_allConsensusAssignments.RData",sep=""))
  }
  
  likelihoods = temp.list$likelihoods; BIC = temp.list$BIC; AIC = temp.list$AIC; DIC = temp.list$DIC
  
  plotScores("log_likelihood", samplename, no,iters, 0, likelihoods, "log-likelihood")
  plotScores("BIC", samplename, no,iters, 0, BIC, "BIC")
  plotScores("AIC", samplename, no,iters, 0, AIC, "AIC")
  plotScores("DIC", samplename, no,iters, 0, DIC, "DIC")
  
  writeScoresTable(likelihoods, "likelihood", samplename, no,iters)
  writeScoresTable(BIC, "BIC", samplename, no,iters)
  writeScoresTable(AIC, "AIC", samplename, no,iters)
  writeScoresTable(DIC, "DIC", samplename, no,iters)
  
  best.index = which.min(BIC)
  best.tree = all.consensus.trees[[best.index]]

  if(!is.null(bin.indices)) {
    best.node.assignments = all.disaggregated.consensus.assignments[[best.index]]
  }else{
    best.node.assignments = all.consensus.assignments[[best.index]]
  }
  
  png(paste("bestConsensusTree_",samplename,"_",no.iters,"iters.png",sep=""),width=no.subsamples*1000,height=1000)
  if(!is.null(bin.indices)) {
    plotTree(annotateTree(best.tree,all.disaggregated.consensus.assignments,annotation),font.size=2,main=paste(samplename,subsamplenames,sep=""))
  }else{
    plotTree(annotateTree(best.tree,all.consensus.assignments,annotation),font.size=2,main=paste(samplename,subsamplenames,sep=""))
  }
  dev.off()
  write.table(best.node.assignments,paste(samplename,"_",no.iters,"iters_",no.iters.burn.in,"_bestConsensusAssignments.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
  
  if(no.subsamples>1){
    plotBestScatter(mutCount, WTCount, kappa, best.tree, samplename, no.iters, subsamplenames, best.node.assignments)
  }	
}

writeScoresTable <- function(scores, score_name, samplename, no,iters) {
  #   write.table(data.frame(tree.index = 1:length(likelihoods),likelihood=likelihoods),paste("consensus_likelihoods_",samplename,"_",no.iters,"iters.txt",sep=""),sep="\t",row.names=F,quote=F)
  data = data.frame(tree.index = 1:length(scores),likelihood=scores)
  colnames(data) = c("tree.index", score_name)
  write.table(data,paste("consensus_",score_name,"_",samplename,"_",no.iters,"iters.txt",sep=""),sep="\t",row.names=F,quote=F)
}

plotBestScatter <- function(mutCount, WTCount, kappa, best.tree, samplename, no.iters, subsamplenames, best.node.assignments) {
  no.subsamples = ncol(mutCount)
  subclonal.fraction = array(NA,dim(mutCount))
  for(i in 1:no.subsamples){
    subclonal.fraction[,i] = mutCount[,i] / ((mutCount[,i]+WTCount[,i])*kappa[,i])
    subclonal.fraction[kappa[,i]==0,i] = NA
  }  
  
  #its hard to distinguish more than 8 different colours
  max.cols = 8
  no.nodes = nrow(best.tree)
  unique.nodes = rownames(best.tree)
  cols = rainbow(min(max.cols,no.nodes))
  png(paste("bestScatterPlots_",samplename,"_",no.iters,"iters.png",sep=""),width=no.subsamples*(no.subsamples-1)*500,height=1000)
  par(mfrow = c(1,no.subsamples*(no.subsamples-1)/2),cex=2)
  plot.data = subclonal.fraction
  plot.data[is.na(plot.data)]=0
  for(i in 1:(no.subsamples-1)){
    for(j in (i+1):no.subsamples){
      plot(plot.data[,i],plot.data[,j],type = "n",xlab = paste(samplename,subsamplenames[i]," subclonal fraction",sep=""), ylab = paste(samplename,subsamplenames[j]," subclonal fraction",sep=""),xlim = c(0,max(plot.data[,i])*1.3))
      for(n in 1:no.nodes){
        points(plot.data[,i][best.node.assignments==unique.nodes[n]],plot.data[,j][best.node.assignments==unique.nodes[n]],col=cols[(n-1) %% max.cols + 1],pch=20 + floor((n-1)/max.cols))
      }
      #legend(max(plot.data[,i]),max(plot.data[,j]),legend = unique.nodes[n],col=cols[(n-1) %% max.cols + 1],pch=20 + floor((n-1)/max.cols))
      legend(max(plot.data[,i])*1.1,max(plot.data[,j]),legend = unique.nodes,col=cols[(0:(no.nodes-1)) %% max.cols + 1],pch=20 + floor((0:(no.nodes-1))/max.cols))
    }
  }
  dev.off()
}
