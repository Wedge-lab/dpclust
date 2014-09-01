source("DensityEstimator.R")
source("PlotTreeWithIgraph.R")

makeValidThetas<-function(tree){
  theta.cols = names(tree)[grep("theta",names(tree))]
  levels = sapply(1:nrow(tree),function(t,i){length(strsplit(rownames(t)[i],":")[[1]])},t=tree)

  if(max(levels)>1){
    for(level in 2:max(levels)){
      level.list = rownames(tree)[levels==level]
      for(anc in unique(tree[level.list,"ancestor"])){
        node.list = level.list[tree$ancestor[match(level.list,rownames(tree))]==anc]
        for(tc in theta.cols){
          if(sum(tree[node.list,tc]) > 0.99 * tree[anc,tc]){
            tree[node.list,tc] = tree[node.list,tc] * 0.99 * tree[anc,tc] / sum(tree[node.list,tc])
          }
        }
      }
    }
  }
  return(tree)
}

#plotConsensusTree<-function(consensus.assignments,samplename,subsamplenames,no.iters, no.iters.burn.in, node.assignments,trees,mutCount,WTCount,kappa,hist.device,density.device,tree.device,optimised.tree.device,tree.population.device,scatter.device,bin.indices=NULL,shrinkage.threshold = 0.1, tree.number=NA){
plotConsensusTree<-function(consensus.assignments, samplename, subsamplenames, no.iters, no.iters.burn.in, node.assignments, trees, mutCount, WTCount, kappa, plot.devs, bin.indices=NULL, shrinkage.threshold=0.1, tree.number=NA){
  start = Sys.time()
  print("Starting plotConcensusTree")
  no.iters.post.burn.in = no.iters-no.iters.burn.in
  no.subsamples = length(subsamplenames)
  no.muts = length(consensus.assignments)
  
  subclonal.fraction = mutCount / ((mutCount+WTCount)*kappa)
  subclonal.fraction[kappa == 0] = NA
  
  print("Getting node info")
  ############################## Get Node Info ############################################################
  temp.unique.nodes = unique(consensus.assignments)
  levels = sapply(1:length(temp.unique.nodes),function(t,i){length(strsplit(t[i],":")[[1]])},t=temp.unique.nodes)
  max.level = max(levels)
  unique.nodes=character(0)
  for(l in 1:max.level){
    nodes = temp.unique.nodes[levels==l]
    unique.nodes=c(unique.nodes,sort(nodes))
  }
  
  print("Getting consensus thetas")
  consensus.thetas=list()
  no.nodes = length(unique.nodes)
  is.single.point=vector(mode="logical",length=no.nodes)

  for(n in 1:no.nodes){
    print(paste(length(trees),n,sep=","))
    IDs = which(consensus.assignments==unique.nodes[n])
    if(length(IDs)==1){
      is.single.point[n]=T
    }
    consensus.theta = array(NA,c(no.subsamples,length(IDs),no.iters.post.burn.in))
    for(x in 1:no.subsamples){
      consensus.theta[x,,] = sapply((no.iters.burn.in+1):no.iters,function(t,a,p,i){t[[i]][match(a[p,i],t[[i]]$label),paste("theta.S",x,sep="")]},t=trees,a=node.assignments,p=IDs)				
    }
    consensus.thetas[[n]] = consensus.theta
  }

  ############################## PLOTTING ############################################################
  print("Starting plotting")
  #### HISTOGRAMS ####
  dev.set(which = plot.devs$hist.device)
  for(n in 1:no.nodes){
    for(p in 1:no.subsamples){
      #hist(consensus.thetas[[n]][p,,],breaks = seq(0,1,0.025),main=paste(samplename,subsamplenames[p]," ",unique.nodes[n],", (",dim(consensus.thetas[[n]])[2]," muts, ",no.nodes," nodes)",sep=""),xlab="subclonal fraction")
      if(is.na(tree.number)){
        hist(consensus.thetas[[n]][p,,],breaks = seq(0,1,0.025),main=paste(samplename,subsamplenames[p]," ",unique.nodes[n],", (",dim(consensus.thetas[[n]])[2]," muts)",sep=""),xlab="subclonal fraction")				
      }else{
        hist(consensus.thetas[[n]][p,,],breaks = seq(0,1,0.025),main=paste("tree #",tree.number,": ",samplename,subsamplenames[p]," ",unique.nodes[n],", (",dim(consensus.thetas[[n]])[2]," muts)",sep=""),xlab="subclonal fraction")				
      }
    }
  }
  ancs = paste(c("Root",sapply(2:no.nodes,function(n,i){spl = strsplit(n[i],":");paste(spl[[1]][1:(length(spl[[1]])-1)],collapse=":")},n=unique.nodes)),":",sep="")
  
  #### DENSITY ####
  dev.set(which = plot.devs$density.device)
  highest.density.thetas = array(NA,c(no.subsamples,no.nodes))
  print(start-Sys.time())
  
  # Estimate theta densities
  print("Estimating densities")
  for(n in 1:no.nodes){
    for(p in 1:no.subsamples){
      main = paste(samplename,subsamplenames[p]," ",unique.nodes[n],", (",dim(consensus.thetas[[n]])[2]," muts)",sep="")
      if(is.single.point[n]){
        temp.cons.thetas = array(consensus.thetas[[n]][p,,],c(1,no.iters.post.burn.in))
      } else {
        temp.cons.thetas = consensus.thetas[[n]][p,,]
      }
      
      if(!is.na(tree.number)){
        main = paste("tree #",tree.number,": ",main,sep="")
      }
      
      highest.density.thetas[p,n] = DensityEstimator(temp.cons.thetas,subclonal.fraction[,p][consensus.assignments==unique.nodes[n]],main=main,density.smooth=1,x.max=1.2)							
    }
  }
  
  #### TREE ####
  dev.set(which = plot.devs$tree.device)
  consensus.tree = as.data.frame(cbind(unique.nodes,t(highest.density.thetas),ancs,NA),stringsAsFactors=F)
  names(consensus.tree) = c("label",paste("theta.S",1:no.subsamples,sep=""),"ancestor","annotation")
  for(n in 1:no.subsamples){
    consensus.tree[,n+1] = as.numeric(consensus.tree[,n+1])
  }
  row.names(consensus.tree) = consensus.tree$label
  plotTree(consensus.tree,main=paste(samplename,subsamplenames,sep=""),font.size=1.25)
  
  #### OPTIMISED TREE ####
  #091013 - average theta, weighted by depth, which may be a better starting point for the whole.tree sampler
  #this tree should be recorded (and plotted)
  print("Building consensus trees")
  for (k in 1:no.subsamples) {
    for(n in 1:no.nodes){
      node.inds = which(consensus.assignments==unique.nodes[n])
      weights = mutCount[node.inds,k] + WTCount[node.inds,k]
      weights = weights/sum(weights)
      average.theta = mean(sapply((no.iters.burn.in+1):no.iters,function(i){sum(trees[[i]][node.assignments[node.inds,i],paste("theta.S",k,sep="")]*weights)}))			
      average.theta[is.na(average.theta)] = 0
      consensus.tree[unique.nodes[n],k+1] = average.theta
    }
  }
  print("average tree:")
  print(consensus.tree)
  
  #adjust thetas to make a valid tree, otherwise whole.tree.slice.sampler will be unable to find a solution
  consensus.tree = makeValidThetas(consensus.tree)
  print("adjusted average tree:")
  print(consensus.tree)
  
  # resample whole tree
  print("Resampling whole tree")
  no.slice.samples = 20
  capped.kappa = kappa
  capped.kappa[capped.kappa>0.999] = 0.999
  for (k in 1:no.subsamples) {
    curr.thetas = consensus.tree[,k+1]
    curr.thetas[curr.thetas<0.001] = 0.001
    curr.thetas[curr.thetas>0.999] = 0.999
    #consensus.tree[,k+1] <- whole.tree.slice.sampler(consensus.tree, curr.thetas, mutCount[,k], mutCount[,k] + WTCount[,k], capped.kappa[,k], consensus.assignments, shrinkage.threshold)
    #this may be slow, but it gives better estimates, plus, potentially, confidence intervals
    theta.samples = array(NA,c(no.slice.samples,no.nodes))
    
    call.whole.tree.slice.sampler = function(l,k1,consensus.tree1,curr.thetas1,mutCount1,WTCount1,capped.kappa1,consensus.assignments1,shrinkage.threshold1) { whole.tree.slice.sampler(consensus.tree1, curr.thetas1, mutCount1[,k1], mutCount1[,k1] + WTCount1[,k1], capped.kappa1[,k1], consensus.assignments1, shrinkage.threshold1) }
    
    theta.samples = t(sapply(1:no.slice.samples, FUN=call.whole.tree.slice.sampler, k=k, consensus.tree1=consensus.tree, curr.thetas1=curr.thetas, mutCount1=mutCount, WTCount1=WTCount, capped.kappa1=capped.kappa, consensus.assignments1=consensus.assignments, shrinkage.threshold1=shrinkage.threshold))
    
    consensus.tree[,k+1] = apply(theta.samples,2,median)
    print("median and stdev of consensus tree resamples:")
    print(cbind(consensus.tree[,k+1], apply(theta.samples,2,sd)))
  }
  dev.set(which=plot.devs$optimised.tree.device)
  plotTree(consensus.tree,main=paste(samplename,subsamplenames,sep=""),font.size=1.25)
  
  # check how many mutations were aggregated in the binning
  print("Checking aggregated mutations")
  if(!is.null(bin.indices)){
    print(paste("no.muts=",no.muts,sep=""))
    print(paste("length(bin.indices)=",length(bin.indices),sep=""))
    mut.table = vector(mode="numeric",length = no.nodes)
    mut.table2 = vector(mode="numeric",length = no.nodes)
    for(n in 1:no.nodes){
      node = unique.nodes[n]
      inds = which(consensus.assignments==node)
      mut.table[n] = sum(sapply(inds,function(b,i){length(b[[i]])},b=bin.indices))
    }
    mut.table = as.table(mut.table)
    names(mut.table) = unique.nodes
    
    disaggregated.consensus.assignments = rep(NA,sum(mut.table))
    for(m in 1:no.muts){
      disaggregated.consensus.assignments[bin.indices[[m]]] = consensus.assignments[m]
    }
  }else{
    mut.table = table(consensus.assignments[consensus.assignments!=""])
  }
  
  #### NO OF MUTS TREE ####
  dev.set(which=plot.devs$tree.population.device)
  population.tree=consensus.tree[,c("label","theta.S1","ancestor","annotation")]
  for(n in 1:no.nodes){
    #population.tree[unique.nodes[n],"theta.S1"] = mut.table[unique.nodes[n],1]
    population.tree[unique.nodes[n],"theta.S1"] = mut.table[unique.nodes[n]]
  }
  plotTree(population.tree,main=samplename,font.size=1.25,plotAsPercentage=F)
  
  if(no.subsamples>1){
    #its hard to distinguish more than 8 different colours
    max.cols = 8
    cols = rainbow(min(max.cols,no.nodes))
    dev.set(which = plot.devs$scatter.device)
    plot.data = subclonal.fraction
    plot.data[is.na(plot.data)]=0
    for(i in 1:(no.subsamples-1)){
      for(j in (i+1):no.subsamples){
        plot(plot.data[,i],plot.data[,j],type = "n",xlab = paste(samplename,subsamplenames[i]," subclonal fraction",sep=""), ylab = paste(samplename,subsamplenames[j]," subclonal fraction",sep=""),xlim = c(0,max(plot.data[,i])*1.2))
        for(n in 1:no.nodes){
          points(plot.data[,i][consensus.assignments==unique.nodes[n]],plot.data[,j][consensus.assignments==unique.nodes[n]],col=cols[(n-1) %% max.cols + 1],pch=20 + floor((n-1)/max.cols))
        }
        #legend(max(plot.data[,i]),max(plot.data[,j]),legend = unique.nodes[n],col=cols[(n-1) %% max.cols + 1],pch=20 + floor((n-1)/max.cols))
        legend(max(plot.data[,i]),max(plot.data[,j]),legend = unique.nodes,col=cols[(0:(no.nodes-1)) %% max.cols + 1],pch=20 + floor((0:(no.nodes-1))/max.cols),cex=0.5)
      }
    }
  }
  
  end_time = Sys.time()
  print(paste("Finished plotConsensusTrees in", as.numeric(end_time-start,units="secs"), "seconds"))
  #returned optimised tree and node assignments	
  if(!is.null(bin.indices)){
    return(list(consensus.tree,consensus.assignments,disaggregated.consensus.assignments))
  }else{
    return(list(consensus.tree,consensus.assignments))
  }
}