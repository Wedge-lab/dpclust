PhasingAssignment = function(mutCount=dataset$mutCount, 
                  WTCount=dataset$WTCount, 
                  no.iters, #input
                  num.timepoints = num.timepoints,
                  error.rate = error.rate,
                  killclu = kill.clu.frac,
                  cellularity=cellularity, 
                  conc_param = 1,
                  C = num_clu_samplers ,
                  all.phasing.index ,
                  all.phasing.reads ,
                  all.phasing.info ,
                  totalCopyNumber=dataset$totalCopyNumber, 
                  mutation.copy.number=dataset$mutation.copy.number,
                  copyNumberAdjustment=dataset$copyNumberAdjustment, 
                  mutationTypes=dataset$mutationType,
                  samplename=samplename, 
                  subsamplesrun=subsamples,
                  output_folder=outdir, 
                  cluster_conc=cluster_conc,
                  mut.assignment.type=mut.assignment.type,
                  most.similar.mut=most.similar.mut,
                  max.considered.clusters=max.considered.clusters,
                  normalCopyNumber=array(2,dim(mutCount))){

  parameter.setup = set.up.gibbs(B= conc_param ,
                                 C = C,
                                 num.timepoints = num.timepoints,
                                 no.iters = no.iters,
                                 cellularity = cellularity,
                                 mutCount = mutCount,
                                 WTCount = WTCount,
                                 totalCopyNumber = totalCopyNumber,
                                 normalCopyNumber = normalCopyNumber,
                                 copyNumberAdjustmen = copyNumberAdjustment)

  S.i = parameter.setup$S.i
  V.h = parameter.setup$V.h
  Pi.h = parameter.setup$Pi.h
  Pr.S = parameter.setup$Pr.S
  A = parameter.setup$A
  B = parameter.setup$B
  mutBurdens = parameter.setup$mutBurdens
  mutCopyNum = parameter.setup$mutCopyNum
  lower = parameter.setup$lower
  upper = parameter.setup$upper
  num.muts = parameter.setup$num.muts
  num.timepoints = parameter.setup$num.timepoints
  alpha = parameter.setup$alpha
  
  phased.pair.number = dim(all.phasing.index)[1]
  
  whole.tree = list()
  tree.win = rep(NA,no.iters)
  #initilize tree structure
  phy.clu.iter = array(NA,c(C,C,no.iters))
  phy.clu.iter[,,1] = initize.tree(C)
  
  clu.drop = 0 #cluster that contains less than specific number of mutations will be dropped
  #core sampler
  iter.count = 50
  for(m in 2:no.iters){
    if(m %% 50 == 0){
      print(paste(m ,"/",no.iters,sep=" "))
    }
    
    kill.index = 0 #indicate whether new cluster is killed during this iter
    
    phasing.pair.pro.temp = round(m/50)
    phasing.pair.used.length = min(round(phasing.pair.pro.temp*phased.pair.number*0.1),phased.pair.number) #how many pairs is used
    used.pairs.index = sample(1:phased.pair.number,phasing.pair.used.length )#which pairs we should use 
    #tree annealing to improve accuracy

    if(num.timepoints>1){
      tem = min(m,500)/500
    }else{
      tem = 1
    } #for Simulated Annealing in multiple sample analysis
    
    
    for(k in 1:num.muts){
      ass.prob = rep(0,C)
      Pr.S[k,1] <- log(V.h[m-1,1])
      Pr.S[k,2:C] = sapply(2:C, FUN=function(j, V) { log(V[j]) + sum(log(1-V[1:(j-1)])) }, V=V.h[m-1,])
      for(t in 1:num.timepoints){
        for(c in 1:C){
          ass.prob[c] <- ass.prob[c] + tem*mutCount[k,t]*log(mutBurdens[c,t,k]) + tem*WTCount[k,t]*log(1-mutBurdens[c,t,k])
        }
      }
      Pr.S[k,] = Pr.S[k,]+ass.prob
      Pr.S[k,is.na(Pr.S[k,])] = 0
      
      
      if(length(which(clu.drop!=0))>0){
        Pr.S[k,] = Pr.S[k,] - max(Pr.S[k,-clu.drop])
        Pr.S[k,] = exp(Pr.S[k,])
        Pr.S[k,clu.drop] = 0
      }else{
        Pr.S[k,] = Pr.S[k,] - max(Pr.S[k,])
        Pr.S[k,] = exp(Pr.S[k,]) 
      }
      Pr.S[k,] = Pr.S[k,] / sum(Pr.S[k,])
    }
    
    S.i[m,] = apply(Pr.S, 1, function(mut) { which(rmultinom(1,1,mut)==1) })
    
    
    if(iter.count > 49 & length(used.pairs.index)!=0){
      for(k in used.pairs.index){
        
        mut1.index = all.phasing.index[k,1]
        mut2.index = all.phasing.index[k,2]
        
        if(runif(1)>0.01){
          S.i[m,c(mut1.index,mut2.index)] = phasing.assign.onebyone(mut1.index,mut2.index,
                                                                    Pi.h, Pr.S, S.i, phy.clu.iter, all.phasing.reads,all.phasing.info, 
                                                                    C,k, m, clu.drop, error.rate, cellularity,num.timepoints)
        }else{
          S.i[m,c(mut1.index,mut2.index)] = phasing.assign.together(mut1.index,mut2.index,
                                                                    Pi.h, Pr.S, S.i, phy.clu.iter, all.phasing.reads,all.phasing.info, 
                                                                    C,k, m, clu.drop, error.rate, cellularity,num.timepoints)
        }
      }
    }
    
    # Update stick-breaking weights
    S.i[m,which(is.na(S.i[m,]))]=0
    V.h[m,1:(C-1)] = sapply(1:(C-1), function(S, curr.m, curr.alpha, h) {rbeta(1, 1+sum(S[curr.m,] == h), curr.alpha+sum(S[curr.m,] > h))}, S=S.i, curr.m=m, curr.alpha=alpha[m-1])
    V.h[m,c(V.h[m,1:(C-1)] == 1,FALSE)] = 0.999 # Need to prevent one stick from taking all the remaining weight
    
    # Get expected number of mutant reads per mutation copy number
    countsPerCopyNum = array(NA,c(num.timepoints,num.muts))
    for (t in 1:num.timepoints) {
      countsPerCopyNum[t,] = (mutCount[,t]+WTCount[,t])*mutationCopyNumberToMutationBurden(copyNumberAdjustment[,t],totalCopyNumber[,t],cellularity[t],normalCopyNumber[,t])
    }
    
    for(t in 1:num.timepoints){
      Pi.h[m,,t] = infer.pi.h(S.i,Pi.h, t, m,C, lower,upper, mutCount,countsPerCopyNum, tree.structure,num.timepoints)
    }
   
    
    for(t in 1:num.timepoints){
      for(c in 1:C){
        mutBurdens[c,t,]=mutationCopyNumberToMutationBurden(Pi.h[m,c,t] * copyNumberAdjustment[,t],totalCopyNumber[,t],cellularity[t],normalCopyNumber[,t])
      }
    }
    
    clu.max = max(unique(S.i[m,]))
    alpha[m] <- rgamma(1, shape=C+A-1, rate=B-sum(log(1-V.h[m,1:(clu.max-1)]))) 

    #drop cluster
    if(m>800){
      drop.result = drop.cluster(S.i,C, m, clu.drop, iter.count, killclu)
      clu.drop = drop.result[[1]]
      iter.count = drop.result[[2]]
    }else{
      clu.drop = 0
    }
    
    S.i[m,which(S.i[m,]%in%clu.drop)]=0 #mut move to unassigned statue
    
    #infer tree
    #step.1 

    if(iter.count <30){
      phy.clu.iter[,,m]=phy.clu.iter[,,m-1]
      tree.win[m] = tree.win[m-1]
      next
    }
    
    non.empty = c(1:C)
    non.empty = non.empty[which(!non.empty%in%clu.drop)]
    
    #new inferrence of tree

    if(length(non.empty)==2){
      clu.order = order(Pi.h[m,non.empty,1],decreasing=T)
      par = non.empty[clu.order[1]]
      dau = non.empty[clu.order[2]]
      tree.strcuture = c(par,dau)
      tree.win.temp=which(sapply(whole.tree,identical,tree.structure )==T)
      if(length(tree.win.temp)==0){
        whole.tree = c(whole.tree,list(tree.structure))
        tree.win.temp= length(whole.tree)
      }
      tree.win[m,1]=tree.win.temp
      
      whole.tree = c(whole.tree,list(c(par,dau)))
      phy.clu.iter[,,m] = 0
      phy.clu.iter[par,dau,m]=1
      phy.clu.iter[dau,par,m]=2
      phy.clu.iter[par,par,m] = phy.clu.iter[dau,dau,m] = 4
      next
    }

    if(num.timepoints>1){
      find.tree = find.tree.multiple(m=m,
                                     V.h,
                                     Pi.h,
                                     S.i,
                                     C,
                                     clu.drop = clu.drop,
                                     num.timepoints = num.timepoints,
                                     non.empty,
                                     all.phasing.index,
                                     all.phasing.reads,
                                     all.phasing.info,
                                     mutBurdens,
                                     mutCount,
                                     WTCount,
                                     error.rate,
                                     cellularity)
    }else{
      find.tree = find.tree.single(m = m,
                                   V.h = V.h,
                                   Pi.h = Pi.h,
                                   S.i = S.i,
                                   C = C,
                                   num.timepoints = num.timepoints,
                                   non.empty = non.empty,
                                   all.phasing.index = all.phasing.index,
                                   all.phasing.reads = all.phasing.reads,
                                   all.phasing.info = all.phasing.info,
                                   mutBurdens = mutBurdens,
                                   mutCount = mutCount,
                                   WTCount = WTCount,
                                   error.rate  = error.rate,
                                   cellularity = cellularity)
    }
    
    phy.clu.iter[,,m] = find.tree[[1]]
    tree.structure = find.tree[[2]]
    discard.clu = find.tree[[3]]
    tree.win.temp=which(sapply(whole.tree,identical,tree.structure )==T)
    if(length(tree.win.temp)==0){
      whole.tree = c(whole.tree,list(tree.structure))
      tree.win.temp= length(whole.tree)
    }
    tree.win[m]=tree.win.temp
  }
  return(list(whole.tree = whole.tree,
              tree.selected = tree.win,
              phylogeny.tree = phy.clu.iter,
              S.i = S.i,
              pi.h = Pi.h,
              V.h = V.h,
              mutBurdens = mutBurdens,
              alpha = alpha,
              mutCount = mutCount,
              WTCount = WTCount,
              N1=mutCount+WTCount))
}

