###############################
#set up gibbs sampler 

set.up.gibbs = function(A = 1,
                        B,
                        C ,
                        num.timepoints,
                        mutCount,
                        no.iters,
                        WTCount,
                        totalCopyNumber,
                        cellularity,
                        normalCopyNumber,
                        copyNumberAdjustment){
  A = 1
  B = conc_param 
  num.muts = NROW(mutCount)
  num.timepoints = t = NCOL(mutCount)
  #setup 
  
  Pi.h = array(NA, c(no.iters,C,num.timepoints))
  mutBurdens = array(NA, c(C,num.timepoints,num.muts))
  Pr.S = array(NA,c(num.muts,C))
  V.h = matrix(1, nrow=no.iters, ncol=C)
  S.i = matrix(NA, nrow= no.iters, ncol=num.muts)
  alpha = rep(NA,  no.iters)
  
  lower = array(NA,num.timepoints)
  upper = array(NA,num.timepoints)
  
  mutCopyNum = array(NA,c(num.muts,num.timepoints))
  
  for (t in 1:num.timepoints) {  
    mutCopyNum[,t] = mutationBurdenToMutationCopyNumber(mutCount[,t]/(mutCount[,t]+WTCount[,t]),totalCopyNumber[,t],cellularity[t],normalCopyNumber[,t]) / copyNumberAdjustment[,t]
    lower[t] = min(mutCopyNum[,t])
    upper[t] = max(mutCopyNum[,t])
    difference = upper[t]-lower[t]
    lower[t] = lower[t]-difference/10
    upper[t] = upper[t]+difference/10
    # Randomise starting positions of clusters
    set.seed(seed)
    Pi.h[1,,t]=runif(C,lower[t],upper[t])
    for(c in 1:C){
      mutBurdens[c,t,]=mutationCopyNumberToMutationBurden(Pi.h[1,c,t] * copyNumberAdjustment[,t], totalCopyNumber[,t], cellularity[t],normalCopyNumber[,t]) 
    }
  }	
  
  
  V.h[1,] = c(rep(0.5,C-1), 1)
  V.h[1: no.iters, C] = rep(1,  no.iters)
  alpha[1] = 1
  a = 1
  S.i[1, ] = sample(c(1:C), num.muts, replace = T)
  
  return(list(S.i=S.i, V.h=V.h, Pi.h=Pi.h, Pr.S = Pr.S, A = A, B = B, alpha = alpha, 
              mutBurdens = mutBurdens, mutCopyNum = mutCopyNum,lower = lower,upper = upper,
              num.muts = num.muts, num.timepoints = num.timepoints))
  
}

###############################
#randomize phylogeny

initize.tree = function(C){
  
  phy.1 = matrix(NA,C,C)
  for(clu.1 in 1:(C-1)){
    phy.1[clu.1,clu.1] = 4
    phy.1[clu.1,(clu.1+1):C] = sample(c(1,2,3),length((clu.1+1):C),replace =T)
  }
  phy.1[C,C] = 4
  for(clu.1 in 2:C){
    for(clu.2 in 1:(clu.1-1)){
      phy.1[clu.1,clu.2] = 3-phy.1[clu.2,clu.1]
      if(phy.1[clu.1,clu.2]==0){
        phy.1[clu.1,clu.2] = 3
      }
    } 
  }
  return(phy.1)
}

########################################
#phasing assignment likelihood calculation
#@ K: index for phasing pairs
phasing.assign.onebyone = function(mut1.index, mut2.index,
                                   Pi.h, Pr.S, S.i, phy.clu.iter, all.phasing.reads, all.phasing.info, 
                                   C, k, m, clu.drop,
                                   error.rate, cellularity,num.timepoints){
  
  mut1.cluster = S.i[m-1,mut1.index]
  mut2.cluster = S.i[m-1,mut2.index]
  
  if(mut1.cluster == 0 | mut2.cluster == 0){
    pha.1.ass = S.i[m,mut1.index]
    pha.2.ass = S.i[m,mut2.index]
    return(c(pha.1.ass,pha.2.ass))
  }#in case we don't use them in burn in period 

  read.mut.mut.temp = all.phasing.reads[k,2,]
  cp.1 = cp.2 = all.phasing.info$Copy_Major[k]
  cp.total = all.phasing.info$Copy_Minor[k] + all.phasing.info$Copy_Major[k]
  if(length(which(read.mut.mut.temp>0))>0 | cp.total==1 ){
    same.chr = TRUE
  }else{
    same.chr = FALSE
  }
  
  for(t in 1:num.timepoints){
      phasing.read = as.numeric(all.phasing.reads[k,,t])
      mut.index.input= mut1.index
      ubb = Pr.S[mut1.index,]
      pha.likeli = rep(NA,C)
      for (c in 1:C){
        j = phy.clu.iter[c,mut2.cluster,m-1]
        if(j==0){
          Pr.S[mut1.index,c] = 0 
          next
        }

        pha.likeli[c] = infer.likelihood.same(Pi.h = Pi.h,m = m-1, t=t,
                                                     error.rate = error.rate,
                                                     phy = j,
                                                     cellularity=cellularity,
                                                     mut1.cluster = c,
                                                     mut2.cluster = mut2.cluster,
                                                     mut.index.input = mut.index.input,
                                                     phasing.read = phasing.read, 
                                                     copy.1 = cp.1, copy.2 = cp.2, copy.total = cp.total) 
      
        if(!same.chr){
          pha.likeli[c] = 0.5*pha.likeli[c]+0.5*infer.likelihood.same(Pi.h = Pi.h,m = m-1,t=t,
                                                                         error.rate = error.rate,
                                                                         phy = 3,
                                                                         cellularity=cellularity,
                                                                         mut1.cluster = c,
                                                                         mut2.cluster = mut2.cluster,
                                                                         mut.index.input = mut.index.input,
                                                                         phasing.read = phasing.read, 
                                                                         copy.1 = cp.1, copy.2 = cp.2, copy.total = cp.total) 
        }
        
        #print(pha.likeli[mut1.index,c,m])
        Pr.S[mut1.index,c] = Pr.S[mut1.index,c]*pha.likeli[c]
      }
      
  }
  
  if(length(which((Pr.S[mut1.index,]==0)))==C){
    Pr.S[mut1.index,]=rep(1,C)
  }
  Pr.S[mut1.index,] = Pr.S[mut1.index,] / sum(Pr.S[mut1.index,])
  Pr.S[mut1.index,clu.drop] = 0
  
  pha.1.ass = which(rmultinom(1,1,Pr.S[mut1.index,])==1)  
  
  #assign 2nd
  
  mut1.cluster = S.i[m-1,mut2.index]
  mut2.cluster = pha.1.ass

  ubb = Pr.S[mut2.index,]
  pha.likeli = rep(NA,C)
  
  for(t in 1:num.timepoints){
    phasing.read = as.numeric(all.phasing.reads[k,c(1,2,4,3),t])
    mut.index.input= mut2.index
    ubb = Pr.S[mut2.index,]
    pha.likeli = rep(NA,C)
    cp.1 = cp.2 = all.phasing.info$Copy_Major[k]
    cp.total = all.phasing.info$Copy_Minor[k] + all.phasing.info$Copy_Major[k]
    for (c in 1:C){
      j = phy.clu.iter[c,mut2.cluster,m-1]
      if(j==0){
        Pr.S[mut2.index,c] = 0 
        next
      }
      pha.likeli[c] = infer.likelihood.same(Pi.h = Pi.h,m = m-1,t=t,
                                                   error.rate = error.rate,
                                                   phy = j,
                                                   cellularity=cellularity,
                                                   mut1.cluster = c,
                                                   mut2.cluster = mut2.cluster,
                                                   mut.index.input = mut.index.input,
                                                   phasing.read = phasing.read, 
                                                   copy.1 = cp.1, copy.2 = cp.2, copy.total = cp.total)   
      if(!same.chr){
        pha.likeli[c] = 0.5*pha.likeli[c]+0.5*infer.likelihood.same(Pi.h = Pi.h,m = m-1,t=t,
                                                                       error.rate = error.rate,
                                                                       phy = 3,
                                                                       cellularity=cellularity,
                                                                       mut1.cluster = c,
                                                                       mut2.cluster = mut2.cluster,
                                                                       mut.index.input = mut.index.input,
                                                                       phasing.read = phasing.read, 
                                                                       copy.1 = cp.1, copy.2 = cp.2, copy.total = cp.total) 
      }
      
      #print(pha.likeli[mut1.index,c,m])
      Pr.S[mut2.index,c] = Pr.S[mut2.index,c]*pha.likeli[c]
    }
    
  }
  
  if(length(which((Pr.S[mut2.index,]==0)))==C){
    Pr.S[mut2.index,]=rep(1,C)
  }
  Pr.S[mut2.index,] = Pr.S[mut2.index,] / sum(Pr.S[mut2.index,])
  Pr.S[mut2.index,clu.drop] = 0
  pha.2.ass = which(rmultinom(1,1,Pr.S[mut2.index,])==1)  
  
  return(c(pha.1.ass, pha.2.ass))
}


########################################
#phasing assignment likelihood calculation
phasing.assign.together = function(mut1.index,mut2.index,
                                   Pi.h, Pr.S, S.i, phy.clu.iter, all.phasing.reads, all.phasing.info, 
                                   C,k, m, clu.drop, error.rate, cellularity,num.timepoints){

  Pr.S[mut1.index,clu.drop] = 0
  Pr.S[mut2.index,clu.drop] = 0
  Prs.together = matrix(NA,C,C)
  
  cp.1 = cp.2 = all.phasing.info$Copy_Major[k]
  cp.total = all.phasing.info$Copy_Minor[k] + all.phasing.info$Copy_Major[k]
  
  read.mut.mut.temp = all.phasing.reads[k,2,]
  
  if(length(which(read.mut.mut.temp>0))>0 | cp.total==1){
    same.chr = TRUE
  }else{
    same.chr = FALSE
  }
  
  for(t in 1:num.timepoints){
    phasing.read = as.numeric(all.phasing.reads[k,,t])
    for(c.1 in 1:C){
      for(c.2 in 1:C){
        j = phy.clu.iter[c.1,c.2,m-1]
        if(j==0){
          Prs.together[c.1,c.2] = 0 
          next
        }
        phalike.together = infer.likelihood.same(Pi.h = Pi.h,m = m-1,t=t,
                                                        cellularity=cellularity,
                                                        error.rate =error.rate,
                                                        phy = j,
                                                        mut1.cluster = c.1,
                                                        mut2.cluster = c.2,
                                                        mut.index.input = mut.index.input,
                                                        phasing.read = phasing.read,
                                                        copy.1 = cp.1,copy.2 = cp.2,copy.total = cp.total)  
        
        if(!same.chr){
          phalike.together = 0.5*phalike.together+0.5*infer.likelihood.same(Pi.h = Pi.h,m = m-1,t=t,
                                                                               cellularity=cellularity,
                                                                               error.rate =error.rate,
                                                                               phy = 3,
                                                                               mut1.cluster = c.1,
                                                                               mut2.cluster = c.2,
                                                                               mut.index.input = mut.index.input,
                                                                               phasing.read = phasing.read,
                                                                               copy.1 = cp.1,copy.2 = cp.2,copy.total = cp.total)
        }
        
        Prs.together[c.1,c.2] = Pr.S[mut1.index,c.1]*Pr.S[mut2.index,c.2]*phalike.together 
      }
    }
  }
  
  prob.tem = as.vector(Prs.together)
  if(length(which((prob.tem==0)))==C*C){
    prob.tem=rep(1,C*C)
  }
  ass.tem =  which(rmultinom(1,1,prob.tem)==1) 
  pha.2.ass = ceiling(ass.tem/C)
  s1.tem = ass.tem%%C
  if(s1.tem==0){
    s1.tem=C
  }
  pha.1.ass = s1.tem
  
  return(c(pha.1.ass, pha.2.ass))
}


########################################
#drop cluster
drop.cluster = function(S.i,C, m, clu.drop, iter.count, killclu ){
  clu.size = sapply(1:C,function(t){length(which(S.i[m,]==t))})
  clu.size.pro = clu.size/sum(clu.size)
  if(0%in%clu.drop){
    clu.drop = clu.drop[-which(clu.drop==0)]
  }
  
  iter.count = iter.count+1
  clu.drop.temp = which(clu.size.pro<killclu) #cluster contains less than 1% mutations
  new.drop = clu.drop.temp[which(!clu.drop.temp%in%clu.drop)]
  
  if(iter.count<50){
    clu.drop = clu.drop
  }else{
    if(length(new.drop)==1){
      clu.drop = clu.drop.temp
      iter.count=0
      kill.index = 1
    }
    if(length(new.drop)>1){
      new.drop = sample(new.drop,1)
      clu.drop = c(clu.drop,new.drop)
      iter.count=0
      kill.index = 1
    }
    if(length(clu.drop)==0){
      clu.drop = 0
    }
  }
  return(c(list(clu.drop), list(iter.count)))
}


########################################
#infer Pi.h
infer.pi.h = function(S.i,Pi.h,t, m,C, lower,upper,mutCount,countsPerCopyNum, tree.structure,num.timepoints){
  if(m>500 & num.timepoints == 1){
    #sampler.count = 0
    #for(sampler.iter in 1:2000){
    stru.constrain=1
    pi.count=0
    while(stru.constrain){
      stru.constrain=0
      pi.count=pi.count+1
      if(pi.count>20000){
        break
      }
      
      ccf.infer = runif(C,lower[t],upper[t])

      for(c in unique(S.i[m,])){
        #040213 - problem if sum(countsPerCopyNum[t,S.i[m,]==c])==0 fixed
        if(sum(countsPerCopyNum[t,S.i[m,]==c])==0){
          ccf.infer[c] = 0
        }else{
          ccf.infer[c] = rgamma(1,shape=sum(mutCount[S.i[m,]==c,t]),rate=sum(countsPerCopyNum[t,S.i[m,]==c]))
        }
      }
      
      if(length(which(is.na(ccf.infer))>0)){
        stru.constrain=1
      }
      for(c.1 in unique(S.i[m,])){
        c.2 = tree.structure[which(tree.structure[,1]==c.1),2]
        if(ccf.infer[c.1]<sum(ccf.infer[c.2])){
          stru.constrain=1
        }        
      }
    }
    
  }else{
    ccf.infer = runif(C,lower[t],upper[t])
    
    for(c in unique(S.i[m,])){
      if(sum(countsPerCopyNum[t,S.i[m,]==c])==0){
        ccf.infer[c] = 0
      }else{
        ccf.infer[c] = rgamma(1,shape=sum(mutCount[S.i[m,]==c,t]),rate=sum(countsPerCopyNum[t,S.i[m,]==c]))
      }
      
    }	
  }
  
  return(ccf.infer)
}

#########################################
#function for calculation likelihood
infer.likelihood.same = function(Pi.h, m,t,phy,cellularity,
                                        mut1.cluster,mut2.cluster,mut.index.input,phasing.read,
                                        copy.1,copy.2,copy.total,error.rate){
  ccf.mut1 = Pi.h[m,mut1.cluster,t]*cellularity[t]
  ccf.mut2 = Pi.h[m,mut2.cluster,t]*cellularity[t]
  
  purity = cellularity[t]
  copy.bear.1 = copy.1
  copy.bear.2 = copy.2
  copy.total = copy.total
  p.mut.mut = p.wt.mut = p.mut.wt = p.wt.wt = 0
  
  if(phy == 4){
    ccf.mut = (ccf.mut1+ccf.mut2)/2
    
    if(copy.bear.1 == copy.bear.2){
      p.mut.wt = p.wt.mut = 0
    }else{
      p.wt.mut = ccf.mut*max(0,copy.bear.2-copy.bear.1)
      p.mut.wt = ccf.mut*max(0,copy.bear.1-copy.bear.2)
    }
    
    p.wt.wt = 1- ccf.mut + ccf.mut*(copy.total-max(copy.bear.1,copy.bear.2))/copy.total
    p.mut.mut = ccf.mut*min(copy.bear.1,copy.bear.2)/copy.total
  }
  if(phy == 3){
    
    p.mut.mut = 0
    p.wt.mut = ccf.mut2*copy.bear.2/copy.total
    p.mut.wt = ccf.mut1*copy.bear.1/copy.total
    p.wt.wt = 1 - p.mut.wt - p.wt.mut
  }
  if(phy == 2){
    p.mut.wt = 0
    p.wt.mut = (ccf.mut2-ccf.mut1)*copy.bear.2/copy.total
    p.mut.mut = ccf.mut1*min(copy.bear.1,copy.bear.2)/copy.total
    p.wt.wt = 1 - p.mut.wt - p.wt.mut-p.mut.mut
  }
  if(phy == 1){
    p.wt.mut = 0
    p.mut.wt = (ccf.mut1-ccf.mut2)*copy.bear.1/copy.total
    p.mut.mut = ccf.mut2*min(copy.bear.1,copy.bear.2)/copy.total
    p.wt.wt = 1 - p.mut.wt - p.wt.mut-p.mut.mut
  }
  
  pro = c(p.wt.wt,p.mut.mut,p.wt.mut,p.mut.wt)
  pro[which(pro<0)] = 0
  pro[which(pro>1)] = 1
  p.wt.wt = pro[1]
  p.mut.mut = pro[2]
  p.wt.mut = pro[3]
  p.mut.wt = pro[4]
  
  p.wt.wt.error = p.wt.wt*(1-error.rate)*(1-error.rate) + p.wt.mut*error.rate + p.mut.wt*error.rate*(1-error.rate) + p.mut.mut*error.rate*error.rate
  p.mut.mut.error = p.wt.wt*(error.rate)*(error.rate) + (p.wt.mut+p.mut.wt)*error.rate*(1-error.rate) + p.mut.mut*(1-error.rate)*(1-error.rate)
  p.wt.mut.error = p.wt.wt*(1-error.rate)*error.rate + p.wt.mut*(1-error.rate)*(1-error.rate) + p.mut.wt*error.rate*error.rate + p.mut.mut*error.rate*(1-error.rate)
  p.mut.wt.error = p.wt.wt*error.rate*(1-error.rate) + p.wt.mut*error.rate*error.rate + p.mut.wt*(1-error.rate)*(1-error.rate) + p.mut.mut*(1-error.rate)*error.rate
  
  pro.error = c(p.wt.wt.error,p.mut.mut.error,p.wt.mut.error,p.mut.wt.error)
  pha.pro = dmultinom(x=phasing.read,prob = pro.error)
  
  return(pha.pro)
  
}


########################################
#function to grow tree

find.tree.single= function(m,
                          V.h,
                          Pi.h,
                          S.i,
                          C,
                          num.timepoints,
                          non.empty,
                          all.phasing.index,
                          all.phasing.reads,
                          all.phasing.info,
                          mutBurdens,
                          mutCount,
                          WTCount,
                          error.rate,
                          cellularity){
  
  CCFs = Pi.h[m,,]
  CCF.order = order(CCFs,decreasing=T)
  CCF.order = CCF.order[which(CCF.order%in%non.empty)] 
  
  tree.temp = grow.tree(CCF.order,Pi.h,m,C,num.timepoints,soft.threshold=0)  
  tree.str.temp = tree.temp$Tree
  tree.phy.temp = tree.temp$Phy 
  
  
  tree.win.index = select.tree.likelihood(m,V.h, Pi.h, S.i, num.timepoints,
                                              used.clu = non.empty,C,mutBurdens,
                                              all.phasing.index,all.phasing.info,all.phasing.reads,
                                              tree.str.temp, tree.phy.temp,sample.tree=F,error.rate=error.rate)
  
 
  
  tree.phy = tree.phy.temp[,,tree.win.index]
  tree.structure =  tree.str.temp[,,tree.win.index]
  discard.clu = c()
  
  return(list(Phy = tree.phy,
              Str = tree.structure,
              discard.clu))
}


find.tree.multiple = function(m,
                              V.h,
                              Pi.h,
                              S.i,
                              C,
                              num.timepoints,
                              non.empty,
                              clu.drop,
                              all.phasing.index,
                              all.phasing.reads,
                              all.phasing.info,
                              mutBurdens,
                              mutCount,
                              WTCount,
                              error.rate,
                              cellularity)
{ 
  clu.num = length(non.empty)
  clu.size = sapply(non.empty,function(t){length(which(S.i[m,]==t))})
  index.matrix = matrix(NA,nrow=2,ncol=clu.num)
  index.matrix[1,] = 1
  index.matrix[2,] = 0
  indi.matrix = expand.grid(as.data.frame(index.matrix))
  indi.clu.num = apply(indi.matrix,1,sum) 
  indi.clu.size = apply(indi.matrix,1,function(t){sum(clu.size*t)}) 
  clu.comb.num = order(indi.clu.num,indi.clu.size,decreasing=F) 
    # if cannot find a tree including all clusters, we drop clusters
  indi.matrix = indi.matrix[clu.comb.num,] # rank the num of combination of all clusters 

  whole.sample.tree = list()
  #finding all tree
  
  for(i.row in 1:dim(indi.matrix)[1]){
    used.clu = non.empty[which(indi.matrix[i.row,]==0)]
    CCFs.sort = get.order(used.clu,m,Pi.h)
    tree.temp = grow.tree(CCFs.sort,Pi.h,m,C,num.timepoints,soft.threshold = -0.05)
    tree.str.temp = tree.temp$Tree
    tree.phy.temp = tree.temp$Phy
    if(length(tree.str.temp)>1){
      break
    }
  }
  #likeli.tree = array(NA,c(length(wholesample.tree.highlikeli.index),num.timepoints))
  
  tree.win.index = select.tree.likelihood(m,V.h, Pi.h, S.i, num.timepoints,
                                          used.clu,C,mutBurdens,
                                          all.phasing.index,all.phasing.info,all.phasing.reads,
                                          tree.str.temp, tree.phy.temp,sample.tree=F,error.rate=error.rate)
  
  tree.phy = tree.phy.temp[,,tree.win.index]
  tree.structure =  tree.str.temp[,,tree.win.index]
  discard.clu = c(non.empty[which(indi.matrix[i.row,]==1)],clu.drop) #not included in trees
  
  return(list(Phy = tree.phy,
              Str = tree.structure,
              discard.clu))
  
}

select.tree.likelihood = function(m, V.h, Pi.h, S.i, num.timepoints,
                                     used.clu,C,mutBurdens,
                                     all.phasing.index,all.phasing.info,all.phasing.reads,
                                     tree.str.temp, tree.phy.temp,sample.tree=T,error.rate)
{
  vh.temp = rep(NA,C)
  vh.temp[1] = log(V.h[m,1])
  vh.temp[2:C] = sapply(2:C, FUN=function(j, V) { log(V[j]) + sum(log(1-V[1:(j-1)])) }, V=V.h[m,])
  
  
  phasing.pair.num = dim(all.phasing.index)[1]
  fre.pair.temp = array(NA,c(phasing.pair.num,C,2))
  phasing.pair.temp = array(NA,c(phasing.pair.num,C,C,4,num.timepoints))
  
  for(i.pair in 1:phasing.pair.num){
    same.chr = FALSE
    mut1.index = all.phasing.index[i.pair,1]
    mut2.index = all.phasing.index[i.pair,2] # the index stored in the last two columns, and is identical across all samples
    read.mut.mut.temp = all.phasing.reads[i.pair,2,]
    
    cp.1 = cp.2 = all.phasing.info$Copy_Major[i.pair]
    cp.total = all.phasing.info$Copy_Minor[i.pair]+all.phasing.info$Copy_Major[i.pair]
    
    if(length(which(read.mut.mut.temp>0))>0 | cp.total == 1){
      same.chr = TRUE
    }
    
    fre.pair.temp[i.pair,used.clu,1] = sapply(used.clu,function(c.1){exp(vh.temp[c.1])})
    fre.pair.temp[i.pair,used.clu,2] = sapply(used.clu,function(c.2){exp(vh.temp[c.2])})
    
    for(t in 1:num.timepoints){
      
      phasing.read = all.phasing.reads[i.pair,1:4,t]
      fre.pair.temp[i.pair,used.clu,1] = fre.pair.temp[i.pair,used.clu,1]*sapply(used.clu,function(c.1){exp(mutCount[mut1.index,t]*log(mutBurdens[c.1,t,mut1.index])+
                                                                                  WTCount[mut1.index,t]*log(1-mutBurdens[c.1,t,mut1.index]))})
      fre.pair.temp[i.pair,used.clu,2] = fre.pair.temp[i.pair,used.clu,2]*sapply(used.clu,function(c.2){exp(mutCount[mut2.index,t]*log(mutBurdens[c.2,t,mut1.index])+
                                                                                WTCount[mut2.index,t]*log(1-mutBurdens[c.2,t,mut1.index]))})    
      
      for(c.1 in used.clu){
        for(phy.iter in 1:4){
            phasing.pair.temp[i.pair,c.1,used.clu,phy.iter,t] = sapply(used.clu,function(c.2){infer.likelihood.same(Pi.h, m,t,phy = phy.iter,cellularity,
                                                                                                                     mut1.cluster=c.1,mut2.cluster=c.2,mut.index.input,phasing.read,
                                                                                                                     copy.1=cp.1,copy.2=cp.2,copy.total=cp.total,error.rate)})
        }
        if(!same.chr){
          phasing.pair.temp[i.pair,c.1,used.clu,,t]=0.5*phasing.pair.temp[i.pair,c.1,used.clu,,t]+0.5*phasing.pair.temp[i.pair,c.1,used.clu,3,t]        
        }
      }
    }    
  }
  
  num.tree = dim(tree.str.temp)[3]
  likeli.tree = rep(NA,num.tree)
  
  for(i.tree in 1:num.tree){
    phy.temp = tree.phy.temp[,,i.tree]
    phasing.cal.temp = array(NA,c(phasing.pair.num,C,C))
    for(c.1 in used.clu){
      for(c.2 in used.clu){
        if(num.timepoints>1){
          phasing.cal.temp[,c.1,c.2] = apply(phasing.pair.temp[,c.1,c.2,phy.temp[c.1,c.2],1:num.timepoints],1,prod)
        }else{
          phasing.cal.temp[,c.1,c.2] = phasing.pair.temp[,c.1,c.2,phy.temp[c.1,c.2],num.timepoints]
          
        }
      }
    }
    tree.pair.temp = sapply(1:phasing.pair.num,function(i){ log(sum(fre.pair.temp[i,used.clu,1]%*%t(fre.pair.temp[i,used.clu,2])* phasing.cal.temp[i,used.clu,used.clu]))})    
    tree.pair.temp = tree.pair.temp[which(tree.pair.temp!=-Inf)]
    likeli.tree[i.tree] = sum(tree.pair.temp)
  }
  
  likeli.tree.multiplesample = likeli.tree - max(likeli.tree)
  likeli.tree.multiplesample= exp(likeli.tree.multiplesample)
  
  if(sample.tree ){
    tree.win.index = sample(c(1:length( likeli.tree.multiplesample)),1,prob= likeli.tree.multiplesample)
  }else{
    tree.win.index = which(likeli.tree.multiplesample == max(likeli.tree.multiplesample))
  }
  
  if(length(tree.win.index)>1){
    tree.win.index = sample(tree.win.index,1)
  }
  
  return(tree.win.index)
}


grow.tree = function(CCF.order,Pi.h,m,C,num.timepoints,soft.threshold){
  tree.start = matrix(NA,nrow = length(CCF.order)-1,ncol = 2)
  tree.start[1,] = c(CCF.order[1],CCF.order[2])
  colnames(tree.start) = c("Ancestors","Descendants")
  
  tree.grow = array(NA,c(nrow = length(CCF.order)-1,ncol = 2,1))
  tree.grow[,,1] = tree.start
  
  
  for(clu.index in 3:length(CCF.order)){
    clu.used = CCF.order[clu.index]#add this new cluster to tree
    j=0
    if(length(dim(tree.grow))==2){
      tree.grow.temp = array(NA,c(length(CCF.order)-1,2,(clu.index-1)))
      tree.grow = array(tree.grow,c(dim(tree.grow),1))
    }else{
      tree.grow.temp = array(NA,c(length(CCF.order)-1,2,dim(tree.grow)[3]*(clu.index-1)))
    }

    for(tree.used.index in 1:dim(tree.grow)[3]){
      tree.used = tree.grow[,,tree.used.index]
      tree.new = tree.used
      clu.tree = unique(c(tree.used))
      clu.tree = clu.tree[-which(is.na(clu.tree))]
      
      for (ance in clu.tree){
        clu.sib = tree.used[which(tree.used[,1]==ance),2]
        if(length(clu.sib) & num.timepoints>1){
          ccf.sib.sum = apply(Pi.h[m,c(clu.sib,clu.used),],2,sum)
        }else{
          ccf.sib.sum = sum(Pi.h[m,c(clu.sib,clu.used),])
        }
        
        indi.diff = Pi.h[m,ance,] - ccf.sib.sum
        
        if(length(which(indi.diff< soft.threshold))){
          next
        }#piegon hole principle, but use soft rules in multiple samples
        
        j=j+1
        tree.new[which(is.na(tree.used[,1]))[1],] = c(ance,clu.used) 
        tree.grow.temp[,,j] = tree.new
      }
    }
    tree.grow = tree.grow.temp[,,1:j]
  }
  if(j==0){
    return(list())
  }
  
  if(length(dim(tree.grow))==2){
    tree.grow = array(tree.grow,c(dim(tree.grow),1))
  }

  phy.tree.infer = array(NA,dim=c(C,C,dim(tree.grow)[3]))
  modify.tree.iter = array(NA,dim = c(dim(tree.start),dim(tree.grow)[3]))
  
  all.clu = c(1:C)
  clu.drop = all.clu[which(!all.clu%in%CCF.order)]
  
  for(i in 1:dim(tree.grow)[3]){
    modify.tree =tree.grow[,,i]
    modify.tree=modify.tree[order(modify.tree[,1],modify.tree[,2]),]
    modify.tree.iter[,,i]=modify.tree 
    
    phy.temp = matrix(NA,C,C)
    tree.str = list()
    for(clu in 1:C){ #find parents and grandparents
      clu.temp = clu
      tree.str.temp = c()
      while(length(which(modify.tree[,2]==clu.temp))){
        clu.temp = modify.tree[which(modify.tree[,2]==clu.temp),1]
        tree.str.temp = c(tree.str.temp,clu.temp)
      }
      tree.str = c(tree.str,list(tree.str.temp)) #tree.structure 
    }
    for(j in 1:C){
      phy.temp[tree.str[[j]],j]=1
      phy.temp[j,tree.str[[j]]]=2
      phy.temp[j,j] = 4
    }
    phy.temp[which(is.na(phy.temp))] = 3
    phy.temp[clu.drop,]=0
    phy.temp[,clu.drop]=0 #tree structure in matrix
    phy.tree.infer[,,i] = phy.temp
  } #find phylogeny matrix
  list(Tree=modify.tree.iter,Phy=phy.tree.infer)
}
get.order = function(used.clu, m, Pi.h){
  used.clu.num = length(used.clu)
  
  ccf.order = c(used.clu[1])
  for(c.index in 2:used.clu.num){
    c = used.clu[c.index]
    indi = rep(F,length(ccf.order))
    for(c.temp.index in 1:length(ccf.order)){
      c.temp = ccf.order[c.temp.index]
      #diff.clu = Pi.h[(m-100):(m),c,] - Pi.h[(m-100):(m),c.temp,]
      #sd.temp = apply(diff.clu,2, sd)
      indi.diff = Pi.h[m,c.temp,] - Pi.h[m,c,] 
      if(length(which(indi.diff< -0.05))){
        indi[c.temp.index] = T
      }
    }
    indi.f = which(!indi)
    if(length(indi.f)==0){
      ccf.order = c(c,ccf.order)
    }else{
      indi.max.f = max(indi.f)
      if(indi.max.f==length(ccf.order)){
        ccf.order = c(ccf.order,c)
      }else{
        ccf.order = c(ccf.order[1:indi.max.f],c,ccf.order[(indi.max.f+1):length(ccf.order)])
      }
      
    }
  }
  
  return(ccf.order)
}


  