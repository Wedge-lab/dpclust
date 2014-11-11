source("Tree_based_DP_Gibbs_sampler.R")
# source("interconvertMutationBurdens.R")
source("CullTree.R")
source("PlotConsensusTree.R")
source("GetConsensusStartPosition.R")

GetConsensusTrees<-function(trees, node.assignments, mutCount, WTCount, kappa=array(0.5,dim(mutCount)), samplename="sample", subsamplenames=1:ncol(mutCount), no.iters=1250, no.iters.burn.in=250, resort.mutations=T, shrinkage.threshold=0.1, init.alpha=0.01, outdir=getwd(), bin.indices=NULL){
  start = Sys.time()
  
  setwd(outdir)
	
	no.subsamples = ncol(mutCount)
	no.muts = nrow(mutCount)
  no.iters.post.burn.in = no.iters-no.iters.burn.in

	# Assemble the current strengths for each pair of mutations
  print("Calculating mutation strengths")
  strengths = get.mut.ass.strengths(no.muts, no.iters, no.iters.post.burn.in, node.assignments)
  ancestor.strengths = strengths$ancestor.strengths
  sibling.strengths = strengths$sibling.strengths
  identity.strengths = strengths$identity.strengths

  # Save the strengths as output
  write.table(ancestor.strengths, file='ancestor.strengths.csv', row.names=F, sep=",")
  write.table(sibling.strengths, file='sibling.strengths.csv', row.names=F, sep=",")
  write.table(identity.strengths, file='identity.strengths.csv', row.names=F, sep=",")
  
  print(Sys.time()-start)

  # Create devices to push plots towards lateron
  plot.devs = create.plotting.devices(samplename, no.iters, no.iters.burn.in, no.subsamples)

  print("Obtaining current agreement")
  consensus.assignments = rep("M:",no.muts)
#   subclonal.fraction = mutCount / ((mutCount+WTCount)*kappa)
#   subclonal.fraction[kappa == 0] = NA
#   start.pos = find.cons.start.pos(subclonal.fraction, identity.strengths, sibling.strengths, ancestor.strengths, 0.1)
#   consensus.assignments = start.pos$assignments
  
  current.agreement = calc.curr.agreement(consensus.assignments, identity.strengths, ancestor.strengths, sibling.strengths)
	fractional.current.agreement = current.agreement/(no.muts*no.muts*no.iters.post.burn.in)

  print(Sys.time()-start)
	print(paste("M:",current.agreement,fractional.current.agreement))
	
	res = do_em(trees=trees,
	            node.assignments=node.assignments, 
	            ancestor.strengths=ancestor.strengths, 
	            sibling.strengths=sibling.strengths, 
              identity.strengths=identity.strengths, 
	            bin.indices=bin.indices, 
	            consensus.assignments=consensus.assignments, 
	            current.agreement=current.agreement, 
	            mutCount=mutCount, 
	            WTCount=WTCount, 
	            kappa=kappa, 
	            no.iters.post.burn.in=no.iters.post.burn.in, 
	            no.iters=no.iters, 
	            no.iters.burn.in=no.iters.burn.in, 
	            subsamplenames=subsamplenames, 
	            resort.mutations=resort.mutations, 
	            plot.devs=plot.devs)
  
	all.consensus.trees = res$trees; 
	all.consensus.assignments = res$all.consensus.assignments
	likelihoods = res$likelihoods
	all.disaggregated.consensus.assignments = res$all.disaggregated
  all.likelihoods = res$all.likelihoods
  all.thetas = res$all.thetas

  print(paste("Finished EM in", as.numeric(Sys.time()-start,units="secs"), "seconds"))

  # There cannot possibly be no consensus trees. If this happens the random walk data is not right. Therefore exit.
  if (length(all.consensus.trees) == 0) {
    print("Did not find any consensus tree. Exit")
    q(save="no", status=1)
  }

  # No more trees to be plotted. Close the devices
  close.plotting.devices(no.subsamples, plot.devs)
	tree.sizes = sapply(1:length(all.consensus.trees),function(t,i){nrow(t[[i]])},t=all.consensus.trees)
  BIC = bic(likelihoods, no.subsamples, tree.sizes, no.muts)
  AIC = aic(likelihoods, no.subsamples, tree.sizes)
  # Calculate the DIC for each tree separately. Usually the DIC takes into account previous theta's, but that's not what we want here.
  DIC = array(NA, ncol(all.likelihoods))
  for (i in 1:ncol(all.likelihoods)) {
    print(dim(all.likelihoods))
    print(length(all.thetas))
    print(length(likelihoods))
    
    # Subset the thetas to contain just those of the current tree
    all.thetas.i = list()
    for (j in 1:no.subsamples) { all.thetas.i[[j]] = all.thetas[[j]][,i] }
    DIC[i] = dic(y=mutCount, n=WTCount+mutCount, kappa, all.likelihoods[,i],all.thetas.i)
  }

	best.BIC.index = which.min(BIC)
	print("likelihoods and BIC")
	print(cbind(likelihoods,BIC,AIC,DIC))
	print(paste("best BIC index=",best.BIC.index,sep=""))
	print("best BIC tree:")
	print(all.consensus.trees[[best.BIC.index]])

  print(paste("Finished GenerateConsensus in", as.numeric(Sys.time()-start,units="secs"), "seconds"))
	
	#save(all.consensus.trees,file=paste(samplename,"_",no.iters,"iters_",no.iters.burn.in,"_allConsensusTrees.RData",sep=""))
	if(!is.null(bin.indices)){
		return(list(all.consensus.trees = all.consensus.trees, all.consensus.assignments = all.consensus.assignments, all.disaggregated.consensus.assignments = all.disaggregated.consensus.assignments, all.likelihoods=all.likelihoods, likelihoods = likelihoods, BIC = BIC, AIC=AIC, DIC=DIC))
		#save(all.consensus.assignments,file=paste(samplename,"_",no.iters,"iters_",no.iters.burn.in,"_allBinnedConsensusAssignments.RData",sep=""))
		#save(all.disaggregated.consensus.assignments,file=paste(samplename,"_",no.iters,"iters_",no.iters.burn.in,"_allConsensusAssignments.RData",sep=""))		
	}else{
		#save(all.consensus.assignments,file=paste(samplename,"_",no.iters,"iters_",no.iters.burn.in,"_allConsensusAssignments.RData",sep=""))
		return(list(all.consensus.trees = all.consensus.trees, all.consensus.assignments = all.consensus.assignments, all.likelihoods=all.likelihoods, likelihoods = likelihoods, BIC = BIC, AIC=AIC, DIC=DIC))

	}	
}

get.mut.ass.strengths = function(no.muts, no.iters, no.iters.post.burn.in, node.assignments) {
  #
  # Calculates for each pair of mutations how often the pair is assigned to:
  #   - The same node (identity)
  #   - Nodes in parent-offspring relation (ancestor)
  #   - Nodes that are siblings (sibling)
  #
  calc.ancestor.strengths <- function(m,ancestor.strengths, node.assignments, no.iters.since.burnin, identity.strengths) {
    #
    # Calculates the strengths for iteration m of the MCMC algorithm
    #
    node.assignments.all = node.assignments[,no.iters.since.burnin]
    node.assignments.m = node.assignments[m,no.iters.since.burnin]
    
    temp.ancestor.or.identity.relationship = younger.direct.descendants(node.assignments.m,node.assignments.all)
    temp.ancestor.strengths = ancestor.strengths[m,] + (temp.ancestor.or.identity.relationship & (node.assignments.m != node.assignments.all))
    temp.identity.strengths = identity.strengths[m,] + (node.assignments.m == node.assignments.all)
    
    return(list(temp.ancestor.strengths, temp.ancestor.or.identity.relationship, temp.identity.strengths))
  }
  
  ancestor.strengths = array(0,c(no.muts,no.muts))
  sibling.strengths = array(0,c(no.muts,no.muts))
  identity.strengths = array(0,c(no.muts,no.muts))
  
  for(i in 1:no.iters.post.burn.in){ # for each iteration past burnin
    ancestor.or.identity.relationship = array(NA,c(no.muts,no.muts))
    
    # i+no.iters-no.iters.post.burn.in = no.iters-no.iters.post.burn.in is where the burn.in stops, i represents the iterations after that point
    res = sapply(1:no.muts, FUN=calc.ancestor.strengths, ancestor.strengths, node.assignments, i+no.iters-no.iters.post.burn.in, identity.strengths)
    ancestor.strengths = do.call(rbind,res[1,])
    ancestor.or.identity.relationship = do.call(rbind,res[2,])
    identity.strengths = do.call(rbind,res[3,])
    
    sibling.strengths = sibling.strengths + as.numeric(!ancestor.or.identity.relationship & !t(ancestor.or.identity.relationship))
  }
  
  return(list(ancestor.strengths=ancestor.strengths, sibling.strengths=sibling.strengths, identity.strengths=identity.strengths))
}

do_em = function(trees,node.assignments,ancestor.strengths, sibling.strengths, identity.strengths, bin.indices, consensus.assignments, current.agreement, mutCount, WTCount, kappa, no.iters.post.burn.in, no.iters, no.iters.burn.in, subsamplenames, resort.mutations, plot.devs) {
  
  start = Sys.time()
	pairwise.agreements = identity.strengths
  
  tree.number=1
	no.muts=nrow(mutCount)
  no.subsamples=length(subsamplenames)

  # Set up storage variables
	likelihoods = numeric(0)
  all.likelihoods = array(NA,c(no.muts,0))
  all.thetas = list()
  for (i in 1:no.subsamples) { all.thetas[[i]] = array(NA,c(no.muts,0)) }
	post.mean.deviance = numeric(0)
	all.consensus.trees = list()
	all.consensus.assignments = list()
	if(!is.null(bin.indices)){
		all.disaggregated.consensus.assignments = list()
	}	
	new.pairwise.agreements = list()
  
	node.added=T
	while(node.added){	
    # Perform EM to try out various layout extensions
    res = do_em_layout(no.muts, consensus.assignments, pairwise.agreements, new.pairwise.agreements, identity.strengths, ancestor.strengths, sibling.strengths)
    new.agreements = res$new.agreements
    new.nodes = res$new.nodes
    muts.to.move = res$muts.to.move
    new.pairwise.agreements = res$new.pairwise.agreements
    
    gc()
    
    print("Moving round1 done")
		
    #
    # pull out the best extension/agreement found above
    # if there is a higher agreement than we have currently
    #   insert the new node into the tree
    #   get optimised tree and node assignments and save them
    #   calculate the likelihood of the new tree
    unique.nodes = unique(consensus.assignments)
    
    best.node = which.max(new.agreements)
		if(new.agreements[best.node] > current.agreement){
			new.node = new.nodes[best.node]
			n = floor((best.node-1)/2) + 1
      # Check if node should be inserted above the current assignment ??
			above = (best.node %% 2 == 0)
			level = length(strsplit(new.node,":")[[1]])
			if(above){
				for(y in which(younger.direct.descendants(unique.nodes[n],unique.nodes))){
					spl = strsplit(unique.nodes[y],":")[[1]]
					if(new.node == unique.nodes[y]){
						consensus.assignments[consensus.assignments == unique.nodes[y]] = paste(paste(c(spl[1:level],1),collapse=":"),":",sep="")						
					}else{
						consensus.assignments[consensus.assignments == unique.nodes[y]] = paste(paste(c(spl[1:level],1,spl[(level+1):length(spl)]),collapse=":"),":",sep="")
					}
				}
			}	
			consensus.assignments[muts.to.move[[best.node]]] = new.node
			current.agreement = new.agreements[best.node]
			pairwise.agreements = new.pairwise.agreements[[best.node]]
			fractional.current.agreement = current.agreement/(no.muts*no.muts*no.iters.post.burn.in)

      temp = plotConsensusTree(consensus.assignments, samplename, subsamplenames, no.iters, no.iters.burn.in, node.assignments, trees, mutCount, WTCount, kappa, plot.devs, bin.indices, tree.number)
			new.consensus.tree = temp[[1]]
      new.consensus.ass = temp[[2]]
      
      tree.number = tree.number+1
			all.consensus.trees[[length(all.consensus.trees)+1]] = new.consensus.tree
			all.consensus.assignments[[length(all.consensus.assignments)+1]] = new.consensus.ass
			if(!is.null(bin.indices)){
				all.disaggregated.consensus.assignments[[length(all.disaggregated.consensus.assignments)+1]] = temp[[3]]
			}
			
      colnames = paste("theta.S", 1:no.subsamples, sep="")
      new.likelihood = log.f.of.y(mutCount, mutCount+WTCount, kappa, new.consensus.tree[new.consensus.ass,colnames])
      # Save the likelihood per mutation in order to calculate the DIC score later
      all.likelihoods = cbind(all.likelihoods, new.likelihood)
      likelihoods = c(likelihoods,sum(new.likelihood))
      for (i in 1:no.subsamples) {
        # Save the theta per subsample in a separate data.frame
        all.thetas[[i]] = cbind(all.thetas[[i]], new.consensus.tree[new.consensus.ass,paste("theta.S",i,sep='')])
      }
      
			#fix tree structure and shuffle mutation assignments
			if(resort.mutations){
        print("In resort mutations")
				unique.nodes = unique(consensus.assignments)

				saved.consensus.assignments = consensus.assignments
				mut.moved=T
				count=1
				#we may need to avoid infinite cycling by checking whether a set of node assignments has been repeated
				#assignments.tracker = matrix(NA, length(consensus.assignments), length(consensus.assignments))
				#assignments.tracker[,1] = consensus.assignments
        #no.assignments.tracked = 1
        
				while(mut.moved){
          if (! count %% 50) {
            print(paste("Round 2: ",count, sep=''))
            print(Sys.time()-start)
          }
          
					count=count+1
					mut.moved=F
					rand.inds = sample(no.muts)
					for(r in rand.inds){
						old.agreement = sum(pairwise.agreements[r,])
						curr.ass = consensus.assignments[r]
						possible.ass = unique.nodes[unique.nodes != curr.ass]
            
            # Get best new agreement
            new.agreements = as.array(sapply(1:length(possible.ass), FUN=get.new.agreements, r, possible.ass, unique.nodes, consensus.assignments, identity.strengths, ancestor.strengths, sibling.strengths))
						new.agreement = max(new.agreements)
						best.index = which.max(new.agreements)
            
            # Change the consensus assignment temporarily and check whether this new set of assignments was already 
            # tried. If it was already tried it should not be tried again
#             consensus.ass.temp = consensus.assignments
# 						consensus.ass.temp[r] = possible.ass[best.index]

            # Check if the new assignment when moving this mut was not seen before
            if(new.agreement > old.agreement) { #} & !is.assignment.known(assignments.tracker, consensus.ass.temp)){
							mut.moved=T
							res = get.desc.and.anc(possible.ass[best.index], unique.nodes)
              anc = res$anc; desc = res$desc
              consensus.assignments[r] = possible.ass[best.index]
# 							pairwise.agreements = move.mut2(r,possible.ass[best.index], unique.nodes, pairwise.agreements, consensus.assignments, identity.strengths, ancestor.strengths, sibling.strengths)
							pairwise.agreements = update.agreements(r, possible.ass[best.index], desc, anc, pairwise.agreements, consensus.assignments, identity.strengths, ancestor.strengths, sibling.strengths)
						}
					}
				}
        
        gc()

  			print("Moving round2 done")
				# Remove empty nodes
# 				save(file=paste("cullTree_",format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".RData", sep=""), consensus.assignments, all.consensus.trees)
				consensus.assignments = cullTree(consensus.assignments, all.consensus.trees)
        #if(ncol(assignments.tracker) == no.assignments.tracked) {
        #  ## The assignments tracker is getting full, extend it with more columns
        #  assignments.tracker = cbind(assignments.tracker, matrix(NA, length(consensus.assignments), length(consensus.assignments)))
        #}
        #assignments.tracker[,no.assignments.tracked+1] = consensus.assignments
				#no.assignments.tracked = no.assignments.tracked+1
        
				temp = plotConsensusTree(consensus.assignments, samplename, subsamplenames, no.iters, no.iters.burn.in, node.assignments, trees, mutCount, WTCount, kappa, plot.devs, bin.indices, tree.number)
        new.consensus.tree = temp[[1]]
        new.consensus.ass = temp[[2]]
        tree.number = tree.number+1
				
				all.consensus.trees[[length(all.consensus.trees)+1]] = new.consensus.tree
				all.consensus.assignments[[length(all.consensus.assignments)+1]] = temp[[2]]
				if(!is.null(bin.indices)){
					all.disaggregated.consensus.assignments[[length(all.disaggregated.consensus.assignments)+1]] = new.consensus.ass
				}
			
        new.likelihood = log.f.of.y(mutCount, mutCount+WTCount, kappa, new.consensus.tree[new.consensus.ass,colnames])
        # Save the likelihood per mutation in order to calculate the DIC score later
        all.likelihoods = cbind(all.likelihoods, new.likelihood)
				likelihoods = c(likelihoods,sum(new.likelihood))
        for (i in 1:no.subsamples) {
          # Save the theta per subsample in a separate data.frame
          all.thetas[[i]] = cbind(all.thetas[[i]], new.consensus.tree[new.consensus.ass,paste("theta.S",i,sep='')])
			  }
				print("finished re-sorting mutations")
			}
		}else{
			node.added = F
		}
	}
	if(is.null(bin.indices)) {
            all.disaggregated.consensus.assignments = list()
        }
	return(list(trees=all.consensus.trees, all.consensus.assignments=all.consensus.assignments, likelihoods=likelihoods, all.disaggregated=all.disaggregated.consensus.assignments, all.likelihoods=all.likelihoods, all.thetas=all.thetas))
}


do_em_layout = function(no.muts, consensus.assignments, pairwise.agreements, new.pairwise.agreements, identity.strengths, ancestor.strengths, sibling.strengths) {
  #
  # Searches for nodes that benefit from having a new node assigned either below or above it.
  # Mutations are then shuffled around to see whether the layout change has improved the agreements.
  # It returns a list that contains an array of new nodes and the associated agreements.
  #
  unique.nodes = unique(consensus.assignments)
  no.nodes = length(unique.nodes)
  
  new.agreements = array(NA,2*no.nodes)
  new.nodes = array(NA,2*no.nodes)
  muts.to.move = list()
  for(n in 1:no.nodes){
    # 
    # for each node
    #   add a new node below current node
    #   while a mutation was moved 
    #     move each mutation to the new node, or move it back when it was already assigned there
    #     calculate agreement
    #     if agreement is better, then move the mutation permanently
    #   
    #   then add new node above current node and perform the shuffling
    #
    #   save which mutations are to be moved
    #
    for(above in c(F,T)){
      new.pairwise.agreements[[2*(n-1)+above+1]] = pairwise.agreements
      new.consensus.assignments = consensus.assignments
      level = length(strsplit(unique.nodes[n],":")[[1]])
      if(above){
        # Add the new node above the current level
        res = add.node.above.curr.level(n,unique.nodes,new.consensus.assignments,level)
        new.node = res[[1]]
        new.consensus.assignments = res[[2]]
      }else{
        # Add new node below current level
        new.node = add.node.below.curr.level(n,unique.nodes,level)
      }
      new.nodes[2*(n-1)+above+1] = new.node
      new.unique.nodes = unique(c(new.consensus.assignments,new.node))
      
      #EM algorithm - iteratively move muts to the new node or back again
      mut.moved=T
      count=1
      
      saved.consensus.assignments = new.consensus.assignments
      #we may need to avoid infinite cycling by checking whether a set of node assignments has been repeated
      while(mut.moved){
        
        if (! count %% 50) {
          print(paste("Round 1: ",count, sep=''))
          print(Sys.time()-start)
        }
        
        count=count+1
        mut.moved=F
        rand.inds = sample(no.muts)
        for(r in rand.inds){
          # Check whether mutation was already assigned to the new node
          if(new.consensus.assignments[r]==new.node){
            # Move mut back to old node
            new.ass = saved.consensus.assignments[r]
          }else{
            # Save the new node as the node to be assigned to
            new.ass = new.node
          }
          res = get.desc.and.anc(new.ass, new.unique.nodes)
          desc = res$desc; anc = res$anc
          
          temp.ass = new.consensus.assignments
          temp.ass[r] = new.ass
          
          # Calculate new agreement by summing the scores for ancestor and siblings
          new.agreement = calc.new.agreement(r, temp.ass, new.ass, desc, anc, identity.strengths, ancestor.strengths, sibling.strengths)
          # Get assignment to same node for this mutation with all others
          old.agreement = sum(new.pairwise.agreements[[2*(n-1)+above+1]][r,])
          
          # When the new agreement is better, move the mutation
          if(new.agreement > old.agreement){
            mut.moved=T
            new.consensus.assignments[r] = new.ass
            new.pairwise.agreements[[2*(n-1)+above+1]] = update.agreements(r, new.ass, desc, anc, new.pairwise.agreements[[2*(n-1)+above+1]], new.consensus.assignments, identity.strengths, ancestor.strengths, sibling.strengths)
            old.agreement = new.agreement
          }
        }
      }
      new.agreements[2*(n-1)+above+1] = sum(new.pairwise.agreements[[2*(n-1)+above+1]])
      #don't move a whole node of mutations - not sure whether this should be allowed
      if(length(unique(new.consensus.assignments))<=no.nodes){
        new.agreements[2*(n-1)+above+1] = 0
      }			
      muts.to.move[[2*(n-1)+above+1]] = which(new.consensus.assignments == new.node)
    }
  }
  return(list(new.agreements=new.agreements, new.nodes=new.nodes, muts.to.move=muts.to.move, new.pairwise.agreements=new.pairwise.agreements))
}

add.node.above.curr.level = function(n1,unique.nodes1,new.consensus.assignments1,level1) {
  # Add the new node above the current level
  for(y in which(younger.direct.descendants(unique.nodes1[n1],unique.nodes1))){
    # move the mutations to the node where they were previously assigned for each direct descendant
    spl = strsplit(unique.nodes1[y],":")[[1]]
    if(unique.nodes1[n1] == unique.nodes1[y]){
      # new top node
      new.consensus.node = paste(paste(c(spl[1:level1],1),collapse=":"),":",sep="")						
    }else{
      # new inbetween node
      new.consensus.node = paste(paste(c(spl[1:level1],1,spl[(level1+1):length(spl)]),collapse=":"),":",sep="")
    }
    new.consensus.assignments1[new.consensus.assignments1 == unique.nodes1[y]] = new.consensus.node
  }
  return(list(unique.nodes1[n1],new.consensus.assignments1))
}

add.node.below.curr.level = function(n1,unique.nodes1,level1) {
  # Add new node below current level
  # Check whether there are any children
  no.children=0
  for(y in which(younger.direct.descendants(unique.nodes1[n1],unique.nodes1))){
    spl = strsplit(unique.nodes1[y],":")[[1]]
    if(length(spl)==level1+1){
      no.children = max(no.children,as.numeric(spl[level1+1]))
    }
  }
  # New node is no.children (current) + 1
  return(paste(unique.nodes1[n1],no.children+1,":",sep=""))
}


is.assignment.known <- function(ass.tracker, consensus.ass.temp) {
  res = sapply(1:ncol(ass.tracker), FUN=function(col, ass.tracker, consensus.ass.temp) { identical(ass.tracker[,col],consensus.ass.temp) }, ass.tracker, consensus.ass.temp)
  return(sum(res) > 0)
}

cullTree = function(consensus.assignments, all.consensus.trees) {
  #
	# cull tree (remove empty nodes), and reassign labels
  #
	temp.tree = all.consensus.trees[[length(all.consensus.trees)]]
	temp.tree$node = 1:nrow(temp.tree)
	temp.list <- cull.tree(temp.tree, consensus.assignments)
	mapping <- temp.list$mapping
	new.node.assignments = vector(length = length(consensus.assignments), mode = "character")
	for(i in 1:nrow(mapping)){
		new.node.assignments[consensus.assignments==mapping$old[i]] = mapping$new[i]
	}
	return(new.node.assignments)
}

update.agreements = function(r1, new.ass1, desc1, anc1, new.pairwise.agreements1, new.consensus.assignments1, identity.strengths1, ancestor.strengths1, sibling.strengths1) {
  #
  # This method reassigns the agreements from identity, sibling and ancestor strenghts after mutations have been reassigned.
  # The resulting matrix then contains pairwise aggreements, which are supposed to go up when the consensus tree becomes more
  # like the trees that are created during MCMC.
  #
  new.pairwise.agreements1[r1,] = NA
  new.pairwise.agreements1[,r1] = NA
  new.pairwise.agreements1[r1,new.consensus.assignments1==new.ass1] = identity.strengths1[r1,new.consensus.assignments1==new.ass1]
  new.pairwise.agreements1[new.consensus.assignments1==new.ass1,r1] = identity.strengths1[new.consensus.assignments1==new.ass1,r1]
  new.pairwise.agreements1[r1,new.consensus.assignments1 %in% desc1] = ancestor.strengths1[r1,new.consensus.assignments1 %in% desc1]
  new.pairwise.agreements1[new.consensus.assignments1 %in% desc1,r1] = ancestor.strengths1[r1,new.consensus.assignments1 %in% desc1]
  new.pairwise.agreements1[new.consensus.assignments1 %in% anc1,r1] = ancestor.strengths1[new.consensus.assignments1 %in% anc1,r1]
  new.pairwise.agreements1[r1,new.consensus.assignments1 %in% anc1] = ancestor.strengths1[new.consensus.assignments1 %in% anc1,r1]
  new.pairwise.agreements1[!(new.consensus.assignments1 %in% anc1) & !(new.consensus.assignments1 %in% desc1) & new.consensus.assignments1 != new.ass1,r1] = sibling.strengths1[!(new.consensus.assignments1 %in% anc1) & !(new.consensus.assignments1 %in% desc1) & new.consensus.assignments1 != new.ass1,r1]
  new.pairwise.agreements1[r1, !(new.consensus.assignments1 %in% anc1) & !(new.consensus.assignments1 %in% desc1) & new.consensus.assignments1 != new.ass1] = sibling.strengths1[r1, !(new.consensus.assignments1 %in% anc1) & !(new.consensus.assignments1 %in% desc1) & new.consensus.assignments1 != new.ass1]
  return(new.pairwise.agreements1)
}

calc.curr.agreement = function(consensus.assignments, identity.strengths, ancestor.strengths, sibling.strengths) {
  #
  # Calculates the current agreement of node assignments by looking at the strengths for each pair of nodes.
  #
  current.agreement = 0
  for(i in unique(consensus.assignments)){
    inds1 = which(consensus.assignments==i)
    for(j in unique(consensus.assignments)){
      inds2 = which(consensus.assignments==j)
      if(i==j){
        current.agreement = current.agreement + sum(identity.strengths[inds1,inds2])
      }else if(younger.direct.descendants(i,j)){
        current.agreement = current.agreement + sum(ancestor.strengths[inds1,inds2])
      }else if(younger.direct.descendants(j,i)){
        current.agreement = current.agreement + sum(ancestor.strengths[inds2,inds1])
      }else{
        current.agreement = current.agreement + sum(sibling.strengths[inds1,inds2])
      }		
    }
  }
  return(current.agreement)
}

get.new.agreements <- function(n, r, possible.ass, unique.nodes, consensus.assignments, identity.strengths, ancestor.strengths, sibling.strengths) {
  #
  # Returns the agreement if a possible assignment was to be performed
  #
  new.ass = possible.ass[n]
  res = get.desc.and.anc(new.ass, unique.nodes)
  desc = res$desc; anc = res$anc
  
  temp.ass = consensus.assignments
  temp.ass[r] = new.ass
  
  new.agreement = calc.new.agreement(r, temp.ass, new.ass, desc, anc, identity.strengths, ancestor.strengths, sibling.strengths)
  return(new.agreement)
}

calc.new.agreement = function(r, temp.ass, new.ass, desc, anc, identity.strengths, ancestor.strengths, sibling.strengths) {
  #
  # Calculates the new agreement when a new assignment is the temp assignment. This is used in order to test
  # whether a possible 
  #
  return(sum(identity.strengths[r,temp.ass==new.ass]) + 
    sum(ancestor.strengths[r,temp.ass %in% desc]) + 
    sum(ancestor.strengths[temp.ass %in% anc,r]) + 
    sum(sibling.strengths[!(temp.ass %in% anc) & !(temp.ass %in% desc) & temp.ass != new.ass,r]))
}

get.desc.and.anc = function(node, unique.nodes) {
  #
  # Returns the descendants and ancestor of the node
  #
  desc = unique.nodes[younger.direct.descendants(node,unique.nodes)]
  desc = desc[desc != node]
  anc = unique.nodes[ancestors(node,unique.nodes)]
  anc = anc[anc != node]

  return(list(desc=desc,anc=anc))
}

create.plotting.devices = function(samplename, no.iters, no.iters.burn.in, no.subsamples) {
  pdf(paste("histograms_",samplename,"_",no.iters,"iters_",no.iters.burn.in,"burnin.pdf",sep=""),height=4,width=4*no.subsamples)
  hist.device=dev.cur()
  par(mfrow=c(2,no.subsamples))	
  pdf(paste("densities_",samplename,"_",no.iters,"iters_",no.iters.burn.in,"burnin.pdf",sep=""),height=4,width=4*no.subsamples)
  density.device=dev.cur()	
  par(mfrow=c(2,no.subsamples))
  pdf(paste("consensus_trees_",samplename,"_",no.iters,"iters_",no.iters.burn.in,"burnin.pdf",sep=""),height=4,width=4*no.subsamples)
  tree.device=dev.cur()
  pdf(paste("consensus_trees_optimised_",samplename,"_",no.iters,"iters_",no.iters.burn.in,"burnin.pdf",sep=""),height=4,width=4*no.subsamples)
  optimised.tree.device=dev.cur()
  pdf(paste("consensus_trees_no_of_muts_",samplename,"_",no.iters,"iters_",no.iters.burn.in,"burnin.pdf",sep=""),height=4,width=4)
  tree.population.device=dev.cur()
  if(no.subsamples>1){
    pdf(paste("consensus_scatter_",samplename,"_",no.iters,"iters_",no.iters.burn.in,"burnin.pdf",sep=""),height=4,width=no.subsamples*(no.subsamples-1)*2)
    scatter.device=dev.cur()
    par(mfrow=c(1,no.subsamples*(no.subsamples-1)/2))
  } else {
    scatter.device = NA
  }
  return(list(hist.device=hist.device, density.device=density.device, tree.device=tree.device, optimised.tree.device=optimised.tree.device, tree.population.device=tree.population.device, scatter.device=scatter.device))
}

# close.plotting.devices = function(no.subsamples, hist.device, density.device, tree.device, optimised.tree.device, tree.population.device, scatter.device) {
close.plotting.devices = function(no.subsamples, plot.devs) {
  dev.off(which = plot.devs$hist.device)
  dev.off(which = plot.devs$density.device)
  dev.off(which = plot.devs$tree.device)
  dev.off(which = plot.devs$optimised.tree.device)
  dev.off(which = plot.devs$tree.population.device)
  if(no.subsamples>1){
    dev.off(which = plot.devs$scatter.device)
  }
}

