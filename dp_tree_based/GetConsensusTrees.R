source("Tree_based_DP_Gibbs_sampler.R")
# source("interconvertMutationBurdens.R")
source("PlotTreeWithIgraph.R")
source("DensityEstimator.R")
source("CullTree.R")

makeValidThetas<-function(tree){
	theta.cols = names(tree)[grep("theta",names(tree))]
	levels = sapply(1:nrow(tree),function(t,i){length(strsplit(rownames(t)[i],":")[[1]])},t=tree)
	print(levels)
	if(max(levels)>1){
		for(level in 2:max(levels)){
			level.list = rownames(tree)[levels==level]
			print(level)
			print(level.list)
			for(anc in unique(tree[level.list,"ancestor"])){
				node.list = level.list[tree$ancestor[match(level.list,rownames(tree))]==anc]
				print(tree$ancestor[match(level.list,rownames(tree))])
				print(anc)
				print(node.list)
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

plotConsensusTree<-function(consensus.assignments,samplename,subsamplenames,no.iters, no.iters.burn.in, node.assignments,trees,mutCount,WTCount,kappa,hist.device,density.device,tree.device,optimised.tree.device,tree.population.device,scatter.device,bin.indices=NULL,shrinkage.threshold = 0.1, tree.number=NA){
  start = Sys.time()
  
  no.iters.post.burn.in = no.iters-no.iters.burn.in
	no.subsamples = length(subsamplenames)
	no.muts = length(consensus.assignments)
	
	subclonal.fraction = array(NA,dim(mutCount))
	for(i in 1:no.subsamples){
		subclonal.fraction[,i] = mutCount[,i] / ((mutCount[,i]+WTCount[,i])*kappa[,i])
		subclonal.fraction[kappa[,i]==0,i] = NA
	}

	##############################Get Node Info############################################################
	temp.unique.nodes = unique(consensus.assignments)
	levels = sapply(1:length(temp.unique.nodes),function(t,i){length(strsplit(t[i],":")[[1]])},t=temp.unique.nodes)
	max.level = max(levels)
	unique.nodes=character(0)
	for(l in 1:max.level){
		nodes = temp.unique.nodes[levels==l]
		unique.nodes=c(unique.nodes,sort(nodes))
	}

	consensus.thetas=list()
	no.nodes = length(unique.nodes)
	is.single.point=vector(mode="logical",length=no.nodes)
  print(Sys.time()-start)
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
  print(Sys.time()-start)
	##############################PLOTTING############################################################		
	dev.set(which = hist.device)
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

	dev.set(which = density.device)
	highest.density.thetas = array(NA,c(no.subsamples,no.nodes))

	# Estimate theta densities
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

	dev.set(which = tree.device)
	consensus.tree = as.data.frame(cbind(unique.nodes,t(highest.density.thetas),ancs,NA),stringsAsFactors=F)
	names(consensus.tree) = c("label",paste("theta.S",1:no.subsamples,sep=""),"ancestor","annotation")
	for(n in 1:no.subsamples){
		consensus.tree[,n+1] = as.numeric(consensus.tree[,n+1])
	}
	row.names(consensus.tree) = consensus.tree$label
	plotTree(consensus.tree,main=paste(samplename,subsamplenames,sep=""),font.size=1.25)
  print(Sys.time()-start)
	print("finished plotting tree")

	#check how many mutations were aggregated in the binning
	if(!is.null(bin.indices)){
		print(paste("no.muts=",no.muts,sep=""))
		print(paste("length(bin.indices)=",length(bin.indices),sep=""))
		mut.table = vector(mode="numeric",length = no.nodes)
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

	#091013 - average theta, weighted by depth, which may be a better starting point for the whole.tree sampler
	#this tree should be recorded (and plotted)
	for (k in 1:no.subsamples) {
		for(n in 1:no.nodes){
			node.inds = which(consensus.assignments==unique.nodes[n])
			weights = mutCount[node.inds,k] + WTCount[node.inds,k]
			weights = weights/sum(weights)
			average.theta = mean(sapply((no.iters.burn.in+1):no.iters,function(i){sum(trees[[i]][node.assignments[node.inds,i],paste("theta.S",k,sep="")]*weights)}))				
			consensus.tree[unique.nodes[n],k+1] = average.theta
		}
	}
  print(Sys.time()-start)
	print("average tree:")
	print(consensus.tree)
	
	#adjust thetas to make a valid tree, otherwise whole.tree.slice.sampler will be unable to find a solution
	consensus.tree = makeValidThetas(consensus.tree)
	print("adjusted average tree:")
	print(consensus.tree)
	
	#resample whole tree
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
	dev.set(which=optimised.tree.device)
	plotTree(consensus.tree,main=paste(samplename,subsamplenames,sep=""),font.size=1.25)
  print(Sys.time()-start)
	print("finished plotting optimised tree")

	dev.set(which=tree.population.device)
	population.tree=consensus.tree[,c("label","theta.S1","ancestor","annotation")]
	for(n in 1:no.nodes){
		#population.tree[unique.nodes[n],"theta.S1"] = mut.table[unique.nodes[n],1]
		population.tree[unique.nodes[n],"theta.S1"] = mut.table[unique.nodes[n]]
	}
	plotTree(population.tree,main=samplename,font.size=1.25,plotAsPercentage=F)
  print(Sys.time()-start)
	print("population tree plotted")

	if(no.subsamples>1){
		#its hard to distinguish more than 8 different colours
		max.cols = 8
		cols = rainbow(min(max.cols,no.nodes))
		dev.set(which = scatter.device)
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
	#returned optimised tree and node assignments	
	if(!is.null(bin.indices)){
		return(list(consensus.tree,consensus.assignments,disaggregated.consensus.assignments))
	}else{
		return(list(consensus.tree,consensus.assignments))
	}
}
	
calc.subclonal.fraction <- function(i,mutCount1, WTCount1, kappa1) {
  fraction = mutCount1[,i] / ((mutCount1[,i]+WTCount1[,i])*kappa1[,i])
  fraction[kappa1[,i]==0] = NA
  return(fraction)
}

# calc.ancestor.strengths <- function(m, node.assignments1, no.iters.since.burnin) {
#   node.assignments.all = node.assignments1[,no.iters.since.burnin]
#   node.assignments.m = node.assignments1[m,no.iters.since.burnin]
#   
#   ancestor.strengths.local = (younger.direct.descendants(node.assignments.m,node.assignments.all) & (node.assignments.m != node.assignments.all))
#   ancestor.or.identity.relationship.local = younger.direct.descendants(node.assignments.m,node.assignments.all)
#   identity.strengths.local = (node.assignments.m == node.assignments.all)
#   return(list(ancestor.strengths.local, ancestor.or.identity.relationship.local, identity.strengths.local))
# }

calc.ancestor.strengths <- function(m,ancestor.strengths, node.assignments, no.iters.since.burnin, identity.strengths) {
	node.assignments.all = node.assignments[,no.iters.since.burnin]
	node.assignments.m = node.assignments[m,no.iters.since.burnin]

	temp.ancestor.or.identity.relationship = younger.direct.descendants(node.assignments.m,node.assignments.all)
	temp.ancestor.strengths = ancestor.strengths[m,] + (temp.ancestor.or.identity.relationship & (node.assignments.m != node.assignments.all))
	temp.identity.strengths = identity.strengths[m,] + (node.assignments.m == node.assignments.all)
	return(list(temp.ancestor.strengths, temp.ancestor.or.identity.relationship, temp.identity.strengths))
}

GetConsensusTrees<-function(trees, node.assignments, mutCount, WTCount, kappa = array(0.5,dim(mutCount)), samplename="sample", subsamplenames = 1:ncol(mutCount), no.iters = 1250, no.iters.burn.in = 250, resort.mutations = T, shrinkage.threshold = 0.1, init.alpha = 0.01, outdir = getwd(),bin.indices = NULL){
  start = Sys.time()
  
  setwd(outdir)
	
	tree.number=1
	
	likelihoods = numeric(0)
	
	no.subsamples = ncol(mutCount)
	no.muts = nrow(mutCount)
 	subclonal.fraction = array(NA,dim(mutCount))
 	#for(i in 1:no.subsamples){
 	#	subclonal.fraction[,i] = mutCount[,i] / ((mutCount[,i]+WTCount[,i])*kappa[,i])
 	#	subclonal.fraction[kappa[,i]==0,i] = NA
 	#}
#  subclonal.fraction = sapply(1:no.subsamples, FUN=calc.subclonal.fraction, mutCount, WTCount, kappa)

	#subclonal.fraction2 = array(NA,dim(mutCount))
	subclonal.fraction = sapply(1:no.subsamples, FUN=function(i,mutCount,WTCount,kappa) { mutCount[,i] / ((mutCount[,i]+WTCount[,i])*kappa[,i]) }, mutCount,WTCount,kappa)
	subclonal.fraction[kappa == 0] = NA

	#print(head(subclonal.fraction))
	#print(head(subclonal.fraction2))

	no.iters.post.burn.in = no.iters-no.iters.burn.in

	ancestor.strengths = array(0,c(no.muts,no.muts))
	sibling.strengths = array(0,c(no.muts,no.muts))
	identity.strengths = array(0,c(no.muts,no.muts))

	#it would be faster to use apply
	for(i in 1:no.iters.post.burn.in){
		ancestor.or.identity.relationship = array(NA,c(no.muts,no.muts))
	#	for(m in 1:no.muts){
	#		ancestor.strengths[m,] = ancestor.strengths[m,] + (younger.direct.descendants(node.assignments[m,i+no.iters-no.iters.post.burn.in],node.assignments[,i+no.iters-no.iters.post.burn.in]) & (node.assignments[m,i+no.iters-no.iters.post.burn.in] != node.assignments[,i+no.iters-no.iters.post.burn.in]))
	#		ancestor.or.identity.relationship[m,] = younger.direct.descendants(node.assignments[m,i+no.iters-no.iters.post.burn.in],node.assignments[,i+no.iters-no.iters.post.burn.in])
	#	  identity.strengths[m,] = identity.strengths[m,] + (node.assignments[m,i+no.iters-no.iters.post.burn.in] == node.assignments[,i+no.iters-no.iters.post.burn.in])
	#	}

		res = sapply(1:no.muts, FUN=calc.ancestor.strengths, ancestor.strengths, node.assignments, i+no.iters-no.iters.post.burn.in, identity.strengths)
		ancestor.strengths = t(do.call(rbind,res[1,]))
		ancestor.or.identity.relationship = t(do.call(rbind,res[2,]))
		identity.strengths = t(do.call(rbind,res[3,]))
		
		sibling.strengths = sibling.strengths + as.numeric(!ancestor.or.identity.relationship & !t(ancestor.or.identity.relationship))
	}

  print(Sys.time()-start) # 3.455465 secs

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
	}


	consensus.assignments = rep("M:",no.muts)
	no.nodes=1
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
	fractional.current.agreement = current.agreement/(no.muts*no.muts*no.iters.post.burn.in)
  print(Sys.time()-start) # 3.525014 secs
	print(paste("M:",current.agreement,fractional.current.agreement))
	
	#initialise all mutations in one node, so the number of pairwise agreements is just the number of times a pair of mutations appear in the same node
	pairwise.agreements = identity.strengths

	all.consensus.trees = list()
	all.consensus.assignments = list()
	if(!is.null(bin.indices)){
		all.disaggregated.consensus.assignments = list()
	}	
	new.pairwise.agreements = list()
	node.added=T
	while(node.added){	
    print("Node added")
		unique.nodes = unique(consensus.assignments)
		no.nodes = length(unique.nodes)
		new.agreements = array(NA,2*no.nodes)
		new.nodes = array(NA,2*no.nodes)
		muts.to.move = list()
		for(n in 1:no.nodes){
			for(above in c(F,T)){
				new.pairwise.agreements[[2*(n-1)+above+1]] = pairwise.agreements
				new.consensus.assignments = consensus.assignments
				level = length(strsplit(unique.nodes[n],":")[[1]])
				if(above){
					for(y in which(younger.direct.descendants(unique.nodes[n],unique.nodes))){
						spl = strsplit(unique.nodes[y],":")[[1]]
						if(unique.nodes[n] == unique.nodes[y]){
							new.consensus.assignments[new.consensus.assignments == unique.nodes[y]] = paste(paste(c(spl[1:level],1),collapse=":"),":",sep="")						
						}else{
							new.consensus.assignments[new.consensus.assignments == unique.nodes[y]] = paste(paste(c(spl[1:level],1,spl[(level+1):length(spl)]),collapse=":"),":",sep="")
						}
					}
					new.node = unique.nodes[n]
				}else{
					spl = strsplit(unique.nodes,":")
					no.children=0
					for(y in which(younger.direct.descendants(unique.nodes[n],unique.nodes))){
						spl = strsplit(unique.nodes[y],":")[[1]]
						if(length(spl)==level+1){
							no.children = max(no.children,as.numeric(spl[level+1]))
						}
					}
					new.node = paste(unique.nodes[n],no.children+1,":",sep="")
				}
				new.nodes[2*(n-1)+above+1] = new.node
				new.unique.nodes = unique(c(new.consensus.assignments,new.node))

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
						old.agreement = sum(new.pairwise.agreements[[2*(n-1)+above+1]][r,])
						if(new.consensus.assignments[r]==new.node){
							new.ass = saved.consensus.assignments[r]
						}else{
							new.ass = new.node
						}
						desc = new.unique.nodes[younger.direct.descendants(new.ass,new.unique.nodes)]
						desc = desc[desc != new.ass]
						anc = new.unique.nodes[ancestors(new.ass,new.unique.nodes)]
						anc = anc[anc != new.ass]
						temp.ass = new.consensus.assignments
						temp.ass[r] = new.ass
						new.agreement = sum(identity.strengths[r,temp.ass==new.ass]) + sum(ancestor.strengths[r,temp.ass %in% desc]) + sum(ancestor.strengths[temp.ass %in% anc,r]) + sum(sibling.strengths[!(temp.ass %in% anc) & !(temp.ass %in% desc) & temp.ass != new.ass,r])
						if(new.agreement > old.agreement){
							mut.moved=T
							new.consensus.assignments[r] = new.ass
							new.pairwise.agreements[[2*(n-1)+above+1]][r,]=NA
							new.pairwise.agreements[[2*(n-1)+above+1]][,r]=NA
							new.pairwise.agreements[[2*(n-1)+above+1]][r,new.consensus.assignments==new.ass] = identity.strengths[r,new.consensus.assignments==new.ass]
							new.pairwise.agreements[[2*(n-1)+above+1]][new.consensus.assignments==new.ass,r] = identity.strengths[new.consensus.assignments==new.ass,r]
							new.pairwise.agreements[[2*(n-1)+above+1]][r,new.consensus.assignments %in% desc] = ancestor.strengths[r,new.consensus.assignments %in% desc]
							new.pairwise.agreements[[2*(n-1)+above+1]][new.consensus.assignments %in% desc,r] = ancestor.strengths[r,new.consensus.assignments %in% desc]
							new.pairwise.agreements[[2*(n-1)+above+1]][new.consensus.assignments %in% anc,r] = ancestor.strengths[new.consensus.assignments %in% anc,r]
							new.pairwise.agreements[[2*(n-1)+above+1]][r,new.consensus.assignments %in% anc] = ancestor.strengths[new.consensus.assignments %in% anc,r]
							new.pairwise.agreements[[2*(n-1)+above+1]][!(new.consensus.assignments %in% anc) & !(new.consensus.assignments %in% desc) & new.consensus.assignments != new.ass,r] = sibling.strengths[!(new.consensus.assignments %in% anc) & !(new.consensus.assignments %in% desc) & new.consensus.assignments != new.ass,r]
							new.pairwise.agreements[[2*(n-1)+above+1]][r, !(new.consensus.assignments %in% anc) & !(new.consensus.assignments %in% desc) & new.consensus.assignments != new.ass] = sibling.strengths[r, !(new.consensus.assignments %in% anc) & !(new.consensus.assignments %in% desc) & new.consensus.assignments != new.ass]
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
    print("Moving round1 done")
		print(Sys.time()-start) # 4.236197 secs
		best.node = which.max(new.agreements)
		if(new.agreements[best.node] > current.agreement){
			new.node = new.nodes[best.node]
			n = floor((best.node-1)/2) + 1
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
			print(paste(new.node,current.agreement,fractional.current.agreement))

			temp = plotConsensusTree(consensus.assignments,samplename,subsamplenames,no.iters, no.iters.burn.in,node.assignments,trees,mutCount,WTCount,kappa,hist.device,density.device,tree.device,optimised.tree.device,tree.population.device,scatter.device,bin.indices,tree.number)
			tree.number = tree.number+1
			all.consensus.trees[[length(all.consensus.trees)+1]] = temp[[1]]
			all.consensus.assignments[[length(all.consensus.assignments)+1]] = temp[[2]]
			if(!is.null(bin.indices)){
				all.disaggregated.consensus.assignments[[length(all.disaggregated.consensus.assignments)+1]] = temp[[3]]
			}
			
			new.likelihood = 0
			colnames = paste("theta.S", 1:no.subsamples, sep="")
			for(i in 1:no.muts){
				lfoy = log.f.of.y(mutCount[i,], mutCount[i,] + WTCount[i,], kappa[i,], temp[[1]][temp[[2]][i],colnames])
				if(!is.nan(lfoy)){
					new.likelihood <- new.likelihood + lfoy
				}
			}
			likelihoods = c(likelihoods,new.likelihood)
			
			#fix tree structure and shuffle mutation assignments
			if(resort.mutations){
				unique.nodes = unique(consensus.assignments)

				saved.consensus.assignments = consensus.assignments
				mut.moved=T
				count=1
				#we may need to avoid infinite cycling by checking whether a set of node assignments has been repeated
				while(mut.moved){
					count=count+1
					mut.moved=F
					rand.inds = sample(no.muts)
					for(r in rand.inds){
						old.agreement = sum(pairwise.agreements[r,])
						curr.ass = consensus.assignments[r]
						possible.ass = unique.nodes[unique.nodes != curr.ass]
						new.agreements = array(NA,length(possible.ass))
						for(n in 1:length(possible.ass)){
							new.ass = possible.ass[n]
							desc = unique.nodes[younger.direct.descendants(new.ass,unique.nodes)]
							desc = desc[desc != new.ass]
							anc = unique.nodes[ancestors(new.ass,unique.nodes)]
							anc = anc[anc != new.ass]
							temp.ass = consensus.assignments
							temp.ass[r] = new.ass
							new.agreements[n] = sum(identity.strengths[r,temp.ass==new.ass]) + sum(ancestor.strengths[r,temp.ass %in% desc]) + sum(ancestor.strengths[temp.ass %in% anc,r]) + sum(sibling.strengths[!(temp.ass %in% anc) & !(temp.ass %in% desc) & temp.ass != new.ass,r])
						}

						new.agreement = max(new.agreements)
						best.index = which.max(new.agreements)
						if(new.agreement > old.agreement){
							mut.moved=T
							new.ass = possible.ass[best.index]
							desc = unique.nodes[younger.direct.descendants(new.ass,unique.nodes)]
							desc = desc[desc != new.ass]
							anc = unique.nodes[ancestors(new.ass,unique.nodes)]
							anc = anc[anc != new.ass]			
							consensus.assignments[r] = new.ass
							pairwise.agreements[r,]=NA
							pairwise.agreements[,r]=NA
							pairwise.agreements[r,consensus.assignments==new.ass] = identity.strengths[r,consensus.assignments==new.ass]
							pairwise.agreements[consensus.assignments==new.ass,r] = identity.strengths[consensus.assignments==new.ass,r]
							pairwise.agreements[r,consensus.assignments %in% desc] = ancestor.strengths[r,consensus.assignments %in% desc]
							pairwise.agreements[consensus.assignments %in% desc,r] = ancestor.strengths[r,consensus.assignments %in% desc]
							pairwise.agreements[consensus.assignments %in% anc,r] = ancestor.strengths[consensus.assignments %in% anc,r]
							pairwise.agreements[r,consensus.assignments %in% anc] = ancestor.strengths[consensus.assignments %in% anc,r]
							pairwise.agreements[!(consensus.assignments %in% anc) & !(consensus.assignments %in% desc) & consensus.assignments != new.ass,r] = sibling.strengths[!(consensus.assignments %in% anc) & !(consensus.assignments %in% desc) & consensus.assignments != new.ass,r]
							pairwise.agreements[r, !(consensus.assignments %in% anc) & !(consensus.assignments %in% desc) & consensus.assignments != new.ass] = sibling.strengths[r, !(consensus.assignments %in% anc) & !(consensus.assignments %in% desc) & consensus.assignments != new.ass]
						}
					}
				}
        print("Moving round2 done")
				print(Sys.time()-start)
				#cull tree (remove empty nodes), and reassign labels
				temp.tree = all.consensus.trees[[length(all.consensus.trees)]]
				temp.tree$node = 1:nrow(temp.tree)
				temp.list <- cull.tree(temp.tree, consensus.assignments)
				mapping <- temp.list$mapping
				new.node.assignments = vector(length = length(consensus.assignments), mode = "character")
				for(i in 1:nrow(mapping)){
					new.node.assignments[consensus.assignments==mapping$old[i]] = mapping$new[i]
				}
				consensus.assignments = new.node.assignments

				temp = plotConsensusTree(consensus.assignments,samplename,subsamplenames,no.iters, no.iters.burn.in,node.assignments,trees,mutCount,WTCount,kappa,hist.device,density.device,tree.device,optimised.tree.device,tree.population.device,scatter.device,bin.indices,tree.number)
				tree.number = tree.number+1
				
				all.consensus.trees[[length(all.consensus.trees)+1]] = temp[[1]]
				all.consensus.assignments[[length(all.consensus.assignments)+1]] = temp[[2]]
				if(!is.null(bin.indices)){
					all.disaggregated.consensus.assignments[[length(all.disaggregated.consensus.assignments)+1]] = temp[[3]]
				}
			
				new.likelihood = 0
				colnames = paste("theta.S", 1:no.subsamples, sep="")
				for(i in 1:no.muts){		
					lfoy = log.f.of.y(mutCount[i,], mutCount[i,] + WTCount[i,], kappa[i,], temp[[1]][temp[[2]][i],colnames])					
					if(!is.nan(lfoy)){
						new.likelihood <- new.likelihood + lfoy
					}					
				}
				likelihoods = c(likelihoods,new.likelihood)
			
				print("finished re-sorting mutations")
			}

		}else{
			node.added = F
		}
	}
  print("EM done")
  print(Sys.time()-start)

	dev.off(which = hist.device)
	dev.off(which = density.device)
	dev.off(which = tree.device)
	dev.off(which = optimised.tree.device)
	dev.off(which = tree.population.device)
	if(no.subsamples>1){
		dev.off(which = scatter.device)
	}

	tree.sizes = sapply(1:length(all.consensus.trees),function(t,i){nrow(t[[i]])},t=all.consensus.trees)
	BIC = -2 * likelihoods + 2 * no.subsamples * tree.sizes * log(no.muts)	
	best.BIC.index = which.min(BIC)
	print("likelihoods and BIC")
	print(cbind(likelihoods,BIC))
	print(paste("best BIC index=",best.BIC.index,sep=""))
	print("best BIC tree:")
	print(all.consensus.trees[[best.BIC.index]])

  print("GenerateConsensus done")
  print(Sys.time()-start)
	
	#png(paste("log_likelihood_plot_",samplename,"_consensus_",no.iters,"iters.png",sep=""),width=1000)
	#par(cex.lab=3,cex.axis=3,lwd=3,mar=c(7,7,5,2))
	#plot(1:no.iters,likelihoods,type="l",col="red",xlab="MCMC iteration", ylab="log-likelihood",main=samplename)
	#dev.off()
	#png(paste("BIC_plot_",samplename,"_consensus_",no.iters,"iters.png",sep=""),width=1000)
	#par(cex.lab=3,cex.axis=3,lwd=3,mar=c(7,7,5,2))
	#plot(1:no.iters,BIC,type="l",col="red",xlab="MCMC iteration", ylab="BIC",main=samplename)
	#dev.off()
	
	#save(all.consensus.trees,file=paste(samplename,"_",no.iters,"iters_",no.iters.burn.in,"_allConsensusTrees.RData",sep=""))
	if(!is.null(bin.indices)){
		return(list(all.consensus.trees = all.consensus.trees, all.consensus.assignments = all.consensus.assignments, all.disaggregated.consensus.assignments = all.disaggregated.consensus.assignments, likelihoods = likelihoods, BIC = BIC))
		#save(all.consensus.assignments,file=paste(samplename,"_",no.iters,"iters_",no.iters.burn.in,"_allBinnedConsensusAssignments.RData",sep=""))
		#save(all.disaggregated.consensus.assignments,file=paste(samplename,"_",no.iters,"iters_",no.iters.burn.in,"_allConsensusAssignments.RData",sep=""))		
	}else{
		#save(all.consensus.assignments,file=paste(samplename,"_",no.iters,"iters_",no.iters.burn.in,"_allConsensusAssignments.RData",sep=""))
		return(list(all.consensus.trees = all.consensus.trees, all.consensus.assignments = all.consensus.assignments, likelihoods = likelihoods, BIC = BIC))

	}	
}
