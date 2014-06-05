library(LIM)
source("CullTree.R")
source("RemoveBranch.R")

younger.siblings <- function(curr.node, other.nodes) {
	# Function to test whether other.nodes are younger siblings of curr.node
	# Returns vector of logical values
	split.curr.node <- unlist(strsplit(curr.node, ":"))
	parent <- paste(paste(split.curr.node[1:(length(split.curr.node)-1)], collapse=":"), ":", sep="")
	family.position <- sapply(other.nodes, FUN=function(x) {if (x == "M:") {c(0,1)} else {a <- unlist(strsplit(x, ":")); c(as.double(a[length(a)]), length(a)) }})
	return(grepl(parent, other.nodes) & family.position[1,] > as.double(split.curr.node[length(split.curr.node)]) & family.position[2,] == length(split.curr.node))
}

older.siblings <- function(curr.node, other.nodes) {
	# Function to test whether other.nodes are younger siblings of curr.node
	# Returns vector of logical values
	split.curr.node <- unlist(strsplit(curr.node, ":"))
	parent <- paste(paste(split.curr.node[1:(length(split.curr.node)-1)], collapse=":"), ":", sep="")
	family.position <- sapply(other.nodes, FUN=function(x) {if (x == "M:") {c(0,1)} else {a <- unlist(strsplit(x, ":")); c(as.double(a[length(a)]), length(a)) }})
	return(grepl(parent, other.nodes) & family.position[1,] < as.double(split.curr.node[length(split.curr.node)]) & family.position[2,] == length(split.curr.node))
}

ancestors <- function(curr.node, other.nodes) {
	# Function to test whether other.nodes are ancestors of curr.node
	# Note: curr.node is counted as an ancestor of itself
	# Returns vector of logical values
	no.nodes = length(other.nodes)
	is.ancestor = vector(mode="logical",length=no.nodes)
	for(i in 1:no.nodes){
		is.ancestor[i] = grepl(other.nodes[i],curr.node)
	}
	return(is.ancestor)
}

younger.descendants <- function(curr.node, other.nodes) {
	# Function to test whether other.nodes are younger siblings, direct descendants or descendants of younger siblings of curr.node
	# Returns vector of logical values, including TRUE for where curr.node == other.nodes[i]
	split.curr.node <- unlist(strsplit(curr.node, ":"))
	parent <- paste(paste(split.curr.node[1:(length(split.curr.node)-1)], collapse=":"), ":", sep="")
	family.position.curr.depth <- sapply(other.nodes, FUN=function(x,y) {if (x == "M:") {0} else {a <- unlist(strsplit(x, ":")); max(0,as.double(a[y]),na.rm=TRUE) }}, y=length(split.curr.node))
	return(grepl(parent, other.nodes) & family.position.curr.depth >= as.double(split.curr.node[length(split.curr.node)]))
}

#this treats older and younger siblings symetrically by excluding younger siblings and descendants of younger siblings
younger.direct.descendants <- function(curr.node, other.nodes) {
	# Function to test whether other.nodes are direct descendants of curr.node
	# Returns vector of logical values, including TRUE for where curr.node == other.nodes[i]
	return(grepl(curr.node, other.nodes))
}

log.f.of.y <- function(y1, n1, kappa1, x) {
	#x=1 and kappa=1 causes problems
	x[x>0.999 & kappa1==1] = 0.999
	#allow kappa = 0, for mutations on deleted chromosomes
	non.zero.inds = which(kappa1!=0)
	if(length(non.zero.inds)>0){
		return(sum(lchoose(n1[non.zero.inds], y1[non.zero.inds]) + y1[non.zero.inds] * log(kappa1[non.zero.inds]*x[non.zero.inds]) + (n1[non.zero.inds]-y1[non.zero.inds]) * log(1-kappa1[non.zero.inds]*x[non.zero.inds])) * length(y1)/length(non.zero.inds))
	}else{
		print("WARNING. ALL KAPPAS ARE ZERO. MUTATION ABSENT FROM ALL SAMPLES")
		print(y1)
		print(n1)
		print(kappa1)
		print(x)
		#return(NA)
		return(NaN)
	}
}

#library(compiler)
#log.f.of.y = cmpfun(log.f.of.y)

get.conflicts <- function(index, conflict.array, node.assignments, whole.tree){
	all.nodes = row.names(whole.tree)
	penalties = vector(mode="numeric",length=length(all.nodes))	
	penalties = rep(1,length(all.nodes))
	for(i in 1:length(node.assignments)){
		if(conflict.array[index,i]!=1){
			#penalise all ancestors of conflicting node, i.e. don't put conflicting variants in the same branch
			ancs = which(ancestors(node.assignments[i],all.nodes))
			spl = strsplit(node.assignments[i],":")[[1]]
			for(anc in ancs){
				if(all.nodes[anc] == node.assignments[i]){
					penalties[anc] = penalties[anc] * conflict.array[index,i]
				}else{
					#incorporate fraction of mutations that are expected to conflict
					expected.fraction=1
					anc.spl = strsplit(all.nodes[anc],":")[[1]]
					for(a in (length(anc.spl)+1):length(spl)){
						expected.fraction = expected.fraction * whole.tree$phi[all.nodes==paste(paste(spl[1:a],collapse=":"),":",sep="")]
					}
					if(is.na(expected.fraction)){
						print("ERROR IN EXPECTED FRACTION")
						print(paste(spl,sep=","))
						print(length(spl))
						print(length(anc.spl))
					}
					penalties[anc] = penalties[anc] * conflict.array[index,i] ^ expected.fraction
				}
			}
		}
	}
	names(penalties) = all.nodes
	return(penalties)
}

update.tree.after.removal <- function(curr.tree, removed.indices, mapping, num.muts, y, n, kappa) {
  new.node.assignments = vector(length = num.muts, mode = "character")
  for(i in 1:nrow(mapping)){
    new.node.assignments[new.assignments==mapping$old[i]] = mapping$new[i]
  }
  new.assignments = new.node.assignments
  
  rand.inds = sample(removed.indices)
  # Reassign mutations that were on the removed branch
  #length should always be greater than zero, because tree has been culled. Why is length sometimes zero???
  if(length(rand.inds)>0){
    for (r in 1:length(rand.inds)) {
      log.probs = vector(mode="numeric",length=nrow(curr.tree))
      for (x in 1:nrow(curr.tree)) {
        thetas = unlist(curr.tree[x, grepl("theta", names(curr.tree))])	
        log.probs[x] = log.f.of.y(y[rand.inds[r],],n[rand.inds[r],],kappa[rand.inds[r],],thetas)
      }
      log.probs = log.probs-max(log.probs,na.rm=T)
      probs = exp(log.probs)
      probs[is.nan(probs)] = 0
      cum.probs = cumsum(probs)
      ass.ind = sum(runif(1,0,max(cum.probs))>cum.probs)+1
      new.assignments[rand.inds[r]] = curr.tree$label[ass.ind]
    }
    
  }else{
    print("WARNING! No removed indices after remove.node.")
  }
    
  return(new.assignments)
}

bic <- function(likelihood, num.samples, num.trees, num.muts) {
  return(-2 * likelihood + 2 * num.samples * num.trees * log(num.muts))
#   return(-2 * complete.likelihood[1] + 2 * num.samples * nrow(trees.n[[1]]) * log(num.muts))
}

tree.struct.dirichlet.gibbs <- function(y, n, kappa, iter=1000, d=1, plot.lambda=FALSE, allowed.ranges=c(min.lambda = 0.05, max.lambda = 0.8, min.alpha = 1E-6, max.alpha=50, min.gamma=0.1, max.gamma=10), shrinkage.threshold = 0.1, init.alpha=0.01, init.lambda=0.5, init.gamma = 0.2, remove.node.frequency = NA, remove.branch.frequency = NA, conflict.array = array(1,c(nrow(y),nrow(y))), parallel=FALSE) {
	# y is a matrix of the number of reads reporting each variant - rows represent each mutation; columns each sample
	# N is a matrix of the number of reads in total across the base in question (in the same order as y obviously!)
	# kappa is a matrix of the expected proportion of reads reporting the variant for a fully clonal mutation at single copy
	# alpha0 is the initial branching parameter - ie alpha(x) = (lambda^x) * alpha0 for a tree of depth x
	# iter is the number of iterations of the Gibbs sampler
	# d is the scale parameter for the beta proposal distribution
	# The theta.S variables are the estimates of the fraction of tumour cells for each node in that sample
	# The shrinkage threshold is the threshold by which to shrink the constraints during the slice sampling for thetas
	
	if (class(y) == "numeric") {y <- as.matrix(y, ncol=1)}
	if (class(n) == "numeric") {n <- as.matrix(n, ncol=1)}
	if (class(kappa) == "numeric") {kappa <- as.matrix(kappa, ncol=1)}
		
	num.muts <- dim(y)[1]
	num.samples <- dim(y)[2]
	
	# Set up data formats for recording iterations
	trees.n <- vector("list", length=iter)  # The set of trees for each iteration
	node.assignments <- matrix(NA, nrow=num.muts, ncol=iter)
	lambda <- alpha0 <- gamma <- rep(NA, iter)  # The hyperparameter for the decay in tree depth - in the interval (0,1]
	
	complete.likelihood <- vector(mode="numeric",length=iter)
	BIC <- vector(mode="numeric",length=iter)
	
	# Initialise
	lambda[1] <- init.lambda
	alpha0[1] <- init.alpha
	gamma[1] <- init.gamma
	
	#start with just root node
	trees.n[[1]] <- data.frame(nu = runif(1), psi = 1, row.names="M:", label="M:", ancestor="Root:", node=1, stringsAsFactors=FALSE)
	trees.n[[1]]$phi <- trees.n[[1]]$psi[1]
	trees.n[[1]]$edge.length <- trees.n[[1]]$nu[1]
	trees.n[[1]] <- cbind(trees.n[[1]], matrix(rep(1, each=num.samples), nrow=1, byrow=TRUE))	
	node.assignments[,1] <- rep("M:",num.muts)
	
	names(trees.n[[1]]) <- c(names(trees.n[[1]])[1:7], paste("theta.S", 1:num.samples, sep=""))
	
	complete.likelihood[1] <- 0
	colnames = paste("theta.S", 1:num.samples, sep="")
	for(i in 1:num.muts){
		lfoy = log.f.of.y(y[i,], n[i,], kappa[i,], trees.n[[1]][node.assignments[i,1],colnames])
		if(!is.nan(lfoy)){
			complete.likelihood[1] <- complete.likelihood[1] + lfoy
		}
	}
# 	BIC[1] = -2 * complete.likelihood[1] + 2 * num.samples * nrow(trees.n[[1]]) * log(num.muts)	
  BIC[1] = bic(complete.likelihood[1], num.samples, nrow(trees.n[[1]]), log(num.muts))
	print(paste("init log likelihood=",complete.likelihood[1],sep=""))
	print(paste("init BIC=",BIC[1],sep=""))
		
	for (m in 2:iter) {
		print(paste("iter",m,sep=" "))
		curr.tree <- trees.n[[m-1]]

		new.assignments = node.assignments[,m-1]
		#remove node randomly, to cover more search space
    do_update = FALSE
		if(!is.na(remove.branch.frequency)){
			if(m %% remove.branch.frequency ==0){
				temp.list = remove.branch(curr.tree,new.assignments, y, n, kappa)
        do_update = TRUE
				
			}else if(!is.na(remove.node.frequency)){
				if(m %% remove.node.frequency ==0){
					temp.list = remove.node(curr.tree,new.assignments, y, n, kappa)
					do_update = TRUE
				}
			}
		}else if(!is.na(remove.node.frequency)){
			if(m %% remove.node.frequency ==0){
				temp.list = remove.node(curr.tree,new.assignments, y, n, kappa)
				do_update = TRUE
			}
		}
    
    if (do_update) {
      new.assignments = update.tree.after.removal(temp.list[[1]], temp.list[[2]], temp.list[[3]], num.muts, y, n, kappa)
    }
			
		#randomise the order of muts
		rand.inds = sample(num.muts)
		node.assignments[,m] = new.assignments
		for (i in 1:num.muts) {
			conflicts = get.conflicts(rand.inds[i], conflict.array, node.assignments[,m],curr.tree)
			temp.assignments <- sample.assignment(y[rand.inds[i],], n[rand.inds[i],], kappa[rand.inds[i],], curr.tree, node.assignments[rand.inds[i],m], lambda[m-1], alpha0[m-1], gamma[m-1], conflicts)
			curr.tree <- temp.assignments[[2]]
			node.assignments[rand.inds[i],m] <- temp.assignments[[3]]
		}
		
		#cull tree (remove empty nodes)
		temp.list <- cull.tree(curr.tree, node.assignments[,m])
		curr.tree <- temp.list$culled.tree
		mapping <- temp.list$mapping
		new.node.assignments = vector(length = num.muts, mode = "character")
		for(i in 1:nrow(mapping)){
			new.node.assignments[node.assignments[,m]==mapping$old[i]] = mapping$new[i]
		}
		node.assignments[,m] = new.node.assignments
		
		print("# Resample theta.S variables")
		for (k in grep("theta", names(curr.tree))) {
			curr.tree[, k] <- whole.tree.slice.sampler(curr.tree, curr.tree[,k], y[,k-7], n[,k-7], kappa[,k-7], node.assignments[,m], shrinkage.threshold)
		}
		
		print("# Resample hyperparameters")
		temp.hypers <- sample.hyperparameters(alpha0[m-1], lambda[m-1], gamma[m-1], allowed.ranges, curr.tree)
		alpha0[m] <- temp.hypers[1]
		lambda[m] <- temp.hypers[2]
		gamma[m] <- temp.hypers[3]
		
		print("# Resample stick lengths")
		curr.tree <- sample.sticks(curr.tree, node.assignments[,m], alpha0[m], lambda[m], gamma[m])
		
		colnames = paste("theta.S", 1:num.samples, sep="")
		complete.likelihood[m] <- 0
		for(i in 1:num.muts){
			lfoy = log.f.of.y(y[i,], n[i,], kappa[i,], curr.tree[node.assignments[i,m],colnames])
			if(!is.nan(lfoy)){
				complete.likelihood[m] <- complete.likelihood[m] + lfoy
			}
		}
		
# 		BIC[m] = -2 * complete.likelihood[m] + 2 * num.samples * nrow(curr.tree) * log(num.muts)	
    BIC[m] = bic(complete.likelihood[m], num.samples, nrow(curr.tree), log(num.muts))
		print(paste("log likelihood=",complete.likelihood[m],sep=""))
		print(paste("BIC=",BIC[m],sep=""))
		
		trees.n[[m]] <- curr.tree	
	}
	return(list(trees.n, node.assignments, alpha0, lambda, gamma, complete.likelihood, BIC))
} 

find.node <- function(u, tree, lambda, depth, alpha0, gamma, conflicts = array(1,nrow(tree))) {
	# Subroutine to find the node of the tree partition that holds the random datum u on (0,1]
	# curr.state is the current state of the mutation
	descend <- function(u1, tree1, curr.node, depth1, lambda1, alpha01, gamma, conflicts) {
		curr.nu = tree1[curr.node,"nu"]
		#adjust root node to allow for conflicts - if there are conflicting variants they should not be placed in the root node
		#doesn't work very well!!
		#if(curr.node=="M:"){curr.nu = curr.nu / conflicts["M:"]}
		if(is.na(u1) | is.na(curr.nu)){
			print("ERROR. u1 or curr.nu is NA")
			print(paste(curr.node,u1, curr.nu, sep=","))
			print(tree1)
		}
		if (u1 < curr.nu) {return(list("Success", tree1[curr.node,"label"], tree1))}
		
		if (depth1 > 10) {warning("Max depth reached", immediate.=FALSE); return(list("Success", tree1[curr.node,"label"], tree1))}
		u1 <- (u1 - tree1[curr.node,"nu"]) / (1-tree1[curr.node, "nu"])
		
		# Define immediate children of current node
		temp.matrix <- tree1[tree1$ancestor == curr.node, ]
		while (dim(temp.matrix)[1] == 0 | u1 > (1-prod(1-temp.matrix$psi))) {
			new.psi <- rbeta(1,1,gamma)
			#don't allow nu = 1 (causes problems with resampling alpha and probably other probs)
			new.nu <- min(rbeta(1,1,alpha01 * (lambda1^depth1)),0.999)
			
			if (dim(temp.matrix)[1] == 0) {# spawn new descendant of current node
				new.phi <- new.psi 
				}
			else {	
				new.phi <- new.psi * prod(1-temp.matrix$psi)
			}
			new.node <- max(tree1$node)+1
			new.node.direction <- dim(temp.matrix)[1]+1
			new.label <- paste(curr.node, new.node.direction, ":", sep="")
			new.edge.length <- new.nu * new.phi * tree1[curr.node, "edge.length"] * (1-tree1[curr.node, "nu"]) / tree1[curr.node, "nu"]
			new.thetas <- sapply(grep("theta", names(tree1)), 
					FUN=function(i,x,y) {max <- x[y,i] - sum(x[x$ancestor == y,i]); runif(1,0,max)}, x=tree1, y=curr.node)
			
			new.df <- data.frame(nu=new.nu, psi=new.psi, label=new.label, ancestor=curr.node, node=new.node, phi=new.phi, edge.length=new.edge.length, 
					t(new.thetas), stringsAsFactors=FALSE)
			names(new.df) <- names(temp.matrix)
			temp.matrix <- rbind(temp.matrix, new.df)
			row.names(temp.matrix) <- temp.matrix$label
			tree1 <- rbind(tree1[tree1$ancestor != curr.node,], temp.matrix)
			#no conflict with new node (it contains no muts)
			conflicts = c(conflicts,1)
			names(conflicts)[length(conflicts)] = new.label
		}
		edges <- 1-cumprod(1-temp.matrix$psi)
		max.edge <- edges[length(edges)]
		edges= edges / conflicts[temp.matrix$label]
		edges = edges * max.edge / edges[length(edges)]
				
		index <- sum(u1 > edges) + 1
		edges <- c(0,edges)
		u1 <- (u1-edges[index]) / (edges[index+1] - edges[index])
		return(list("Fail", u1, tree1, temp.matrix$label[index], depth1 + 1, lambda1, alpha01, gamma, conflicts))
	}
	
	result <- list("Fail", u, tree, "M:", 0, lambda, alpha0, gamma, conflicts)
	while (result[[1]] == "Fail") {
		result <- descend(result[[2]], result[[3]], result[[4]], result[[5]], result[[6]], result[[7]], result[[8]], result[[9]])
	}
	
	return(result)
}

sample.assignment <- function(y, n, kappa, tree1, curr.assignment, lambda, alpha0, gamma, conflicts = vector(mode="character",length=0)) {
	# y is the vector of the number of reads reporting the variant in each of the samples
	# n is the vector of the total number of reads across the base in each of the samples
	# kappa is the vector of expected proportion of reads at that base
		 
	 max.u <- 1
	 min.u <- 0
	 curr.thetas <- unlist(tree1[curr.assignment, grepl("theta", names(tree1))])
	 	
	 p.slice <- log(runif(n=1)) + log.f.of.y(y,n,kappa,curr.thetas)
	 	 
	 output <- list("Fail", min.u, max.u)
	 old.tree <- tree1
	 
	 while (output[[1]] == "Fail") {
	 	curr.u <- runif(1,min.u,max.u)
	
	 	new.assignment <- find.node(curr.u, tree1, lambda, 0, alpha0, gamma, conflicts)
	 	
	 	tree1 <- new.assignment[[3]]
	 	new.node <- new.assignment[[2]]
	 	thetas.new.asst <- unlist(tree1[new.node, grepl("theta", names(tree1))])
	 	curr.p <- log.f.of.y(y,n,kappa,thetas.new.asst)
	 	if (curr.p == "NaN") {return(list("Success", old.tree, curr.assignment))}
	 	if (curr.p > p.slice) {output <- list("Success", tree1, new.node)}
	 	else {if (max.u - min.u < 10E-6) {
	 			#print("Slice sampler shrank down: keep current state")
	 			output <- list("Success", old.tree, curr.assignment)
	 			}
	 		else {
	 			if (new.node < curr.assignment) {
	 				min.u <- curr.u
	 				output <- list("Fail")
	 			}	
	 			else {
	 				max.u <- curr.u
	 				output <- list("Fail")
	 			}
	 		}
	 	}
	 	#if FAIL, don't add node
	 	tree1<-old.tree
	 }
	 
	 return(output)
	 
}

bounded.slice.sampler <- function(log.fn, A, B, curr.x, ...) {
	# log.fn(x=curr.x, ...) is the log of the function that is proportional to the pdf of the variable to be sampled
	# A is the lower bound of the interval in which to slice sample; B the upper bound
	# ... are the other arguments to be passed to log.fn()
	# curr.x is the current value of the variable
	
	p.slice <- log(runif(n=1)) + log.fn(x=curr.x, ...)
	min.u <- A
	max.u <- B
	
	output <- list("Fail")
	
	while (output[[1]] == "Fail") {
		possible.x <- runif(n=1, min.u, max.u)
		if (max.u - min.u < 10E-6) {
	 		#print("Slice sampler shrank down: keep current state")
	 		output <- list("Success", curr.x)
		}	
	 	else {
	 		if (p.slice < log.fn(x = possible.x, ...)) {
	 			output <- list("Success", possible.x)
	 		}
	 		else {
	 			if (possible.x < curr.x) {min.u <- possible.x}
	 			else {max.u <- possible.x}
	 		}
	 	}		
	 }
	return(output) 
}

whole.tree.slice.sampler <- function(curr.tree, curr.thetas, y, n, kappa, curr.assignments, shrinkage.threshold) {
	# y is a vector of y values for data points; n and kappa similar
	# curr.assignments is a vector of assignments for each data point
	# curr.thetas is a vector with the current theta values
	
	names(curr.thetas) <- row.names(curr.tree)
	x.by.assignment <- curr.thetas[curr.assignments]
	num.nodes <- length(curr.thetas)
	ancs <- unique(curr.tree$ancestor[curr.tree$ancestor != "Root:"])
	
	gradient.log.f.of.xi <- function(y1, n1, kappa1, xi) {
		sum(y1 / xi - kappa1*(n1-y1) / (1-kappa1*xi))
	}	
	
	p.slice <- log(runif(n=1)) + log.f.of.y(y, n, kappa, x.by.assignment)
	
	# Initialise constraints matrix
	constraints <- data.frame(matrix(0, nrow=num.nodes, ncol=length(ancs)), row.names=row.names(curr.tree))
	names(constraints) <- ancs
	for (i in ancs) {
		constraints[curr.tree$label[curr.tree$ancestor == i], i] <- -1
		constraints[i,i] <- 1
	}
	constraints.matrix <- rbind(t(as.matrix(constraints)), diag(1,num.nodes,num.nodes), diag(-1,num.nodes,num.nodes))
	min.u <- rep(0,num.nodes)
	max.u <- rep(-1,num.nodes)
	names(min.u) <- names(max.u) <- row.names(curr.tree)
	
	exact.equality.for.M <- rep(0,num.nodes)
	exact.equality.for.M[row.names(curr.tree) == "M:"] <- 1
	
	# Run the slice sampler
	# max.repeats introduced, to prevent infinite looping
	#max.repeats = 10000
	max.repeats = 1000
	repeat.count=0
	#debug!!!
	#save(exact.equality.for.M,constraints.matrix,ancs,curr.thetas,file=paste("beforeRepeat_",format(Sys.time(), "%X"),"RData",sep=""))
	#print("whole.tree.slice.sampler before repeat")
	repeat {
		possible.theta <- xsample(E=exact.equality.for.M, F=1, G=constraints.matrix, H=c(rep(0, length(ancs)), min.u, max.u), iter=500, x0=curr.thetas)$X[500,]
		#xsample fails if x0 is not a valid solution, but re-calling xsample doesn't help!
		#possible.theta = try(xsample(E=exact.equality.for.M, F=1, G=constraints.matrix, H=c(rep(0, length(ancs)), min.u, max.u), iter=500, x0=curr.thetas)$X[500,],T)
		#count=0
		#while(class(possible.theta)=="try-error"){
		#	possible.theta = try(xsample(E=exact.equality.for.M, F=1, G=constraints.matrix, H=c(rep(0, length(ancs)), min.u, max.u), iter=500, x0=curr.thetas)$X[500,],T)
		#	count = count+1
		#	if(count %% 100 == 0){
		#		print(paste("xsample iters=",count,sep=""))
		#	}
		#}
		
		if (max(-max.u - min.u) < 1E-6) {
	 		print("Slice sampler shrank down: keep current state")
	 		return(curr.thetas)
		}else if(repeat.count>=max.repeats){
	 		print("max repeats reached: keep current state")
	 		#print(paste("max.u=",max.u,", min.u=",min.u,sep=""))
	 		return(curr.thetas)
	 	}else {
	 		possible.x.by.assignment <- possible.theta[curr.assignments]
	 		if (p.slice < log.f.of.y(y,n,kappa,possible.x.by.assignment)) {
	 			return(possible.theta)
	 		}
	 		else {
	 			# Shrink the constraints on the thetas by calculating the gradient of log(f(theta_i))
	 			curr.gradients <- sapply(1:num.nodes, function(curr.ass,node.names,y1,n1,kappa1,thetas,i) {gradient.log.f.of.xi(y1[curr.ass == node.names[i]], 
	 				n1[curr.ass == node.names[i]], kappa1[curr.ass == node.names[i]], thetas[i])}, curr.ass=curr.assignments, node.names=row.names(curr.tree), y1=y, n1=n, 
	 				kappa1=kappa,thetas=possible.theta)
	 			curr.gradients[row.names(curr.tree) == "M:"] <- 0
	 			shrinkage.axis <- (abs(curr.gradients) * (-max.u-min.u) > shrinkage.threshold)
	 			shrinkage.axis[which.max(abs(curr.gradients) * (-max.u-min.u))] <- TRUE
	 			min.u[shrinkage.axis & possible.theta < curr.thetas] <- possible.theta[shrinkage.axis & possible.theta < curr.thetas] 
	 			max.u[shrinkage.axis & possible.theta > curr.thetas] <- -possible.theta[shrinkage.axis & possible.theta > curr.thetas]
	 		}
	 	}
	 	repeat.count=repeat.count + 1
	 	if(repeat.count %% 50 ==0){
			print(paste("whole.tree.slice.sampler: ",repeat.count," iters",sep=""))
		}
	}
}

sample.hyperparameters <- function(curr.alpha, curr.lambda, curr.gamma, allowed.ranges, curr.tree) {
	curr.tree$depth <- sapply(curr.tree$label, FUN = function(x) {sum(gregexpr(":",x)[[1]]>0)}) - 1
	log.alpha.fn <- function(x, tree, lambda) {sum(dbeta(tree$nu, 1, x*(lambda^tree$depth), log=TRUE))}
	log.lambda.fn <- function(x, tree, alpha) {sum(dbeta(tree$nu, 1, alpha * (x^tree$depth), log=TRUE))}
	log.gamma.fn <- function(x, tree) {sum(dbeta(tree$nu, 1, x, log=TRUE))}
	
	new.alpha <- bounded.slice.sampler(log.alpha.fn, tree=curr.tree, allowed.ranges["min.alpha"], allowed.ranges["max.alpha"], curr.x=curr.alpha, lambda=curr.lambda)
	new.lambda <- bounded.slice.sampler(log.lambda.fn, tree=curr.tree, allowed.ranges["min.lambda"], allowed.ranges["max.lambda"], curr.x=curr.lambda, alpha=new.alpha[[2]])
	new.gamma <- bounded.slice.sampler(log.gamma.fn, tree=curr.tree, allowed.ranges["min.gamma"], allowed.ranges["max.gamma"], curr.x=curr.gamma)
	return(c(alpha=new.alpha[[2]], lambda=new.lambda[[2]], gamma=new.gamma[[2]])) 
}

sample.sticks <- function(curr.tree, curr.assignments, alpha, lambda, gamma) {
	curr.tree$depth <- sapply(curr.tree$label, FUN = function(x) {sum(gregexpr(":",x)[[1]]>0)}) - 1
	curr.tree$num.assignments <- sapply(curr.tree$label, function(a,x) {sum(a == x)}, a=curr.assignments)
	curr.tree <- curr.tree[order(curr.tree$depth, curr.tree$label),]
  
	for (i in 1:dim(curr.tree)[1]) {
		Ne.descendants <- sum(curr.tree$num.assignments[grepl(curr.tree$label[i], curr.tree$label)])
		Ne <- curr.tree$num.assignments[i]
		curr.tree$nu[i] <- rbeta(n=1, Ne+1, Ne.descendants + alpha*(lambda^curr.tree$depth[i]))
		if (curr.tree$label[i] != "M:") {
			# Need to count the number of data points assigned to nodes that are descendants of younger siblings of current node
			Ne.younger.sibs.descendants <- sum(curr.tree$num.assignments[grepl(curr.tree$ancestor[i], curr.tree$label) & younger.descendants(curr.tree$label[i], curr.tree$label) & !grepl(curr.tree$label[i], curr.tree$label)])
			curr.tree$psi[i] <- rbeta(n=1, Ne.descendants+1, gamma + Ne.younger.sibs.descendants)
			if (substring(curr.tree$label[i], first=nchar(curr.tree$label[i])-2) == ":1:") {# Oldest sibling
				curr.tree$phi[i] <- curr.tree$psi[i]
			}
			else {
				curr.tree$phi[i] <- curr.tree$psi[i] * prod(1-curr.tree$psi[older.siblings(curr.tree$label[i], curr.tree$label)])
			}
			curr.tree$edge.length[i] <- curr.tree$nu[i] * curr.tree$phi[i] * curr.tree$edge.length[curr.tree$label == curr.tree$ancestor[i]] * (1-curr.tree$nu[curr.tree$label == curr.tree$ancestor[i]]) /  curr.tree$nu[curr.tree$label == curr.tree$ancestor[i]]

		} 
		
	}
	curr.tree$edge.length[curr.tree$label == "M:"] <- curr.tree$nu[curr.tree$label == "M:"]
	return(curr.tree[,!is.element(names(curr.tree), c("depth", "num.assignments"))])
}


# Not used - sd11
# sample.stick.orders <- function(curr.tree, lambda, alpha, gamma, curr.assign) {
# 	curr.tree$depth <- sapply(curr.tree$label, FUN = function(x) {sum(gregexpr(":",x)[[1]]>0)}) - 1
# 	curr.tree <- curr.tree[order(curr.tree$depth, curr.tree$label),]
# 
# 	descend <- function(tree, depth, node, lambda1, alpha1, gamma1, curr.assignments) {
# 		children <- tree$label[tree$ancestor == node]
# 		if (length(children) == 0) {return(list(tree, curr.assignments))}
# 		curr.children <- children
# 		all.weights <- tree[children, "phi"]
# 		new.order <- c()
# 
# 		while (TRUE) {
# 			if (length(curr.children) == 0) {break}
# 			else {
# 				u <- runif(n=1)
# 				while (TRUE) {
# 					sub.indices <- (1:length(children))[!is.element(1:length(children), new.order)] # Indices not already assigned to new.order
# 					sub.weights <- c(all.weights[sub.indices], 1-sum(all.weights))
# 					sub.weights <- sub.weights / sum(sub.weights)
# 					index <- sum(u > cumsum(sub.weights)) + 1
# 					if (index == length(sub.indices)+1) {
# 						new.psi <- rbeta(1,1,lambda1)
# 						new.nu <- rbeta(1,1,alpha1 * (lambda1^depth))
# 						new.phi <- new.psi * prod(1-tree$psi[tree$ancestor == node])
# 						new.node <- max(tree$node)+1
# 						new.node.direction <- length(children)+1
# 						new.label <- paste(node, new.node.direction, ":", sep="")
# 						new.edge.length <- new.nu * new.phi * tree[node, "edge.length"] * (1-tree[node, "nu"]) / tree[node, "nu"]
# 						new.thetas <- t(sapply(grep("theta", names(tree)), FUN=function(i,x,y) {max <- x[y,i] - sum(x[x$ancestor == y,i]); runif(1,0,max)}, x=tree, y=node))
# 						new.df <- data.frame(nu=new.nu, psi=new.psi, label=new.label, ancestor=node, node=new.node, phi=new.phi, edge.length=new.edge.length, 
# 							new.thetas, depth = depth+1, stringsAsFactors=FALSE)
# 						names(new.df) <- names(tree)
# 						if(new.label %in% tree$label){
# 							print("ERROR. new.label=")
# 							print(new.label)
# 							print(tree$label)
# 						}
# 						tree <- rbind(tree, new.df)
# 						row.names(tree) <- tree$label
# 						tree <- tree[order(tree$depth, tree$label),]
# 						all.weights <- tree$phi[tree$ancestor == node]
# 						children <- c(children, new.label)
# 					}
# 					else {
# 						index <- sub.indices[index]
# 						break
# 					}
# 				}
# 			}
# 			new.order <- c(new.order, index)
# 			curr.children <- curr.children[curr.children != children[index]]
# 		}
#     
# 		new.order <- c(new.order, (1:max(new.order))[!is.element(1:max(new.order), new.order)])
# 		reordered.children <- gsub("M","T",children)
# 		# Reorder children
# 		for (k in 1:length(new.order)) {
# 			now.child <- children[new.order[k]]
# 			now.child.split <- strsplit(reordered.children[new.order[k]], ":")[[1]]
# 			now.child.split[length(now.child.split)] <- k
# 			replacement.child <- paste(paste(now.child.split, collapse=":"),":",sep="")
# 			tree$label <- gsub(now.child, replacement.child, tree$label)
# 			tree$ancestor <- gsub(now.child, replacement.child, tree$ancestor)
# 			curr.assignments <- gsub(now.child, replacement.child, curr.assignments)
# 		}
# 		tree$label <- gsub("T", "M", tree$label)
# 		tree$ancestor <- gsub("T", "M", tree$ancestor)
# 		curr.assignments <- gsub("T", "M", curr.assignments)
# 		row.names(tree) <- tree$label 
# 		tree <- tree[order(tree$depth, tree$label),]
# 		
# 		for (k in tree$label[tree$ancestor == node]) {
# 			temp <- descend(tree, depth+1, k, lambda1, alpha1, gamma1, curr.assignments)
# 			tree <- temp[[1]]
# 			curr.assignments <- temp[[2]]
# 		}
# 		
# 		for (k in tree$label[tree$ancestor == node]) {
# 			if (substring(tree[k, "label"], first=nchar(tree[k, "label"])-2) == ":1:") {# Oldest sibling
# 				tree[k, "phi"] <- tree[k, "psi"]
# 				tree[k, "edge.length"] <- tree[k, "nu"] * tree[k, "phi"] * tree[node, "edge.length"] * (1-tree[node, "nu"]) / tree[node, "nu"]
# 			}
# 			else {
# 				tree[k, "phi"] <- tree[k, "psi"] * prod(1-tree$psi[older.siblings(k, tree$label)])
# 				tree[k, "edge.length"] <- tree[k, "nu"] * tree[k, "phi"] * tree[node, "edge.length"] * (1-tree[node, "nu"]) / tree[node, "nu"]
# 			}
# 		}
# 		return(list(tree, curr.assignments))
# 		
# 	}	
# 	a <- descend(curr.tree, 0, "M:", lambda, alpha, gamma,curr.assign)
# 	return(list(a[[1]][,!grepl("depth",names(a[[1]]))], a[[2]]))
# }

## Compile some internal functions
library(compiler)
bic = cmpfun(bic)
log.f.of.y = cmpfun(log.f.of.y)