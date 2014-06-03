remove.node <- function(tree1, curr.assignments, y, n) {
	saved.tree = tree1
	#remove a node that has a low likelihood
	log.probs = vector(mode="numeric",length=nrow(tree1))
	for (r in 1:nrow(tree1)) {
		thetas = unlist(tree1[r, grepl("theta", names(tree1))])
		inds = match(tree1$label[r],curr.assignments)
		if(length(inds)>0){
			for(i in inds){			
				log.probs[r] = log.probs[r] + log.f.of.y(y[i,],n[i,],kappas[i,],thetas)
			}
			log.probs[r] = log.probs[r] / length(inds)
		}else{
			log.probs[r] = NA
		}
	}
	
	log.probs[is.na(log.probs)] = min(log.probs)
	log.probs = min(log.probs) - log.probs
	probs = exp(log.probs)
	#don't remove root
	probs[1]=0	
	cum.probs = cumsum(probs)
	remove.index = sum(runif(1,0,max(cum.probs))>cum.probs)+1	
	
	remove.node <- tree1$label[remove.index]
	
	node.assignments = which(curr.assignments == remove.node)
		
	spl = unlist(strsplit(remove.node,":"))
	remove.depth = length(spl)
	remove.pos = as.integer(spl[length(spl)])
	ancestor.index = which(tree1$label==tree1$ancestor[remove.index])

	younger.indices = which(younger.direct.descendants(remove.node, tree1$label) & tree1$label!=remove.node)		
	younger.indirect.indices = which(younger.descendants(remove.node, tree1$label) & tree1$label!=remove.node)		
	younger.indirect.indices = younger.indirect.indices[!(younger.indirect.indices %in% younger.indices)]		
	younger.descendants = tree1$label[younger.indices]

	no.to.shift = 0
	if(length(younger.indices)>0)
	{
		for(i in 1:length(younger.indices)){
			ind = younger.indices[i]
			tree1$nu[ind] = tree1$nu[ind] * tree1$nu[remove.index]
			tree1$psi[ind] = tree1$psi[ind] * tree1$psi[remove.index]
			tree1$phi[ind] = tree1$psi[ind] * tree1$phi[remove.index]
			tree1$edge.length[ind] = tree1$nu[ind] * tree1$psi[ind] * tree1$edge.length[ancestor.index]*(1-tree1$nu[ancestor.index])/tree1$nu[ancestor.index]			

			spl2 = unlist(strsplit(tree1$label[ind],":"))
			no.to.shift = max(no.to.shift,as.integer(spl2[remove.depth+1])-1)
			tree1$label[ind] = paste(c(spl[1:(remove.depth-1)],remove.pos+as.integer(spl2[remove.depth+1])-1),collapse=":")
			if(length(spl2)>(remove.depth+1)){
				tree1$label[ind] = paste(c(tree1$label[ind],spl2[(remove.depth+2):length(spl2)]),collapse=":")
			}
			tree1$label[ind] = paste(tree1$label[ind],":",sep="")
			spl3 = unlist(strsplit(tree1$label[ind],":"))
			tree1$ancestor[ind] = paste(paste(spl3[1:(length(spl3)-1)],collapse=":"),":",sep="")
		}
		if(no.to.shift>0 & length(younger.indirect.indices)>0){	
			younger.indirect.descendants = tree1$label[younger.indirect.indices]
			for(i in 1:length(younger.indirect.indices)){
				ind = younger.indirect.indices[i]
				ind2 = younger.indirect.indices[i]
				spl2 = unlist(strsplit(tree1$label[ind2],":"))
				tree1$label[ind2] = paste(c(spl[1:(remove.depth-1)],as.integer(spl2[remove.depth])+no.to.shift),collapse=":")
				if(length(spl2)>(remove.depth)){
					tree1$label[ind2] = paste(c(tree1$label[ind2],spl2[(remove.depth+1):length(spl2)]),collapse=":")
				}
				tree1$label[ind2] = paste(tree1$label[ind2],":",sep="")
				spl3 = unlist(strsplit(tree1$label[ind2],":"))
				tree1$ancestor[ind2] = paste(paste(spl3[1:(length(spl3)-1)],collapse=":"),":",sep="")
			}
		}
	}else if(length(younger.indirect.indices)>0){
		younger.sibling.indices = which(younger.siblings(remove.node, tree1$label) & tree1$label!=remove.node)
		younger.sibling.descendants = tree1$label[younger.sibling.indices]

		for(i in 1:length(younger.indirect.indices)){
			ind = younger.indirect.indices[i]
			spl2 = unlist(strsplit(tree1$label[ind],":"))
			tree1$label[ind] = paste(c(spl[1:(remove.depth-1)],as.integer(spl2[remove.depth])-1),collapse=":")
			if(length(spl2)>remove.depth){
				tree1$label[ind] = paste(c(tree1$label[ind],spl2[(remove.depth+1):length(spl2)]),collapse=":")
			}
			tree1$label[ind] = paste(tree1$label[ind],":",sep="")
			spl3 = unlist(strsplit(tree1$label[ind],":"))
			tree1$ancestor[ind] = paste(paste(spl3[1:(length(spl3)-1)],collapse=":"),":",sep="")

			#adjust psi for younger siblings
			if(ind %in% younger.sibling.indices){
				tree1$psi[ind] = tree1$psi[ind] * (1-tree1$psi[remove.index])
			}
		}
	}
	tree1 <- tree1[-remove.index,]
	old.labels = row.names(tree1)
	row.names(tree1) = tree1$label
	
	df = data.frame(old=old.labels,new=tree1$label,stringsAsFactors=F)	
	return(list(tree1,node.assignments,df))
}

remove.specific.node <- function(tree1, remove.node) {
	saved.tree = tree1
	
	remove.index = which(tree1$label==remove.node)
	
	spl = unlist(strsplit(remove.node,":"))
	remove.depth = length(spl)
	remove.pos = as.integer(spl[length(spl)])
	ancestor.index = which(tree1$label==tree1$ancestor[remove.index])

	younger.indices = which(younger.direct.descendants(remove.node, tree1$label) & tree1$label!=remove.node)		
	younger.indirect.indices = which(younger.descendants(remove.node, tree1$label) & tree1$label!=remove.node)		
	younger.indirect.indices = younger.indirect.indices[!(younger.indirect.indices %in% younger.indices)]		
	younger.descendants = tree1$label[younger.indices]

	no.to.shift = 0
	if(length(younger.indices)>0)
	{
		for(i in 1:length(younger.indices)){
			ind = younger.indices[i]
			tree1$nu[ind] = tree1$nu[ind] * tree1$nu[remove.index]
			tree1$psi[ind] = tree1$psi[ind] * tree1$psi[remove.index]
			tree1$phi[ind] = tree1$psi[ind] * tree1$phi[remove.index]
			tree1$edge.length[ind] = tree1$nu[ind] * tree1$psi[ind] * tree1$edge.length[ancestor.index]*(1-tree1$nu[ancestor.index])/tree1$nu[ancestor.index]			

			spl2 = unlist(strsplit(tree1$label[ind],":"))
			no.to.shift = max(no.to.shift,as.integer(spl2[remove.depth+1])-1)
			tree1$label[ind] = paste(c(spl[1:(remove.depth-1)],remove.pos+as.integer(spl2[remove.depth+1])-1),collapse=":")
			if(length(spl2)>(remove.depth+1)){
				tree1$label[ind] = paste(c(tree1$label[ind],spl2[(remove.depth+2):length(spl2)]),collapse=":")
			}
			tree1$label[ind] = paste(tree1$label[ind],":",sep="")
			spl3 = unlist(strsplit(tree1$label[ind],":"))
			tree1$ancestor[ind] = paste(paste(spl3[1:(length(spl3)-1)],collapse=":"),":",sep="")
		}
		if(no.to.shift>0 & length(younger.indirect.indices)>0){	
			younger.indirect.descendants = tree1$label[younger.indirect.indices]
			for(i in 1:length(younger.indirect.indices)){
				ind = younger.indirect.indices[i]
				ind2 = younger.indirect.indices[i]
				spl2 = unlist(strsplit(tree1$label[ind2],":"))
				tree1$label[ind2] = paste(c(spl[1:(remove.depth-1)],as.integer(spl2[remove.depth])+no.to.shift),collapse=":")
				if(length(spl2)>(remove.depth)){
					tree1$label[ind2] = paste(c(tree1$label[ind2],spl2[(remove.depth+1):length(spl2)]),collapse=":")
				}
				tree1$label[ind2] = paste(tree1$label[ind2],":",sep="")
				spl3 = unlist(strsplit(tree1$label[ind2],":"))
				tree1$ancestor[ind2] = paste(paste(spl3[1:(length(spl3)-1)],collapse=":"),":",sep="")
			}
		}
	}else if(length(younger.indirect.indices)>0){
		younger.sibling.indices = which(younger.siblings(remove.node, tree1$label) & tree1$label!=remove.node)
		younger.sibling.descendants = tree1$label[younger.sibling.indices]

		for(i in 1:length(younger.indirect.indices)){
			ind = younger.indirect.indices[i]
			spl2 = unlist(strsplit(tree1$label[ind],":"))
			tree1$label[ind] = paste(c(spl[1:(remove.depth-1)],as.integer(spl2[remove.depth])-1),collapse=":")
			if(length(spl2)>remove.depth){
				tree1$label[ind] = paste(c(tree1$label[ind],spl2[(remove.depth+1):length(spl2)]),collapse=":")
			}
			tree1$label[ind] = paste(tree1$label[ind],":",sep="")
			spl3 = unlist(strsplit(tree1$label[ind],":"))
			tree1$ancestor[ind] = paste(paste(spl3[1:(length(spl3)-1)],collapse=":"),":",sep="")

			#adjust psi for younger siblings
			if(ind %in% younger.sibling.indices){
				tree1$psi[ind] = tree1$psi[ind] * (1-tree1$psi[remove.index])
			}
		}
	}
	tree1 <- tree1[-remove.index,]
	old.labels = row.names(tree1)
	row.names(tree1) = tree1$label
	
	df = data.frame(old=old.labels,new=tree1$label,stringsAsFactors=F)
	return(list(tree1,df))
}

remove.branch <- function(tree1, curr.assignments, y, n, kappas) {
	saved.tree = tree1
				
	#remove a branch that has a low likelihood
	log.probs = vector(mode="numeric",length=nrow(tree1))
	for (r in 1:nrow(tree1)) {
		inds = grep(tree1$label[r],curr.assignments)
		if(length(inds)>0){
			for(i in inds){
				thetas = unlist(tree1[match(curr.assignments[i],row.names(tree1)), grepl("theta", names(tree1))])
				log.probs[r] = log.probs[r] + log.f.of.y(y[i,],n[i,],kappas[i,],thetas)
			}
			log.probs[r] = log.probs[r] / length(inds)
		}else{
			log.probs[r] = NA
		}
	}
	log.probs[is.na(log.probs)] = min(log.probs,na.rm=T)
	log.probs = min(log.probs) - log.probs
	probs = exp(log.probs)
	#don't remove root
	probs[1]=0
	cum.probs = cumsum(probs)
	remove.index = sum(runif(1,0,max(cum.probs))>cum.probs)+1
	
	
	remove.node <- tree1$label[remove.index]
	spl = unlist(strsplit(remove.node,":"))
	remove.depth = length(spl)
	remove.pos = as.integer(spl[length(spl)])
	ancestor.index = which(tree1$label==tree1$ancestor[remove.index])

	younger.indices = which(younger.direct.descendants(remove.node, tree1$label) & tree1$label!=remove.node)
	
	remove.nodes = c(remove.node,tree1$label[younger.indices])
	node.assignments = which(curr.assignments %in% remove.nodes)
	
	younger.indirect.indices = which(younger.descendants(remove.node, tree1$label) & tree1$label!=remove.node)		
	younger.indirect.indices = younger.indirect.indices[!(younger.indirect.indices %in% younger.indices)]		
			
	if(length(younger.indirect.indices)>0){
		younger.sibling.indices = which(younger.siblings(remove.node, tree1$label) & tree1$label!=remove.node)
		younger.sibling.descendants = tree1$label[younger.sibling.indices]

		for(i in 1:length(younger.indirect.indices)){
			ind = younger.indirect.indices[i]
			spl2 = unlist(strsplit(tree1$label[ind],":"))
			tree1$label[ind] = paste(c(spl[1:(remove.depth-1)],as.integer(spl2[remove.depth])-1),collapse=":")
			if(length(spl2)>remove.depth){
				tree1$label[ind] = paste(c(tree1$label[ind],spl2[(remove.depth+1):length(spl2)]),collapse=":")
			}
			tree1$label[ind] = paste(tree1$label[ind],":",sep="")
			spl3 = unlist(strsplit(tree1$label[ind],":"))
			tree1$ancestor[ind] = paste(paste(spl3[1:(length(spl3)-1)],collapse=":"),":",sep="")

			#adjust psi for younger siblings
			if(ind %in% younger.sibling.indices){
				tree1$psi[ind] = tree1$psi[ind] * (1-tree1$psi[remove.index])
			}
		}
	}
	tree1 <- tree1[-c(remove.index,younger.indices),]
	old.labels = row.names(tree1)
	row.names(tree1) = tree1$label
	
	df = data.frame(old=old.labels,new=tree1$label,stringsAsFactors=F)	
	return(list(tree1,node.assignments,df))
}