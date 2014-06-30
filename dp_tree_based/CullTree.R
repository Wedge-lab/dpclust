#cull.tree removes empty nodes (those with no assigned mutations)
cull.tree <- function(tree1, curr.assignments) {
	saved.tree = tree1
	# tree1 is the current tree
	# curr.assignment is the vector of current assignments of the tree
	nodes.in.tree = createInventory(tree1, curr.assignments)
  
	#use node numbers - they don't change
	original.keep.nodes = tree1$node[tree1$label %in%nodes.in.tree$label[nodes.in.tree$keep.node]]

  # Replace with descendants
  print("Before descendants")
# TODO this doesn't work
# 	sapply(which(nodes.in.tree$replace.with.descendants), FUN=replaceWithDescendants, nodes.in.tree, tree1)
  
  #040112 - not sure that psi and phi values are set correctly
	for(i in which(nodes.in.tree$replace.with.descendants)){
		#name in nodes.in.tree
		remove.node <- nodes.in.tree$label[i]
		#name in tree1
		remove.node = tree1[remove.node,"label"]
		
		spl = unlist(strsplit(remove.node,":"))
		remove.depth = length(spl)
		remove.pos = as.integer(spl[length(spl)])
		remove.index = which(tree1$label==remove.node)
		ancestor.index = which(tree1$label==tree1$ancestor[remove.index])

		#just find nodes in culled tree
		younger.indices2 = which(younger.direct.descendants(remove.node, tree1$label) & tree1$label!=remove.node)		
		younger.indirect.indices2 = which(younger.descendants(remove.node, tree1$label) & tree1$label!=remove.node)		
		younger.indirect.indices2 = younger.indirect.indices2[!(younger.indirect.indices2 %in% younger.indices2)]		
		younger.descendants = tree1$label[younger.indices2]
		younger.indices = match(row.names(tree1)[younger.indices2], nodes.in.tree$label)
		no.to.shift = 0
		for(j in 1:length(younger.indices2)){	
			ind = younger.indices[j]
			ind2 = younger.indices2[j]
			tree1$nu[ind2] = tree1$nu[ind2] * tree1$nu[remove.index]
			tree1$psi[ind2] = tree1$psi[ind2] * tree1$psi[remove.index]
			tree1$phi[ind2] = tree1$psi[ind2] * tree1$phi[remove.index]
			tree1$edge.length[ind2] = tree1$nu[ind2] * tree1$psi[ind2] * tree1$edge.length[ancestor.index]*(1-tree1$nu[ancestor.index])/tree1$nu[ancestor.index]			
			
			spl2 = unlist(strsplit(tree1$label[ind2],":"))
			no.to.shift = max(no.to.shift,as.integer(spl2[remove.depth+1])-1)
			tree1$label[ind2] = paste(c(spl[1:(remove.depth-1)],remove.pos+as.integer(spl2[remove.depth+1])-1),collapse=":")
			if(length(spl2)>(remove.depth+1)){
				tree1$label[ind2] = paste(c(tree1$label[ind2],spl2[(remove.depth+2):length(spl2)]),collapse=":")
			}
			tree1$label[ind2] = paste(tree1$label[ind2],":",sep="")
			spl3 = unlist(strsplit(tree1$label[ind2],":"))
			tree1$ancestor[ind2] = paste(paste(spl3[1:(length(spl3)-1)],collapse=":"),":",sep="")
			#pass on replace.with.indices property to descendants
			if(nodes.in.tree$replace.with.siblings[i] & length(spl2) == (length(spl)+1) & !nodes.in.tree$keep.node[ind] & !nodes.in.tree$replace.with.descendants[ind] & !nodes.in.tree$replace.with.siblings[ind]){
				nodes.in.tree$replace.with.siblings[ind] = T
			}
		}
		if(no.to.shift>0 & length(younger.indirect.indices2)>0){	
			younger.indirect.descendants = tree1$label[younger.indirect.indices2]
			younger.indirect.indices = match(younger.indirect.descendants,nodes.in.tree$label)
			for(j in 1:length(younger.indirect.indices)){
				ind = younger.indirect.indices[j]
				ind2 = younger.indirect.indices2[j]
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
		#remove nodes one at a time
		tree1 <- tree1[-remove.index,]		
	}
    # Replace with siblings
    print("Before siblings")
  # TODO this doesn't work
# 	  sapply(which(nodes.in.tree$replace.with.siblings & !nodes.in.tree$replace.with.descendants), FUN=replaceWithSiblings, nodes.in.tree, tree1)
	for(i in which(nodes.in.tree$replace.with.siblings & !nodes.in.tree$replace.with.descendants)){		
		#name in nodes.in.tree
		remove.node <- nodes.in.tree$label[i]
		#name in tree1
		remove.node = tree1[remove.node,"label"]
		
		spl = unlist(strsplit(remove.node,":"))
		remove.depth = length(spl)
		remove.pos = as.integer(spl[length(spl)])
		remove.index = which(tree1$label==remove.node)
		ancestor.index = which(tree1$label==tree1$ancestor[remove.index])
		## Same till here
		younger.indices2 = which(younger.direct.descendants(remove.node, tree1$label) & tree1$label!=remove.node)	
		younger.indirect.indices2 = which(younger.descendants(remove.node, tree1$label) & tree1$label!=remove.node)		
		younger.indirect.indices2 = younger.indirect.indices2[!(younger.indirect.indices2 %in% younger.indices2)]
		younger.indirect.descendants = tree1$label[younger.indirect.indices2]
		younger.indirect.indices = match(younger.indirect.descendants,nodes.in.tree$label)
    ## Only difference is usage of indirect in variable names - check if used below
		younger.sibling.indices2 = which(younger.siblings(remove.node, tree1$label) & tree1$label!=remove.node)
		younger.sibling.descendants = tree1$label[younger.sibling.indices2]
		younger.sibling.indices = match(younger.sibling.descendants,nodes.in.tree$label)
    ## These are extra
		
		#use indices from culled tree (otherwise some renamed nodes are missed)
		if(length(younger.indirect.indices2)>0){
			for(j in 1:length(younger.indirect.indices2)){
				ind2 = younger.indirect.indices2[j]
				spl2 = unlist(strsplit(tree1$label[ind2],":"))
				tree1$label[ind2] = paste(c(spl[1:(remove.depth-1)],as.integer(spl2[remove.depth])-1),collapse=":")
				if(length(spl2)>remove.depth){
					tree1$label[ind2] = paste(c(tree1$label[ind2],spl2[(remove.depth+1):length(spl2)]),collapse=":")
				}
				tree1$label[ind2] = paste(tree1$label[ind2],":",sep="")
				spl3 = unlist(strsplit(tree1$label[ind2],":"))
				tree1$ancestor[ind2] = paste(paste(spl3[1:(length(spl3)-1)],collapse=":"),":",sep="")

				#adjust psi for younger siblings
				if(ind2 %in% younger.sibling.indices2){
					tree1$psi[ind2] = tree1$psi[ind2] * (1-tree1$psi[remove.index])
				}
			}
		}
		#remove nodes one at a time
		tree1 <- tree1[-remove.index,]
	}
  print("Before wrapup")
	tree1<-tree1[tree1$node %in% original.keep.nodes,]
	row.names(tree1) = tree1$label
	#get mapping between old and new labels
	old.labels = saved.tree$label[match(tree1$node,saved.tree$node)]
	df = data.frame(old=old.labels,new=tree1$label,stringsAsFactors=F)
  print("Done")
	return(list(culled.tree=tree1,mapping=df))
}


createInventory <- function(tree1, curr.assignments) {
  # Crate an inventory of which nodes are to be kept and which ones should be replaced by a sibling
  nodes.with.data <- unique(curr.assignments)
  nodes.in.tree <- data.frame(label=sort(tree1$label), keep.node=rep(FALSE,dim(tree1)[1]), replace.with.siblings=rep(FALSE,dim(tree1)[1]), replace.with.descendants=rep(FALSE,dim(tree1)[1]), stringsAsFactors=FALSE)
  # Keep nodes that have data assigned to them and the root
  nodes.in.tree$keep.node[is.element(nodes.in.tree$label, nodes.with.data)] <- TRUE
  nodes.in.tree$keep.node[1] <- TRUE
  # Now create inventory of which nodes are to be replaced
  res = sapply(2:dim(nodes.in.tree)[1], FUN=createInventoryInner, nodes.in.tree, nodes.with.data)
  res = rbind(c(F,F), as.data.frame(t(res)))
  nodes.in.tree$replace.with.siblings = unlist(res$replace.with.siblings)
  nodes.in.tree$replace.with.descendants = unlist(res$replace.with.descendants)
  return(nodes.in.tree)
}

createInventoryInner <- function(i, nodes.in.tree, nodes.with.data) {
  with.descendants = F
  with.siblings = F
  temp.node <- nodes.in.tree$label[i]
  younger.relatives <- nodes.in.tree$label[younger.descendants(temp.node, nodes.in.tree$label)]
  younger.descendants <- nodes.in.tree$label[younger.direct.descendants(temp.node, nodes.in.tree$label)]
  if (sum(is.element(younger.relatives, nodes.with.data)) > 0){
    if(!(nodes.in.tree$keep.node[i])){  			
      if (sum(is.element(younger.descendants, nodes.with.data)) > 0){
        with.descendants = T
      }
      if(sum(is.element(younger.descendants, nodes.with.data)) < sum(is.element(younger.relatives, nodes.with.data))){
        with.siblings = T
      }
    }
  }
  return(list(replace.with.descendants=with.descendants, replace.with.siblings=with.siblings))
}

replaceWithDescendants <- function(i, nodes.in.tree, tree1) {
  print("In descendants")
  #040112 - not sure that psi and phi values are set correctly
  #name in nodes.in.tree
  remove.node <- nodes.in.tree$label[i]
  #name in tree1
  remove.node = tree1[remove.node,"label"]

  spl = unlist(strsplit(remove.node,":"))
  remove.depth = length(spl)
  remove.pos = as.integer(spl[length(spl)])
  remove.index = which(tree1$label==remove.node)
  ancestor.index = which(tree1$label==tree1$ancestor[remove.index])

  #just find nodes in culled tree
  younger.indices2 = which(younger.direct.descendants(remove.node, tree1$label) & tree1$label!=remove.node)    
  younger.indirect.indices2 = which(younger.descendants(remove.node, tree1$label) & tree1$label!=remove.node)		
  younger.indirect.indices2 = younger.indirect.indices2[!(younger.indirect.indices2 %in% younger.indices2)]		
  younger.descendants = tree1$label[younger.indices2]
  younger.indices = match(row.names(tree1)[younger.indices2], nodes.in.tree$label)

  no.to.shift = 0
  for(j in 1:length(younger.indices2)){	
    ind = younger.indices[j]
    ind2 = younger.indices2[j]
    tree1$nu[ind2] = tree1$nu[ind2] * tree1$nu[remove.index]
    tree1$psi[ind2] = tree1$psi[ind2] * tree1$psi[remove.index]
    tree1$phi[ind2] = tree1$psi[ind2] * tree1$phi[remove.index]
    tree1$edge.length[ind2] = tree1$nu[ind2] * tree1$psi[ind2] * tree1$edge.length[ancestor.index]*(1-tree1$nu[ancestor.index])/tree1$nu[ancestor.index]			
    
    spl2 = unlist(strsplit(tree1$label[ind2],":"))
    no.to.shift = max(no.to.shift,as.integer(spl2[remove.depth+1])-1)
    tree1$label[ind2] = paste(c(spl[1:(remove.depth-1)],remove.pos+as.integer(spl2[remove.depth+1])-1),collapse=":")
    if(length(spl2)>(remove.depth+1)){
      tree1$label[ind2] = paste(c(tree1$label[ind2],spl2[(remove.depth+2):length(spl2)]),collapse=":")
    }
    tree1$label[ind2] = paste(tree1$label[ind2],":",sep="")
    spl3 = unlist(strsplit(tree1$label[ind2],":"))
    tree1$ancestor[ind2] = paste(paste(spl3[1:(length(spl3)-1)],collapse=":"),":",sep="")
    #pass on replace.with.indices property to descendants
    if(nodes.in.tree$replace.with.siblings[i] & length(spl2) == (length(spl)+1) & !nodes.in.tree$keep.node[ind] & !nodes.in.tree$replace.with.descendants[ind] & !nodes.in.tree$replace.with.siblings[ind]){
      nodes.in.tree$replace.with.siblings[ind] = T
    }
  }

  if(no.to.shift>0 & length(younger.indirect.indices2)>0){	
    younger.indirect.descendants = tree1$label[younger.indirect.indices2]
    younger.indirect.indices = match(younger.indirect.descendants,nodes.in.tree$label)
    for(j in 1:length(younger.indirect.indices)){
      ind = younger.indirect.indices[j]
      ind2 = younger.indirect.indices2[j]
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
  #remove nodes one at a time
  tree1 <- tree1[-remove.index,]	
}

replaceWithSiblings <- function(i, nodes.in.tree, tree1) {
  #name in nodes.in.tree
  remove.node <- nodes.in.tree$label[i]
  #name in tree1
  remove.node = tree1[remove.node,"label"]
  
  spl = unlist(strsplit(remove.node,":"))
  remove.depth = length(spl)
  remove.pos = as.integer(spl[length(spl)])
  remove.index = which(tree1$label==remove.node)
  ancestor.index = which(tree1$label==tree1$ancestor[remove.index])

  younger.indices2 = which(younger.direct.descendants(remove.node, tree1$label) & tree1$label!=remove.node)	
  younger.indirect.indices2 = which(younger.descendants(remove.node, tree1$label) & tree1$label!=remove.node)		
  younger.indirect.indices2 = younger.indirect.indices2[!(younger.indirect.indices2 %in% younger.indices2)]
  younger.indirect.descendants = tree1$label[younger.indirect.indices2]
  younger.indirect.indices = match(younger.indirect.descendants,nodes.in.tree$label)

  younger.sibling.indices2 = which(younger.siblings(remove.node, tree1$label) & tree1$label!=remove.node)
  younger.sibling.descendants = tree1$label[younger.sibling.indices2]
  younger.sibling.indices = match(younger.sibling.descendants,nodes.in.tree$label)
  
  #use indices from culled tree (otherwise some renamed nodes are missed)
  if(length(younger.indirect.indices2)>0){
    for(j in 1:length(younger.indirect.indices2)){
      ind2 = younger.indirect.indices2[j]
      spl2 = unlist(strsplit(tree1$label[ind2],":"))
      tree1$label[ind2] = paste(c(spl[1:(remove.depth-1)],as.integer(spl2[remove.depth])-1),collapse=":")
      if(length(spl2)>remove.depth){
        tree1$label[ind2] = paste(c(tree1$label[ind2],spl2[(remove.depth+1):length(spl2)]),collapse=":")
      }
      tree1$label[ind2] = paste(tree1$label[ind2],":",sep="")
      spl3 = unlist(strsplit(tree1$label[ind2],":"))
      tree1$ancestor[ind2] = paste(paste(spl3[1:(length(spl3)-1)],collapse=":"),":",sep="")
      
      #adjust psi for younger siblings
      if(ind2 %in% younger.sibling.indices2){
        tree1$psi[ind2] = tree1$psi[ind2] * (1-tree1$psi[remove.index])
      }
    }
  }
  #remove nodes one at a time
  tree1 <- tree1[-remove.index,]
}
