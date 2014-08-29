annotateTree<-function(tree,node.assignments,mutation.annotation){
	tree$annotation = NA
	if(length(grep("mutation.annotation",ls()))>0){
		for(ind in which(mutation.annotation!="")){
			node = node.assignments[ind]
			if(is.na(tree[node,"annotation"])){
				tree[node,"annotation"] = mutation.annotation[ind]
			}else{
				tree[node,"annotation"] = paste(tree[node,"annotation"],mutation.annotation[ind],sep=",")
			}
		}
	}
	return(tree)
}
