library(igraph)
plotTree<-function(tree,main="",font.size=0.75,plotAsPercentage=T){
	theta.cols = grep("theta",names(tree))
	labels <-rownames(tree)
	depth <- c(0,t(sapply(1:nrow(tree),function(n,i){length(strsplit(n[i],":")[[1]])},n=rownames(tree))))
	no.samples = length(theta.cols)
	par(mfrow=c(1,no.samples),mar=c(2,4,2,4),cex.main=3)
	cols = rainbow(no.samples)
	for(t in 1:no.samples){
		tc=theta.cols[t]
		if(plotAsPercentage){
			tree$tree.label = paste(round(100*tree[,tc],1),"%",sep="")
		}else{
			tree$tree.label = tree[,tc]
		}
		tree$tree.label[!is.na(tree$annotation)] = paste(tree$tree.label[!is.na(tree$annotation)],"(",tree$annotation[!is.na(tree$annotation)],")",sep="")
		
		el<-as.matrix(tree[,c("ancestor","label")])
		
		gr <- graph.edgelist(el,directed=F)
		
		V(gr)$name = c("MRCA",tree$tree.label)
		
		V(gr)$color <- cols[t]
		V(gr)$label.font <- 2 #bold
		V(gr)$label.color <- "black"
		#edge.width and edge.color are not working
		V(gr)$edge.width <- 2
		if("col" %in% names(tree)){
			V(gr)$edge.color <- c("black",tree$col)
		}
		E(gr)$width <- 2
		E(gr)$color <- "black"
		lay=layout.sugiyama(gr,layers=depth,attributes="all")
		plot(gr,layout=lay$layout,vertex.size=36,main="",vertex.label.cex=font.size)
		if(length(main)==1){
			if(t==floor((no.samples+1)/2)){
				title(main,cex.main=font.size)
			}
		}else{
			title(main[t],cex.main=font.size)
		}
	}
}
