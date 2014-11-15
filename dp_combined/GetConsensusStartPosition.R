# # outdir = getwd()
# args=commandArgs(TRUE)
# run = as.integer(args[1])
# datpath = toString(args[2])
# purity_file = toString(args[3])
# identfile = toString(args[4])
# siblingfile = toString(args[5])
# ancestorfile = toString(args[6])
# 
# library(ggplot2)
# library(reshape2)
# 
# source("~/repo/dirichlet/dp_combined/LoadData.R")
# source("~/repo/dirichlet/dp_combined/interconvertMutationBurdens.R")
# 
# ###################### Set data variables ##########################
# # datpath = "/lustre/scratch110/sanger/sd11/dirichlet/simulated/Data/"
# #datafiles = c("simulated_2_sample1_3clusters_150muts_no_cna.txt", "simulated_2_sample2_3clusters_150muts_no_cna.txt", "simulated_2_sample3_3clusters_150muts_no_cna.txt")
# #cellularity = c(1, 1, 1)
# # identfile = "identity.strengths.csv"
# # siblingfile = "sibling.strengths.csv"
# # ancestorfile = "ancestor.strengths.csv"
# 
# sample2purity = read.table(purity_file, header=T, stringsAsFactors=F)
# samplename = unique(sample2purity$sample)[run]
# datafiles = sample2purity[sample2purity$sample==samplename,]$datafile
# subsamples = sample2purity[sample2purity$sample==samplename,]$subsample
# cellularity = sample2purity[sample2purity$sample==samplename,]$cellularity
# 
# read.in.data = function(filename) {
#   return(as.data.frame(read.table(filename, sep=",", header=T)))
# }

library(reshape2)
library(ggplot2)

transform.data = function(df) {
  #
  # A symmetric data.frame is mirrored through its diagonal
  # This method melts the df and then selects the lower half
  # which is returned in its melted form.
  #
  df = as.data.frame(df)  
  # Add the column names as an additional column
  df$id = factor(colnames(df), levels=colnames(df))
  # Melt the data creating colname-colname-value tripples
  df.m = melt(df)
  # As the input matrix is symmetric, select only the lower half
  df.lower = subset(df.m[lower.tri(df),], variable != id)
  return(df.lower)
}

createHeatmap = function(dat, x, y, value, xlab, ylab) {
  p = ggplot(dat) + 
    aes_string(x=x, y=y) + 
    geom_tile(aes_string(fill=value), colour="white") + 
    scale_fill_gradient(low="red", high="green") + 
    theme_bw() + 
    theme(axis.text.y = element_blank(), axis.text.x = element_blank())
  return(p)
}

createPng = function(p, filename, width, height) {
  png(filename=filename, width=width, height=height)
  print(p)
  dev.off()
}


get.best.cut = function(clustering, no.muts, min.frac=0.1) {
  #
  # Returns the best cut that maximises the number of nodes that have at least
  # the specified minimum fraction of mutations assigned to them.
  #
  best.cut = 1
  best.count = 1
  clusters = NULL
  for (i in 1:30) {
    tree.cut = cutree(clustering, k=i)
    unique.clusters = unique(tree.cut)
    selected.clusters = c()
    
    # Investigate the cut at this level
    cut.count = 0
    for (j in 1:length(unique.clusters)) {
      cluster = unique.clusters[j]
      cl.count = sum(tree.cut==cluster)
      if (cl.count > floor(min.frac*no.muts)) {
        cut.count = cut.count + 1
        selected.clusters = c(selected.clusters, j)
      }
    }
    
    # Keep if the current cut yields a larger number of nodes with more than 10% of mutations assigned to it
    if (cut.count > best.count) {
      best.cut = i
      best.count = cut.count
      clusters = selected.clusters
    }
  }
  return(list(best.cut=best.cut, clusters=clusters))
}

create.empty.tree = function(no.thetas) {
  #
  # Creates an empty tree with a column for each theta
  #
  theta.colnames = paste("theta.S",1:no.thetas, sep="")
  curr.tree = data.frame(t(rep(NA,no.thetas+4)))
  colnames(curr.tree) = c("label", theta.colnames, "ancestor", "annotation", "clusterid")
  curr.tree = curr.tree[-1,] # Remove the dummy NA row
  return(curr.tree)
}

add.node = function(curr.tree, clusterid, thetas, ancestorid) {
  if (nrow(curr.tree)==0) {
    # First node
    node = t(data.frame(c("M:", thetas, "Root:", NA, clusterid)))
    node.label = "M:"
    
  } else {
    ancestor = curr.tree[curr.tree$clusterid==ancestorid, ]
    anc.label = toString(ancestor$label)
    no.sibs = nrow(curr.tree[curr.tree$ancestor==anc.label,])
    node.label = paste(anc.label, no.sibs+1, ":", sep="")
    node = t(data.frame(c(node.label, thetas, anc.label, NA, clusterid)))
    
  }
  
  colnames(node) = colnames(curr.tree)
  rownames(node) = node.label
  curr.tree = rbind(curr.tree, node)
  
  return(curr.tree)
}

get.next.cluster = function(tree.cut, curr.node, clusters, sibs, anc) {
  #
  # Finds the next cluster that should be added to the tree by
  # looking at the sibling and ancestor agreements. The cluster
  # that has the highest agreement with the current node is returned.
  #
  
  # node : find most likely sibling
  agreement.sibs = 0
  most.likely.sib = -1
  for (i in 1:length(clusters)) {
    a = mean(rowSums(sibs[tree.cut==curr.node,tree.cut==clusters[i]]))
    if (a > agreement.sibs) {
      agreement.sibs = a
      most.likely.sib = clusters[i]
    }
  }
  
  # node : find most likely child
  agreement.anc = 0
  most.likely.anc = -1
  for (i in 1:length(clusters)) {
    a = mean(rowSums(anc[tree.cut==curr.node,tree.cut==clusters[i]]))
    if (a > agreement.anc) {
      agreement.anc = a
      most.likely.anc = clusters[i]
    }
  }
  
  if (agreement.anc >= agreement.sibs) {
    return(list(relation="anc", clusterid=most.likely.anc, strength=agreement.anc))
  } else {
    return(list(relation="sibs", clusterid=most.likely.sib, strength=agreement.sibs)) 
  }
}


## m=matrix(data=sample(rnorm(100,mean=0,sd=2)), ncol=10)
## this function makes a graphically appealing heatmap (no dendrogram) using ggplot
## whilst it contains fewer options than gplots::heatmap.2 I prefer its style and flexibility

# ggheat=function(m, rescaling='none', clustering='none', labCol=T, labRow=T, border=FALSE, 
#                 heatscale= c(low='blue',high='red'))
# {
#   ## the function can be be viewed as a two step process
#   ## 1. using the rehape package and other funcs the data is clustered, scaled, and reshaped
#   ## using simple options or by a user supplied function
#   ## 2. with the now resahped data the plot, the chosen labels and plot style are built
#   
#   require(reshape2)
#   require(ggplot2)
#   
#   ## you can either scale by row or column not both! 
#   ## if you wish to scale by both or use a differen scale method then simply supply a scale
#   ## function instead NB scale is a base funct
#   
#   if(is.function(rescaling))
#   { 
#     m=rescaling(m)
#   } 
#   else 
#   {
#     if(rescaling=='column') 
#       m=scale(m, center=T)
#     if(rescaling=='row') 
#       m=t(scale(t(m),center=T))
#   }
#   
#   ## I have supplied the default cluster and euclidean distance- and chose to cluster after scaling
#   ## if you want a different distance/cluster method-- or to cluster and then scale
#   ## then you can supply a custom function 
#   
#   if(is.function(clustering)) 
#   {
#     m=clustering(m)
#   }else
#   {
#     if(clustering=='row')
#       m=m[hclust(dist(m))$order, ]
#     if(clustering=='column')  
#       m=m[,hclust(dist(t(m)))$order]
#     if(clustering=='both') {
#       rows = hclust(dist(m))
#       cols = hclust(dist(t(m)))
#       m=m[rows$order ,cols$order]
#     }
#   }
#   ## this is just reshaping into a ggplot format matrix and making a ggplot layer
#   
#   rows=dim(m)[1]
#   cols=dim(m)[2]
#   melt.m=cbind(rowInd=rep(1:rows, times=cols), colInd=rep(1:cols, each=rows) ,melt(m))
#   g=ggplot(data=melt.m)
#   
#   ## add the heat tiles with or without a white border for clarity
#   
#   if(border==TRUE)
#     g2=g+geom_rect(aes(xmin=colInd-1,xmax=colInd,ymin=rowInd-1,ymax=rowInd, fill=value),colour='white')
#   if(border==FALSE)
#     g2=g+geom_rect(aes(xmin=colInd-1,xmax=colInd,ymin=rowInd-1,ymax=rowInd, fill=value))
#   
#   ## add axis labels either supplied or from the colnames rownames of the matrix
#   
#   if(labCol==T) 
#     g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=colnames(m))
#   if(labCol==F) 
#     g2=g2+scale_x_continuous(breaks=(1:cols)-0.5, labels=rep('',cols))
#   
#   if(labRow==T) 
#     g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rownames(m))	
#   if(labRow==F) 
#     g2=g2+scale_y_continuous(breaks=(1:rows)-0.5, labels=rep('',rows))	
#   
#   ## get rid of grey panel background and gridlines
#   
#   g2=g2+opts(panel.grid.minor=theme_line(colour=NA), panel.grid.major=theme_line(colour=NA),
#              panel.background=theme_rect(fill=NA, colour=NA))
#   
#   ## finally add the fill colour ramp of your choice (default is blue to red)-- and return
#   return(g2+scale_fill_continuous("", heatscale[1], heatscale[2]))
#   
# }

## NB because ggheat returns an ordinary ggplot you can add ggplot tweaks post-production e.g. 
## data(mtcars)
## x= as.matrix(mtcars)
## ggheat(x, clustCol=T)+ opts(panel.background=theme_rect(fill='pink'))
# ggheat(x, clustering='both', rescaling='row', heatscale=c(low='red', high='yellow'))
# 
# 
# 
# for (infile in c(identfile, siblingfile, ancestorfile)) {
#   dat = read.in.data(infile)
#   dat.m = transform.data(dat)
#   p = createHeatmap(dat.m, x="id", y="variable", value="value",xlab="Mutation", ylab="Mutation")
#   createPng(p, filename=paste(infile, ".png", sep=""), width=1000, height=1000)
# }




find.cons.start.pos = function(subclonal.fraction, ident, sibs, anc, min.frac) {
  #
  # Finds a start position for GetConsensusTrees. This allows the EM there to
  # not start from all mutations assigned to one node (i.e. fully clonal), but
  # already starting from a layout that is somewhat in the right direction. This
  # could help GetConsensusTrees to converge faster.
  #
  # ident, sibs and anc are the identity, sibling and ancestor strengths
  # min.frac is the minimum fraction of mutations that need to be assigned to a
  # cluster of mutations in order to consider it as a node on the starting tree
  #
#   if(0) {
  ###################### Cluster ##########################
  print("Clustering")
  rows = hclust(dist(ident))
  cols = hclust(dist(t(ident)))
  save(file="find.cons.start.tree_clustering.RData", rows, cols)

  ###################### Create figures ##########################
  print("Creating figures")
  dat.m = transform.data(ident[rows$order, cols$order])
  p = createHeatmap(dat.m, x="id", y="variable", value="value",xlab="Mutation", ylab="Mutation")
  createPng(p, filename=paste("find.cons.start.tree_ident", ".png", sep=""), width=1000, height=1000)

  dat.m = transform.data(sibs[rows$order, cols$order])
  p = createHeatmap(dat.m, x="id", y="variable", value="value",xlab="Mutation", ylab="Mutation")
  createPng(p, filename=paste("find.cons.start.tree_sibs", ".png", sep=""), width=1000, height=1000)

  dat.m = transform.data(anc[rows$order, cols$order])
  p = createHeatmap(dat.m, x="id", y="variable", value="value",xlab="Mutation", ylab="Mutation")
  createPng(p, filename=paste("find.cons.start.tree_anc", ".png", sep=""), width=1000, height=1000)
#   }
#   load("getconsensus_clustering.RData")
  ###################### Obtain best cut ##########################
  print("Getting best cut")
  res = get.best.cut(rows, nrow(ident), min.frac=min.frac)
  best.cut = res$best.cut
  clusters = res$clusters
  tree.cut = cutree(rows, k=best.cut)
  
  ###################### Obtain theta estimates from subclonal.fractions ##########################
  #clusters = unique(tree.cut)
  thetas = list()
  for (i in 1:length(clusters)) {
    cluster = clusters[i]
    if (cluster %in% selected.clusters) { # A cluster is in selected.clusters if it contains more than min.frac mutations
      print(i)
      print(cluster)
      print(subclonal.fraction[tree.cut==cluster,])
      print(sum(tree.cut==cluster))

      # TODO: Perhaps replace this by a draw from a gamma distribution?
      thetas[[i]] = colSums(subclonal.fraction[tree.cut==cluster,])/sum(tree.cut==cluster)
    }
  }
  
  ###################### Find first node ##########################
  print("Finding first node")
  thetas.sum = unlist(lapply(thetas, sum))
  first.node = which(thetas.sum==max(thetas.sum))
  no.thetas = ncol(subclonal.fraction)
  
  # Create a new tree and add the node to the tree
  curr.tree = create.empty.tree(no.thetas)
  curr.tree = add.node(curr.tree, first.node, thetas[[first.node]], NA)
  
  ###################### Iteratively add most likely node ##########################
  print("Adding clusters as nodes")
  # Add the first node to the tree
#   prev.node = first.node
  clusters = clusters[clusters!=first.node]
  
  # Iteratively grow the tree by finding the node with the strongest relation to any node
  # and add it to the tree. Once added the cluster is removed from the list. Stop once all 
  # clusters are on the tree.
  #
  # We're using integers here to keep track of nodes and clusters. Clusters are assigned
  # an index (the order in which they appear in the clusters variable). We keep this index
  # when the cluster is added to the node. Therefore clusterid == nodeid.
  while(length(clusters) > 0) {
    best.strength = 0
    best.node = 0
    best.cluster = NULL
    print("Searching for next cluster")
    for (nodeid in curr.tree$clusterid) { # iterate over all nodes in the tree
      print(paste("Looking at node ",nodeid))
      next.cluster = get.next.cluster(tree.cut, nodeid, clusters, sibs, anc)
      if (next.cluster$strength > best.strength) { # if this cluster is better than best, keep it
        best.strength = next.cluster$strength
        best.cluster = next.cluster
        best.node = nodeid # save the id of the node that this cluster has the strongest strength with
        
        print("Next best cluster")
        print(best.cluster)
        print(best.node)
      }
    }
    
    # Add the cluster as node to the tree
    if (next.cluster$relation == "anc") {
      # Direct ancestor, best.node is the ancestor
      curr.tree = add.node(curr.tree, best.cluster$clusterid, thetas[[best.cluster$clusterid]], best.node)
    } else {
      # Sibling, ancestor of best.node is the ancestor
      curr.tree = add.node(curr.tree, best.cluster$clusterid, thetas[[best.cluster$clusterid]], as.numeric(curr.tree[curr.tree$clusterid==best.node,]$clusterid))
    }
    # Remove the cluster from the list
    clusters = clusters[clusters!=best.cluster$clusterid]
#     prev.node = best.cluster$clusterid
  }
  # Set mutation assignments
  assignments = curr.tree$label[tree.cut]
  
  ###################### Save results ##########################
  print(curr.tree)
  write.table(file="find.cons.start.tree.txt", curr.tree, row.names=F, sep="\t", quote=F)
  write.table(file="find.cons.start.assignments.txt", assignments, row.names=F, sep="\t", quote=F)
  return(list(tree=curr.tree, assignments=assignments)) 
}


# ###################### Load the data ##########################
# dataset = load.data(datpath,
#                     "",
#                     datafiles, 
#                     cellularity=cellularity, 
#                     Chromosome="chr", 
#                     WT.count="WT.count", 
#                     mut.count="mut.count", 
#                     subclonal.CN="subclonal.CN", 
#                     no.chrs.bearing.mut="no.chrs.bearing.mut", 
#                     mutation.copy.number="mutation.copy.number", 
#                     subclonal.fraction="subclonal.fraction", 
#                     data_file_suffix="")
# 
# ident = read.in.data(identfile)
# sibs = read.in.data(siblingfile)
# anc = read.in.data(ancestorfile)
# 
# find.cons.start.pos(dataset$subclonal.fraction, ident, sibs, anc, 0.15)
# 
# q(save="no")
