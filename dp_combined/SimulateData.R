source("interconvertMutationBurdens.R")
source("PlotTreeWithIgraph.R")
#
# Testcases
# - getSimpleSubclones with 0.5, 0.3 and 0.1 should yield 3 clusters with allele frequencies around those numbers
# 
#
#
#
#

# Datasets generated:
# cellularity = c(1,1,1)
# no.muts = 50
# sim = simulate.muts(100, 50, 3, mut.frac.of.cells=c(1,1,1))
# ds = sim.muts2dataset(sim, 100, cellularity)
# writeDataset("", "simulated_1", c("1", "2", "3"), "dp_input.txt", ds, cellularity)
# 
# ds2 = sim.muts2dataset(simulate.muts(25, 50, 3, mut.frac.of.cells=c(1,1,0)), 25, cellularity)
# ds3 = sim.muts2dataset(simulate.muts(25, 50, 3, mut.frac.of.cells=c(1,0,1)), 25, cellularity)
# combined = appendSubcloneData(list(ds, ds2, ds3))
# writeDataset("", "simulated_2", c("1", "2", "3"), "dp_input.txt", combined, cellularity)

# Not used
# sim1 = simulate.muts(5, 50, 3, mut.frac.of.cells=c(1,1,1), mut.cn=c(0,2,3))
# sim2 = simulate.muts(5, 50, 3, mut.frac.of.cells=c(1,1,1), mut.cn=c(3,0,2))
# sim3 = simulate.muts(5, 50, 3, mut.frac.of.cells=c(1,1,1), mut.cn=c(2,3,0))
# sim4 = simulate.muts(5, 50, 3, mut.frac.of.cells=c(1,1,1), wt.cn=c(2,3,0))
# sim5 = simulate.muts(5, 50, 3, mut.frac.of.cells=c(1,1,1), wt.cn=c(0,2,3))
# sim6 = simulate.muts(5, 50, 3, mut.frac.of.cells=c(1,1,1), wt.cn=c(3,0,2))
# sim7 = simulate.muts(5, 50, 3, mut.frac.of.cells=c(1,1,1), wt.subcl.cn=c(0.3,0,0))
# sim8 = simulate.muts(5, 50, 3, mut.frac.of.cells=c(1,1,1), wt.subcl.cn=c(0,0.1,0.1))
# sim9 = simulate.muts(5, 50, 3, mut.frac.of.cells=c(1,1,1), mut.cn=c(0,0.6,0.2))
# sim10 = simulate.muts(5, 50, 3, mut.frac.of.cells=c(1,1,1), mut.cn=c(0.1,0,0.4))
# 
# mut.cn = matrix(c(rpois(15*3,1.5), rep(1, 35*3)), ncol=3, byrow=T)
# mut.cn[,1] = mut.cn[,1]+1
# 
# mut.cn = matrix(c(0,2,3, 3,0,2, 2,3,0, rep(1, ))
# 
# ds2 = simulate.subclone(c(1,1,1), c(50,50,50), 3, mut.cn=mut.cn



# simulate.subclone = function(no.muts, depth.per.cn, no.subsamples, mut.cn=matrix(1, ncol=no.subsamples, nrow=sum(no.muts)), mut.subcl.cn=matrix(1, ncol=no.subsamples, nrow=sum(no.muts)), wt.cn=matrix(1, ncol=no.subsamples, nrow=sum(no.muts)), wt.subcl.cn=matrix(1, ncol=no.subsamples, nrow=sum(no.muts)), mut.frac.of.cells=matrix(1, ncol=no.subsamples, nrow=sum(no.muts)), cellularity=rep(1, no.subsamples), mut.align.bias=1) {
#   simmed.ds = list()
#   for (i in 1:sum(no.muts)) {
#     # simulate a mutation with the given specs and transform it into a dataset
#     simmed.ds[[i]] = sim.muts2dataset(simulate.muts(1, 
#                                                     depth.per.cn[i], 
#                                                     no.subsamples,
#                                                     mut.cn=mut.cn[i,],
#                                                     mut.subcl.cn=mut.subcl.cn[i,],
#                                                     wt.cn=wt.cn[i,],
#                                                     wt.subcl.cn=wt.subcl.cn[i,],
#                                                     mut.frac.of.cells=mut.frac.of.cells[i,],
#                                                     cellularity=cellularity,
#                                                     mut.align.bias=mut.align.bias), 
#                                       1, 
#                                       cellularity)
#   }
#   # Concatenate all the single mutation datasets
#   ds = appendSubcloneData(simmed.ds)
#   return(ds)
# }


################################### Functions that generate the various sample sets ###################################
create.testing.samples = function() {
  no.subsamples = 3
  #if(0) {
  # Affected by nothing, 3 simple subclones across 3 samples
  ds = sim.dataset(no.subclones=3, 
                    no.muts=c(100,25,25), 
                    no.subsamples=no.subsamples, 
                    cov=c(50,50,50), 
                    mut.frac.of.cells=matrix(c(c(1,0.5,0.3), c(1,0,1), c(1,1,0)), ncol=no.subsamples),
                    sample.cn.param=c(-1,-1,-1), # pois parameter that affects CN values per sample
                    sample.cn.frac=c(0,0,0), # frac of muts affected by CN per sample
                    sample.cn.del.frac=c(0,0,0), # frac of muts that are affected by a deletion
                    sample.subcl.cn.frac=c(0,0,0), # frac of muts affected by subclonal CN per sample
                    cellularity=rep(1, no.subsamples)) #rnorm(3,mean=0.7,sd=0.2))
  writeDataset("", "simulated_001", c("01", "02", "03"), "dp_input.txt", ds, ds$cellularity)

  # Affected by clonal copynumber
  ds = sim.dataset(no.subclones=3, 
                    no.muts=c(100,25,25), 
                    no.subsamples=no.subsamples, 
                    cov=c(50,50,50), 
                    mut.frac.of.cells=matrix(c(c(1,0.5,0.3), c(1,0,1), c(1,1,0)), ncol=no.subsamples),
                    sample.cn.param=c(1.45,1.55,1.5), # pois parameter that affects CN values per sample
                    sample.cn.frac=c(0.25,0.25,0.25), # frac of muts affected by CN per sample
                    sample.cn.del.frac=c(0,0,0), # frac of muts that are affected by a deletion
                    sample.subcl.cn.frac=c(0,0,0), # frac of muts affected by subclonal CN per sample
                    cellularity=rep(1, no.subsamples)) #rnorm(3,mean=0.7,sd=0.2))
  writeDataset("", "simulated_002", c("01", "02", "03"), "dp_input.txt", ds, ds$cellularity)
  
  # Affected by clonal delections
  ds = sim.dataset(no.subclones=3, 
                    no.muts=c(100,25,25), 
                    no.subsamples=no.subsamples, 
                    cov=c(50,50,50), 
                    mut.frac.of.cells=matrix(c(c(1,0.5,0.3), c(1,0,1), c(1,1,0)), ncol=no.subsamples),
                    sample.cn.param=c(-1,-1,-1), # pois parameter that affects CN values per sample
                    sample.cn.frac=c(0,0,0), # frac of muts affected by CN per sample
                    sample.cn.del.frac=c(0.09,0.11,0.1), # frac of muts that are affected by a deletion
                    sample.subcl.cn.frac=c(0,0,0), # frac of muts affected by subclonal CN per sample
                    cellularity=rep(1, no.subsamples)) #rnorm(3,mean=0.7,sd=0.2))
  writeDataset("", "simulated_003", c("01", "02", "03"), "dp_input.txt", ds, ds$cellularity)
  
  # Affected by sublonal copynumber
  ds = sim.dataset(no.subclones=3, 
                    no.muts=c(100,25,25), 
                    no.subsamples=no.subsamples, 
                    cov=c(50,50,50), 
                    mut.frac.of.cells=matrix(c(c(1,0.5,0.3), c(1,0,1), c(1,1,0)), ncol=no.subsamples),
                    sample.cn.param=c(-1,-1,-1), # pois parameter that affects CN values per sample
                    sample.cn.frac=c(0,0,0), # frac of muts affected by CN per sample
                    sample.cn.del.frac=c(0,0,0), # frac of muts that are affected by a deletion
                    sample.subcl.cn.frac=c(0.11,0.09,0.1), # frac of muts affected by subclonal CN per sample
                    cellularity=rep(1, no.subsamples)) #rnorm(3,mean=0.7,sd=0.2))
  writeDataset("", "simulated_004", c("01", "02", "03"), "dp_input.txt", ds, ds$cellularity)
  
  # Affected by lower cellularity
  cellularity = rnorm(no.subsamples,mean=0.7,sd=0.2)
  cellularity[cellularity > 1] = 1
  ds = sim.dataset(no.subclones=3, 
                    no.muts=c(100,25,25), 
                    no.subsamples=no.subsamples, 
                    cov=c(50,50,50), 
                    mut.frac.of.cells=matrix(c(c(1,0.5,0.3), c(1,0,1), c(1,1,0)), ncol=no.subsamples),
                    sample.cn.param=c(-1,-1,-1), # pois parameter that affects CN values per sample
                    sample.cn.frac=c(0,0,0), # frac of muts affected by CN per sample
                    sample.cn.del.frac=c(0,0,0), # frac of muts that are affected by a deletion
                    sample.subcl.cn.frac=c(0,0,0), # frac of muts affected by subclonal CN per sample
                    cellularity=cellularity)
  writeDataset("", "simulated_005", c("01", "02", "03"), "dp_input.txt", ds, ds$cellularity)
  
  # All above combined
  cellularity = rnorm(no.subsamples,mean=0.7,sd=0.2)
  cellularity[cellularity > 1] = 1
  ds = sim.dataset(no.subclones=3, 
                    no.muts=c(100,25,25), 
                    no.subsamples=no.subsamples, 
                    cov=c(50,50,50), 
                    mut.frac.of.cells=matrix(c(c(1,0.5,0.3), c(1,0,1), c(1,1,0)), ncol=no.subsamples),
                    sample.cn.param=c(1.45,1.55,1.5), # pois parameter that affects CN values per sample
                    sample.cn.frac=c(0.25,0.25,0.25), # frac of muts affected by CN per sample
                    sample.cn.del.frac=c(0.09,0.11,0.1), # frac of muts that are affected by a deletion
                    sample.subcl.cn.frac=c(0.11,0.09,0.1), # frac of muts affected by subclonal CN per sample
                    cellularity=cellularity)
  writeDataset("", "simulated_006", c("01", "02", "03"), "dp_input.txt", ds, ds$cellularity)

  # This sample was created to test a few things, but it turned out to be a fundamental bug that affected all samples.
  # They therefore needed to be regenerated anyway. This sample is the same as 006.
  cellularity = rnorm(no.subsamples,mean=0.7,sd=0.2)
  cellularity[cellularity > 1] = 1
  ds = sim.dataset(no.subclones=3, 
                    no.muts=c(100,25,25), 
                    no.subsamples=no.subsamples, 
                    cov=c(50,50,50), 
                    mut.frac.of.cells=matrix(c(c(1,0.5,0.3), c(1,0,1), c(1,1,0)), ncol=no.subsamples),
                    sample.cn.param=c(1.45,1.55,1.5), # pois parameter that affects CN values per sample
                    sample.cn.frac=c(0.25,0.25,0.25), # frac of muts affected by CN per sample
                    sample.cn.del.frac=c(0.09,0.11,0.1), # frac of muts that are affected by a deletion
                    sample.subcl.cn.frac=c(0.11,0.09,0.1), # frac of muts affected by subclonal CN per sample
                    cellularity=cellularity)
  writeDataset("", "simulated_007", c("01", "02", "03"), "dp_input.txt", ds, ds$cellularity)
#   }
}



create.random.trees.small = function() {
#   trees = list()
#   for (i in 1:250) { trees[[i]] = sim.tree() }
#   no.nodes = unlist(lapply(trees, nrow))
#   no.subsamples = unlist(lapply(trees, FUN=function(x) { sum(grepl("theta", colnames(x))) }))
#   hist(no.subsamples)
#   hist(no.nodes)
  allowed.no.muts = rep(c(25,50,100,250,500,1000), 10) #,2500,5000,7500,10000,15000,20000,25000,30000)
  
  
  #
  # CN modelled as a poisson distribution. Lambda parameters learned from Prostate data:
  #       > print(summary(unlist(lambda)))
  #         Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #         0.9618  1.1020  1.3000  1.3800  1.4440  2.7860  
  #
  
  trees = list()
  for (i in 1:10) { trees[[i]] = sim.tree() }

  datasets = list()
  for (i in 1:length(trees)) {
    
    tree = trees[[i]]
    no.subclones = nrow(tree)
    no.subsamples = sum(grepl("theta", colnames(tree)))
    cov = rep(50, no.subsamples)
    no.muts = sample(allowed.no.muts)[1:no.subclones]
    mut.frac.of.cells = tree[,grepl("theta", colnames(tree))]
    print(mut.frac.of.cells)
    if (!is.data.frame(mut.frac.of.cells)) { mut.frac.of.cells = as.data.frame(matrix(mut.frac.of.cells, ncol=1)) }
    print(mut.frac.of.cells)
    sample.cn.param = rep(1.3, no.subsamples) # Average lambda taken from prostate samples
    sample.cn.frac = rep(0.55, no.subsamples) # Average taken from prostate samples
    sample.cn.del.frac = rep(0.05, no.subsamples) # Average taken from prostate samples
    sample.subcl.cn.frac = rep(0.14, no.subsamples) # Average taken from prostate samples
    cellularity = rnorm(no.subsamples,mean=0.7,sd=0.2)
    cellularity[cellularity > 1] = 1    
    
    datasets[[i]] = sim.dataset(no.subclones=no.subclones, no.muts=no.muts, no.subsamples=no.subsamples, cov=cov, mut.frac.of.cells=mut.frac.of.cells, sample.cn.param=sample.cn.param, sample.cn.frac=sample.cn.frac, sample.cn.del.frac=sample.cn.del.frac, sample.subcl.cn.frac=sample.subcl.cn.frac, cellularity=cellularity)
    
  }

  
  writeDataset("", "simulated_a", "simulated_a", datasets, trees)
}

################################### Master functions ###################################
sim.tree = function(no.subsamples=NULL, max.subsamples=10, max.nodes=30, min.node.theta=0.01, min.theta.gap=0.05, max.tries=20) {
  #
  # Simulate a full tree
  #
  
  create.node = function(label, parent.label, node.number, thetas) {
    # Create empty dataframe in right format
    no.subsamples = length(thetas)
    new.node = as.data.frame(t(data.frame(rep(NA, 3+no.subsamples))))
    colnames(new.node) = c("label", "ancestor", "node", paste("theta.S", 1:no.subsamples, sep=""))
    row.names(new.node) = label
    # Fill the various columns, making sure theta's are numeric
    theta.cols = grepl("theta", colnames(new.node))
    new.node[1,] = c(label, as.character(parent.label), node.number, thetas)
    new.node[1,1:3] = as.character(new.node[1,1:3])
    new.node[,theta.cols] = as.numeric(new.node[,theta.cols])
    return(new.node)
  }
  
  
  # Choose number of sequencing samples randomly
  if (is.null(no.subsamples)) { no.subsamples = round(runif(1, 1, max.subsamples)) }
  
  # Create a tree consisting of 1 to 30 nodes
  no.nodes = round(runif(1, 1, max.nodes))
  
  # Create data.frame with the first node
  tree = create.node("M:", "root", 1, rep(1, no.subsamples))
  theta.cols = grepl("theta", colnames(tree)) # save which columns have thetas for later use

  for (i in 2:no.nodes) {
    # Add nodes
    not.node.added = T
    count = 1
    while(not.node.added) {
      if (count <= max.tries) {
        count = count + 1
        
        # 1: randomly draw level
        tree.levels = unlist(lapply(as.character(tree$label), nchar))/2
        level = round(runif(1, 1, max(tree.levels)))
        
        # 2: randomly draw branch, find parent and children
        nodes.at.level = which(tree.levels == level)
        branch = round(runif(1, 1, length(nodes.at.level)))
        parent = tree[nodes.at.level[branch],]
        children = tree[tree$ancestor %in% parent$label,]
        
        # 3: randomly draw thetas, at least MIN.THETA.GAP away from the parent
        max.new.theta = parent[,theta.cols]
        if (nrow(children) > 0) { # If parent already has children, make sure children combined do not add up to more than parent
          if (no.subsamples == 1) {
            max.new.theta = max.new.theta - children[,theta.cols]
          } else {
            max.new.theta = max.new.theta - colSums(children[,theta.cols])
          }
        }
        # Draw a new number between the minimum theta and the parent theta
        if (all(max.new.theta-min.theta.gap > min.node.theta)) {
          new.node.theta = runif(no.subsamples, rep(min.node.theta, no.subsamples), as.numeric(max.new.theta-min.theta.gap))
          not.node.added = F
        
          # 4: create node and add it
          label = paste(parent$label, nrow(children)+1, ":", sep="")
          new.node = create.node(label, parent$label, nrow(tree)+1, new.node.theta)
          tree = rbind(tree,new.node)
        }
      } else {
        warning("Reached max iterations, no tree to be found")
        return(tree)
      }
    }
  }
  return(tree)
}


sim.dataset = function(no.subclones, no.muts, no.subsamples, cov=rep(50, no.subsamples), mut.frac.of.cells=matrix(rep(1, no.subclones*no.subsamples), ncol=no.subsamples), sample.cn.param=rep(-1, no.subsamples), sample.cn.frac=rep(0, no.subsamples), sample.cn.del.frac=rep(0, no.subsamples), sample.subcl.cn.frac=rep(0, no.subsamples), cellularity=rep(1, no.subsamples)) {
  #
  # Simulate a dataset by iteratively creatig subclones for all subsamples.
  # Returns a list that contains all the various bits of information, including mutCount, WTcount, mutation.copy.number, etc.
  # Also stored are cellularity, expected coverage and the copynumber parameters.
  #
  dataset = list()
  for (i in 1:no.subclones) {
    subclone = simulate.subclone(no.muts[i], 
                                 no.subsamples=no.subsamples, 
                                 cov=cov, 
                                 mut.frac.of.cells=mut.frac.of.cells[i,],
                                 mut.cn.lambda=sample.cn.param,
                                 wt.cn.lambda=sample.cn.param,
                                 sample.cn.frac=sample.cn.frac, 
                                 sample.cn.del.frac=sample.cn.del.frac, 
                                 sample.subcl.cn.frac=sample.subcl.cn.frac,
                                 cellularity=cellularity)

    dataset[[i]] = mergeColumns(subclone)
    dataset[[i]]$subcloneid = matrix(rep(i, no.muts[i]), ncol=1)
  }
  
  dataset = appendSubcloneData(dataset)
  dataset$cellularity = cellularity
  dataset$cov = cov
  dataset$mut.frac.of.cells = mut.frac.of.cells
  dataset$sample.cn.param = sample.cn.param
  dataset$sample.cn.frac = sample.cn.frac
  dataset$sample.cn.del.frac = sample.cn.del.frac
  
  return(dataset)
}


################################### Helper functions ###################################
simulate.subclone = function(no.muts, no.subsamples, cov, mut.frac.of.cells, mut.cn.lambda=rep(-1, no.subsamples), wt.cn.lambda=rep(-1, no.subsamples), cellularity=rep(1, no.subsamples), sample.cn.frac=rep(0, no.subsamples), sample.cn.del.frac=rep(0, no.subsamples), sample.subcl.cn.frac=rep(0, no.subsamples)) { #.sd=0.2
  
  create.cn = function(no.muts, cn.lambda, cn.frac, del.frac) {
    if (cn.lambda > -1) { 
      # Draw copynumber from pois distribution+1
      cn = rpois(no.muts, cn.lambda)+1
      
      # CN should affect a fixed percentage of these mutations
      sel = rbinom(no.muts, 1, 1-cn.frac)==1
      cn[sel] = 1 # Reset to 1
      
      # A set percentage should be affected by a loss
      sel = rbinom(no.muts, 1, del.frac)==1
      cn[sel] = 0 # Reset to 0
      
    } else {
      # No CN changes
      cn = rep(1, no.muts)
    }
    return(cn)
  }
  
  create.subcl.cn = function(no.muts, subcl.cn.frac) {
    if (subcl.cn.frac > 0) {
      # Draw subclonal copynumber from uniform [0,1]
      subcl.cn = runif(no.muts)
      
      # subclonal CN should affect a fixed percentage of these mutations
      sel = rbinom(no.muts, 1, 1-subcl.cn.frac)==1
      subcl.cn[sel] = 0 # Reset to 0
      
    } else {
      # No subclonal CN
      subcl.cn = rep(0, no.muts)
    }
    return(subcl.cn)
  }
  
  # For each subclone, create a CN profile, draw subclonal CN, then generate the mutations and finally morph this into dataset format
  subclone = list()
  for (i in 1:no.subsamples) {
    mut.cn = create.cn(no.muts, mut.cn.lambda[i], sample.cn.frac[i], sample.cn.del.frac[i])
    mut.subcl.cn = create.subcl.cn(no.muts, sample.subcl.cn.frac[i])
    wt.cn = create.cn(no.muts, wt.cn.lambda[i], sample.cn.frac[i], sample.cn.del.frac[i])
    wt.subcl.cn = create.subcl.cn(no.muts, sample.subcl.cn.frac[i])
    
    muts = generate.mut(no.muts=no.muts, depth.per.cn=cov[i], mut.cn=mut.cn, mut.subcl.cn=mut.subcl.cn, wt.cn=wt.cn, wt.subcl.cn=wt.subcl.cn, mut.frac.of.cells=as.numeric(mut.frac.of.cells[i]), cellularity=cellularity[i], mut.align.bias=1)
    subclone[[i]] = sim.muts2dataset(muts, no.muts, cellularity[i])
  }
  return(subclone)
}

generate.mut = function(no.muts, depth.per.cn, mut.cn, mut.subcl.cn, wt.cn, wt.subcl.cn, mut.frac.of.cells, cellularity, mut.align.bias) {
  # Tumour cells carying mut
  mutCount = cellularity*(mut.cn+mut.subcl.cn)*mut.frac.of.cells*depth.per.cn*mut.align.bias
  # Tumour cells carying WT - 1+1-mut.frac.of.cells here accounts for 1*tumour cells WT and (1-mut.frac.of.cells)*tumour cells not mut (assuming heterozygous muts)
  WTcount.1 = cellularity*(wt.cn+wt.subcl.cn)*(1+1-mut.frac.of.cells)*depth.per.cn*(2-mut.align.bias)
  # Stromal cells
  WTcount.2 = rep((1-cellularity)*2*depth.per.cn*(2-mut.align.bias), no.muts) #(1-mut.frac.of.cells)*
  
  WTcount = WTcount.1 + WTcount.2
  WTcount[WTcount==0] = 1 # Setting WTcount to 1 in order to obtain a finite, but very small AF below 
  
  totalCopyNumber = mut.cn+mut.subcl.cn+wt.cn+wt.subcl.cn
  copyNumberAdjustment = mut.cn #+mut.subcl.cn #add mut.subcl.cn here to achieve cleaner data. without it some mutations are in too high fraction of cells
  # Depth takes on a poisson distribution
  ndepth = rpois(no.muts,mutCount+WTcount)
  # mutCount a binomial
  nmut = rbinom(no.muts, round(ndepth), mutCount/(mutCount+WTcount))
  nmut[ndepth==0] = 0
  
  return(list(mutCount=nmut, depth=ndepth, true.mutCount=mutCount, true.WTcount=WTcount, totalCopyNumber=totalCopyNumber, copyNumberAdjustment=copyNumberAdjustment))
}


sim.muts2dataset = function(sim, no.muts, cellularity) {
  #
  # Converts simulated mutations into a dataset.
  #
  no.subsamples = ncol(sim$mutCount)
  dataset = list()
  dataset$mutCount = sim$mutCount
  dataset$WTCount = sim$depth-sim$mutCount
  dataset$WTCount[dataset$WTCount < 0] = 0
  dataset$copyNumberAdjustment = sim$copyNumberAdjustment #matrix(rep(1, no.subsamples*no.muts), ncol=no.subsamples)
  dataset$totalCopyNumber = sim$totalCopyNumber #matrix(rep(2, no.subsamples*no.muts), ncol=no.subsamples)
  dataset$mutation.copy.number = mutationBurdenToMutationCopyNumber(sim$mutCount/(sim$depth), sim$totalCopyNumber, cellularity, rep(2,length(sim$mutCount)))
  dataset$kappa = mutationCopyNumberToMutationBurden(1,dataset$totalCopyNumber,cellularity) * dataset$copyNumberAdjustment
  dataset$subclonal.fraction = sim$mutCount / (sim$depth*dataset$kappa)
  dataset$subclonal.fraction[dataset$kappa == 0] = 0
  return(dataset)
}

appendSubcloneData <- function(list_of_subclone_datasets) {
  #
  # rbinds a list of datasets together into one
  #
  combined = do.call(cbind, list_of_subclone_datasets)
  mutCount = do.call(rbind, combined['mutCount',])
  WTCount = do.call(rbind, combined['WTCount',])
  copyNumberAdjustment = do.call(rbind, combined['copyNumberAdjustment',])
  totalCopyNumber = do.call(rbind, combined['totalCopyNumber',])
  mutation.copy.number = do.call(rbind, combined['mutation.copy.number',])
  kappa = do.call(rbind, combined['kappa',])
  subclonal.fraction = do.call(rbind, combined['subclonal.fraction',])
  subcloneid = do.call(rbind, combined['subcloneid',])
  return(list(mutCount=mutCount, WTCount=WTCount, kappa=kappa, copyNumberAdjustment=copyNumberAdjustment, totalCopyNumber=totalCopyNumber, mutation.copy.number=mutation.copy.number, subclonal.fraction=subclonal.fraction, subcloneid=subcloneid))
}

mergeColumns <- function(list_of_datasets) {
  #
  # cbinds a list of datasets together into one
  #
  combined = do.call(cbind, list_of_datasets)
  mutCount = do.call(cbind, combined['mutCount',])
  WTCount = do.call(cbind, combined['WTCount',])
  copyNumberAdjustment = do.call(cbind, combined['copyNumberAdjustment',])
  totalCopyNumber = do.call(cbind, combined['totalCopyNumber',])
  mutation.copy.number = do.call(cbind, combined['mutation.copy.number',])
  kappa = do.call(cbind, combined['kappa',])
  subclonal.fraction = do.call(cbind, combined['subclonal.fraction',])
  return(list(mutCount=mutCount, WTCount=WTCount, kappa=kappa, copyNumberAdjustment=copyNumberAdjustment, totalCopyNumber=totalCopyNumber, mutation.copy.number=mutation.copy.number, subclonal.fraction=subclonal.fraction, subcloneid=list_of_datasets[[1]]$subcloneid))
}

simulatedDataToInputFile <- function(dataset, outfile_names) {
  no.muts = nrow(dataset$mutCount)
  no.subsamples = ncol(dataset$mutCount)
  column_names = c("chr","pos","WT.count", "mut.count", "subclonal.CN","mutation.copy.number","subclonal.fraction","no.chrs.bearing.mut","subclone.id")
  
  for (i in 1:no.subsamples) {
    dat = cbind(array(rep(1,no.muts)),
                array(1:no.muts),
                dataset$WTCount[,i],
                dataset$mutCount[,i],
                dataset$totalCopyNumber[,i],
                dataset$mutation.copy.number[,i],
                dataset$subclonal.fraction[,i],
                dataset$copyNumberAdjustment[,i],
                dataset$subcloneid[,1])

    colnames(dat) = column_names
    write.table(dat,file=outfile_names[i],quote=F,row.names=F, sep="\t")
  }
}


writeDataset = function(outpath, datasetname, samplename_prefix, datasets, trees) {
  
  dataset_index = data.frame()
  
  for (i in 1:length(datasets)) {
    ds = datasets[[i]]
    tree = trees[[i]]
    
    # Collect various parameters for later use
    no.subsamples = sum(grepl("theta", colnames(tree)))
    no.clusters = nrow(tree)
    no.muts = nrow(ds$mutCount)
    
    # Create a samplename where the id is prefixed with zeros when appropriate for properly sorted files on disk
    sampleid = formatC(i, width=4, format="d", flag="0")
    samplename = paste(samplename_prefix, sampleid, sep="_")
    
    # Create subsamplenames where the ids are prefixed with zeros when appropriate for properly sorted files on disk
    subsampleids = formatC(1:no.subsamples, width=2, format="d", flag="0") 
    subsamplenames = paste(samplename, subsampleids, sep="_")
    datafiles = paste(subsamplenames, ".txt", sep="")
    
    # Collapse the thetas into a comma separated list per subsample
    mut.frac.of.cells = c()
    for (i in 1:ncol(ds$mut.frac.of.cells)) {
      mut.frac.of.cells = c(mut.frac.of.cells, paste(ds$mut.frac.of.cells[,i],sep="", collapse=","))
    }

    # Create index for this dataset
    dataset_index_sam = data.frame(sample=rep(samplename, no.subsamples), 
                                   subsample=1:no.subsamples, 
                                   datafile=datafiles, 
                                   cellularity=ds$cellularity, 
                                   no.muts=rep(no.muts, no.subsamples),
                                   no.clusters=rep(no.clusters, no.subsamples),
                                   coverage=ds$cov, 
                                   mut.frac.of.cells=mut.frac.of.cells,
                                   sample.cn.param=ds$sample.cn.param,
                                   sample.cn.frac=ds$sample.cn.frac,
                                   sample.cn.del.frac=ds$sample.cn.del.frac)

    dataset_index = rbind(dataset_index, dataset_index_sam)
    
    # Write the subsamples to disk
    simulatedDataToInputFile(ds, datafiles)
    
    # Write the tree
    write.table(tree, file=paste(outpath, samplename, ".tree.txt", sep=""),quote=F,row.names=F, sep="\t")
    
    # Create a png of the tree
    row.names(tree) = tree$label
    png(paste(outpath, samplename, "_expected_tree.png", sep=""), width=no.subsamples*1000, height=1000)
    plotTree(tree, font.size=5)
    dev.off()
  }  
  
  # Write the index
  write.table(dataset_index,file=paste(outpath,datasetname,".txt",sep=""),quote=F,row.names=F, sep="\t")
}

mutationCopyNumberToMutationBurden<-function(copyNumber,totalCopyNumber,cellularity,normalCopyNumber = rep(2,length(copyNumber))){
  burden = copyNumber*cellularity/(cellularity*totalCopyNumber+normalCopyNumber*(1-cellularity))
  burden[is.nan(burden)|(burden<0.000001)]=0.000001
  burden[burden>0.999999]=0.999999
  return(burden)  
}

mutationBurdenToMutationCopyNumber<-function(burden,totalCopyNumber,cellularity,normalCopyNumber = rep(2,length(burden))){
  mutCopyNumber = burden/cellularity*(cellularity*totalCopyNumber+normalCopyNumber*(1-cellularity))
  mutCopyNumber[is.nan(mutCopyNumber)]=0
  return(mutCopyNumber)
}


# infile = "/lustre/scratch110/sanger/sd11/dirichlet/simulated/Data/simulated_a_0010.tree.txt"
# tree = read.table(infile, header=T)
# row.names(tree) = tree$label
# png("simulated_a_0010_expected_tree.png", width=sum(grepl("theta", colnames(tree)))*1000, height=1000)
# plotTree(tree, font.size=5)
# dev.off()
