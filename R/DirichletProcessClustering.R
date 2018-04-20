
#' @param mutCount
#' @param WTCount
#' @param totalCopyNumber
#' @param copyNumberAdjustment
#' @param mutation.copy.number
#' @param cellularity
#' @param output_folder
#' @param no.iters
#' @param no.iters.burn.in
#' @param subsamplesrun
#' @param samplename
#' @param conc_param
#' @param cluster_conc
#' @param mut.assignment.type
#' @param most.similar.mut
#' @param mutationTypes
#' @param min.frac.snvs.cluster
#' @param max.considered.clusters
DirichletProcessClustering <- function(mutCount, WTCount, totalCopyNumber, copyNumberAdjustment, mutation.copy.number, cellularity, output_folder, no.iters, no.iters.burn.in, subsamplesrun, samplename, conc_param, cluster_conc, mut.assignment.type, most.similar.mut, mutationTypes, min.frac.snvs.cluster, max.considered.clusters) {
  #
  # Run the regular Dirichlet Process based method. Will perform clustering using the given data. The method
  # decides automatically whether the 1D or nD method is run based on the number of samples given at the input.
  # The number of samples is determined through the number of columns of the input.
  #
  # The nD method will yield a series of figures in which each sample is plotted against each other sample. 
  # The 1D method yields a density plot for just the single sample.
  #
  print(dim(mutCount))
  if(!file.exists(output_folder)){
    dir.create(output_folder)
  }
  GS.data = subclone.dirichlet.gibbs(mutCount=mutCount,
                                    WTCount=WTCount,
                                    totalCopyNumber=totalCopyNumber,
                                    copyNumberAdjustment=copyNumberAdjustment,
                                    cellularity=cellularity,
                                    iter=no.iters,
                                    conc_param=conc_param,
                                    cluster_conc=cluster_conc,
                                    C=max.considered.clusters)
  
  save(file=paste(output_folder, "/", samplename, "_gsdata.RData", sep=""), GS.data)
  
  # nD dataset, plot sample versus sample
  if (ncol(mutCount) > 1) {
    ########################
    # Plot density and Assign mutations to clusters - nD
    ########################
    for (i in 1:(length(subsamplesrun)-1)) {
      for (j in (i+1):length(subsamplesrun)) {
        imageFile = paste(output_folder,"/",samplename,subsamplesrun[i],subsamplesrun[j],"_iters",no.iters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_2D_binomial.png",sep="")
        density = Gibbs.subclone.density.est(mutation.copy.number[,c(i,j)]/copyNumberAdjustment[,c(i,j)],
                                             GS.data,
                                             imageFile, 
                                             post.burn.in.start = no.iters.burn.in, 
                                             post.burn.in.stop = no.iters, 
                                             samplenames = paste(samplename,subsamplesrun[c(i,j)],sep=""),
                                             indices=c(i,j))  
        save(file=paste(output_folder,"/",samplename, subsamplesrun[i], subsamplesrun[j], "_densityoutput.RData", sep=""), GS.data, density)
      }
    }
    
    # Assign mutations to clusters using one of the different assignment methods
    opts = list(samplename=samplename, subsamplenames=subsamplesrun, no.iters=no.iters, no.iters.burn.in=no.iters.burn.in, no.iters.post.burn.in=no.iters-no.iters.burn.in, outdir=output_folder)
    print("Assigning mutations to clusters")
    if (mut.assignment.type == 1) {
      consClustering = multiDimensionalClustering(mutation.copy.number=mutation.copy.number, 
                                                  copyNumberAdjustment=copyNumberAdjustment, 
                                                  GS.data=GS.data, 
                                                  density.smooth=0.01, 
                                                  opts=opts)
    } else if (mut.assignment.type == 2) {
      consClustering = mutation_assignment_em(GS.data=GS.data,
                                              mutCount=mutCount, 
                                              WTCount=WTCount, 
                                              subclonal.fraction=mutation.copy.number/copyNumberAdjustment, 
                                              node.assignments=GS.data$S.i, 
                                              opts=opts)
      
    } else if (mut.assignment.type == 3) {  
      warning("binom mut assignment not implemented for multiple timepoints")
      q(save="no", status=1)
    
    } else if (mut.assignment.type == 4) {  
      warning("MPEAR mut assignment not implemented for multiple timepoints")
      q(save="no", status=1)
      
    } else {
      warning(paste("Unknown mutation assignment type", mut.assignment.type, sep=" "))
      q(save="no", status=1)
    }
    return(consClustering)
    
    # 1D dataset, plot just the single density
  } else {
    ########################
    # Plot density and Assign mutations to clusters - 1D
    ########################
    # 1D dataset, plot just the single density
    wd = getwd()
    setwd(output_folder)
    res = Gibbs.subclone.density.est.1d(GS.data, 
                                            paste(samplename,"_DirichletProcessplot.png", sep=''), 
                                            samplename=samplename,
                                            post.burn.in.start=no.iters.burn.in, 
                                            post.burn.in.stop=no.iters,
                                            y.max=15,
					                                  x.max=NA, 
                                            mutationCopyNumber=mutation.copy.number, 
                                            no.chrs.bearing.mut=copyNumberAdjustment)
    density = res$density
    polygon.data = res$polygon.data
    
    # Assign mutations to clusters using one of the different assignment methods
    opts = list(samplename=samplename, subsamplenames=subsamplesrun, no.iters=no.iters, no.iters.burn.in=no.iters.burn.in, no.iters.post.burn.in=no.iters-no.iters.burn.in, outdir=output_folder)

    print("Assigning mutations to clusters")
    if (mut.assignment.type == 1) {
      subclonal.fraction = mutation.copy.number / copyNumberAdjustment
      subclonal.fraction[is.nan(subclonal.fraction)] = 0
      consClustering = oneDimensionalClustering(samplename, subclonal.fraction, GS.data, density, no.iters, no.iters.burn.in)
      setwd(wd) 
    } else if (mut.assignment.type == 2) {
      setwd(wd) # set the wd back earlier. The oneD clustering and gibbs sampler do not play nice yet and need the switch, the em assignment doesnt
      consClustering = mutation_assignment_em(GS.data=GS.data,
                                              mutCount=mutCount, 
                                              WTCount=WTCount, 
                                              subclonal.fraction=mutation.copy.number/copyNumberAdjustment, 
                                              node.assignments=GS.data$S.i, 
                                              opts=opts)
      
    } else if (mut.assignment.type == 3) {  
      consClustering = mutation_assignment_binom(clustering_density=density,
                                                 mutCount=mutCount, 
                                                 WTCount=WTCount, 
                                                 copyNumberAdjustment=copyNumberAdjustment, 
                                                 tumourCopyNumber=totalCopyNumber,
                                                 normalCopyNumber=array(2, dim(mutCount)),
                                                 cellularity=cellularity)
      
    } else if (mut.assignment.type == 4) {  
      consClustering = mutation_assignment_mpear(GS.data=GS.data, 
                                                 no.iters=no.iters, 
                                                 no.iters.burn.in=no.iters.burn.in, 
                                                 min.frac.snvs.cluster=min.frac.snvs.cluster, 
                                                 dataset=dataset, 
                                                 samplename=samplename, 
                                                 outdir=output_folder,
                                                 density=density)
      
    } else {
      warning(paste("Unknown mutation assignment type", mut.assignment.type, sep=" "))
      q(save="no", status=1)
    }

    print(head(consClustering$best.node.assignments))
    print("")
    print(consClustering$cluster.locations)
    print("")
    print(head(consClustering$all.likelihoods))
    print("")
    print(head(consClustering$best.assignment.likelihoods))
    
    # Make a second set of figures with the mutation assignments showing
    # Replot the data with cluster locations
    plot1D(density=density, 
           polygon.data=polygon.data, 
           pngFile=paste(output_folder, "/", samplename, "_DirichletProcessplot_with_cluster_locations.png", sep=""), 
           density.from=0, 
           x.max=1.5, 
           mutationCopyNumber=mutation.copy.number, 
           no.chrs.bearing.mut=copyNumberAdjustment,
           samplename=samplename,
           cluster.locations=consClustering$cluster.locations,
           mutation.assignments=consClustering$best.node.assignments)
    
    plot1D_2(density=density, 
             polygon.data=polygon.data, 
             pngFile=paste(output_folder, "/", samplename, "_DirichletProcessplot_with_cluster_locations_2.png", sep=""), 
             density.from=0, 
             x.max=1.5, 
             mutationCopyNumber=mutation.copy.number, 
             no.chrs.bearing.mut=copyNumberAdjustment,
             samplename=samplename,
             cluster.locations=consClustering$cluster.locations,
             mutation.assignments=consClustering$best.node.assignments,
             mutationTypes=mutationTypes)
    
    setwd(wd)
    return(consClustering)
  }
  
}

write.strengths.table = function(dat, removed_indices, filename) {
  #
  # Adds in an empty column/row for each of the removed_indices and writes
  # the subsequent matrix to disk
  #
  write.table(add.muts.back.in(dat, removed_indices),filename,sep="\t",row.names=F,quote=F,col.names=F)
}

add.muts.back.in = function(dat, removed_indices, def.value=0) {
  #
  # Adds in empty columns and rows for mutations that were removed.
  # Mutations are added sequentially, so this method expects indices
  # of removed mutations in the original (full) matrix. The empty
  # mutations will be added in the place where they were removed,
  # keeping the order in tact.
  #
  for (i in removed_indices) { 
    if (i==1) {
      dat = cbind(rep(def.value, nrow(dat)), dat)
      dat = rbind(rep(def.value, ncol(dat)), dat)
    } else if (i >= ncol(dat)) {
      	dat = cbind(dat, rep(def.value, nrow(dat)))
      	dat = rbind(dat, rep(def.value, ncol(dat)))
    } else {
        dat = cbind(dat[,1:(i-1)], rep(def.value, nrow(dat)), dat[,i:ncol(dat)])
        dat = rbind(dat[1:(i-1),], rep(def.value, ncol(dat)), dat[i:nrow(dat),])
    }
  }
  return(dat)
}

