#run dp
RunDP <- function(analysis_type, datpath, cndatafiles = NA, co_cluster_cna = F, sample.snvs.only = T, remove.snvs = F ,generate_cluster_ordering = F,
                  run_params, sample_params, advanced_params, outdir, cna_params=NULL, mutphasingfiles=NULL) { 
  
  #####################################################################################
  # Unpack parameters
  #####################################################################################
  attach(run_params)
  attach(sample_params)
  attach(advanced_params)
  if (!is.null(cna_params)) {
    attach(cna_params)
  }
  
  #####################################################################################
  # Check input
  #####################################################################################
  # Check whether a supported analysis_type was supplied
  ##################
  supported_commands = c('nd_dp', 'replot_1d', 'replot_nd', 'reassign_muts_1d', 'reassign_muts_nd',"phasing_ass") 
  ###################
  if (!(analysis_type %in% supported_commands)) {
    print(paste("Type of analysis", analysis_type, "unknown."))
    print(paste(c("Specify either ", supported_commands)), sep=" ")
    q(save="no", status=1)
  }
  # Check whether the mut.assignment.type is supported
  supported_mut.assignment.methods = c(1,2,3,4)
  if (!(mut.assignment.type %in% supported_mut.assignment.methods)) {
    print(paste("Type of mutation assignment method", mut.assignment.type, "unknown."))
    print(paste(c("Specify either ", supported_mut.assignment.methods)), sep=" ")
    q(save="no", status=1)
  }
  
  #####################################################################################
  # Setup
  #####################################################################################
  set.seed(seed)
  
  # Create the output directory
  if (!file.exists(outdir)) { dir.create(outdir) }
  
  # Path for output files
  outfiles.prefix = file.path(outdir, paste(samplename, "_", no.iters, "iters_", no.iters.burn.in, "burnin", sep=""))
  
  #####################################################################################
  # Loading data
  #####################################################################################
  print("Loading data...")
  # Load data from disk if there was already a dataset object, otherwise create a new one
  if (file.exists(paste(outdir, "/dataset.RData", sep=""))) {
    
    load(file.path(outdir, "dataset.RData"))
    cndata = dataset$cndata
    mutphasing = dataset$mutphasing
    
    if (!is.null(cndata)) {
      cndata_params = list()
      cndata_params$cndata = cndata
      cndata_params$add.conflicts = add.conflicts
      cndata_params$cna.conflicting.events.only = cna.conflicting.events.only
      cndata_params$num.clonal.events.to.add = num.clonal.events.to.add
      cndata_params$min.cna.size = min.cna.size
    } else {
      cndata_params = NULL
    }
    
    
  } else {
    list_of_datafiles = paste(datpath, datafiles, sep="/")
    cndatafiles = paste(datpath, cndatafiles, sep="")  
    # Note that the phase column is not used
    dataset = load.data(list_of_datafiles, 
                        cellularity=cellularity, 
                        Chromosome="chr", 
                        position="end",
                        WT.count="WT.count", 
                        mut.count="mut.count", 
                        subclonal.CN="subclonal.CN", 
                        no.chrs.bearing.mut="no.chrs.bearing.mut", 
                        mutation.copy.number="mutation.copy.number", 
                        subclonal.fraction="subclonal.fraction", 
                        phase=NULL, # disabled for now, as the data is not used
                        is.male=is.male,
                        is.vcf=F, # reading of VCF input files os disabled
                        ref.genome.version="hg19", # reading of VCF input files is disabled, this parameter is not used
                        supported_chroms=supported_chroms)
    
    if (co_cluster_cna & !is.na(cndatafiles)) {
      cndata = load.cn.data(cndatafiles)
      cndata_params = list()
      cndata_params$cndata = cndata
      cndata_params$add.conflicts = add.conflicts
      cndata_params$cna.conflicting.events.only = cna.conflicting.events.only
      cndata_params$num.clonal.events.to.add = num.clonal.events.to.add
      cndata_params$min.cna.size = min.cna.size
    } else {
      cndata = NULL
      cndata_params = NULL
    }
    
    if (!is.null(mutphasingfiles)) {
      print("Loading mutation phasing info")
      mutphasing = NULL
      #for (infile in mutphasingfiles) {
       # load(infile)
      #}
    } else {
      mutphasing = NULL
    }
    
    save(file=file.path(outdir, "dataset.RData"), dataset)
  }
  
  
  #####################################################################################
  # Loading CN data and deal with the aftermath of that
  #####################################################################################
  # Unpack the copy number inclusion parameters
  if (!is.null(cndata_params)) {
    cndata = cndata_params$cndata
    add.conflicts = cndata_params$add.conflicts
    cna.conflicting.events.only = cndata_params$cna.conflicting.events.only
    num.clonal.events.to.add = cndata_params$num.clonal.events.to.add
    min.cna.size = cndata_params$min.cna.size
  }
  
  # Check if co-clustering of copy number data is in order
  resave.dataset = F # A boolean that keeps track of whether the dataset should be saved again. Set this to TRUE if the dataset changes.
  if (!is.null(dataset$cndata)) {
    # In case of a rerun, pull out the cndata
    cndata = dataset$cndata
  } else if (!is.null(cndata)) {
    dataset = add.in.cn.as.snv.cluster(dataset, 
                                       cndata, 
                                       cellularity=cellularity,
                                       add.conflicts=add.conflicts, 
                                       conflicting.events.only=cna.conflicting.events.only, 
                                       num.clonal.events.to.add=num.clonal.events.to.add,
                                       min.cna.size=min.cna.size)
    resave.dataset = T
  }
  
  #####################################################################################
  # Add in mutation phasing
  #####################################################################################
  # Check for mutationphasing info
  if (!is.null(dataset$mutphasing)) {
    mutphasing = dataset$mutphasing
  } else if (!is.null(mutphasing)) {
    dataset = add.mutphasing(dataset, mutphasing, add.conflicts=add.conflicts)
    resave.dataset = T
  }
  
  # Perform sampling
  if (!is.na(num_muts_sample) & num_muts_sample!="NA") {
    if (is.na(dataset$full.data)) {
      dataset = sample_mutations(dataset, num_muts_sample, sample.snvs.only=sample.snvs.only, remove.snvs=remove.snvs)
      most.similar.mut = dataset$most.similar.mut
      resave.dataset = T
    }       
    most.similar.mut = dataset$most.similar.mut
  } else {
    most.similar.mut = NA
  }
  dataset$cndata = cndata
  # The dataset object was modified, so save it
  if (resave.dataset) { save(file=file.path(outdir, "dataset.RData"), dataset) }
  
  if (analysis_type == 'nd_dp') {
    print("Running DPClust...")
    ##############################
    # nD DP clustering
    ##############################
    clustering = DirichletProcessClustering(mutCount=dataset$mutCount, 
                                            WTCount=dataset$WTCount, 
                                            no.iters=no.iters, 
                                            no.iters.burn.in=no.iters.burn.in, 
                                            cellularity=cellularity, 
                                            totalCopyNumber=dataset$totalCopyNumber, 
                                            mutation.copy.number=dataset$mutation.copy.number,
                                            copyNumberAdjustment=dataset$copyNumberAdjustment, 
                                            mutationTypes=dataset$mutationType,
                                            samplename=samplename, 
                                            subsamplesrun=subsamples,
                                            output_folder=outdir, 
                                            conc_param=conc_param, 
                                            cluster_conc=cluster_conc,
                                            mut.assignment.type=mut.assignment.type,
                                            most.similar.mut=most.similar.mut,
                                            max.considered.clusters=max.considered.clusters)
    
  } else if (analysis_type == "replot_1d") {
    print("Running Remaking plots...")
    ##############################
    # Replot 1D clustering
    ##############################
    density = read.table(file.path(outdir, paste(samplename, "_DirichletProcessplotdensity.txt", sep="")), header=T)
    polygon.data = read.table(file.path(outdir, paste(samplename, "_DirichletProcessplotpolygonData.txt", sep="")), header=T)
    load(file=paste(outfiles.prefix, "_bestConsensusResults.RData", sep=""))
    replot_1D(outdir=outdir,
              outfiles.prefix=outfiles.prefix,
              samplename=samplename, 
              dataset=dataset, 
              clustering=clustering, 
              density=density, 
              polygon.data=polygon.data)
    
  } else if (analysis_type == "replot_nd") {  
    print("Remaking plots...")
    ##############################
    # Replot nD clustering
    ##############################
    load(file=file.path(outdir, paste(samplename, "_gsdata.RData", sep="")))
    replot_nD(outdir=outdir, 
              outfiles.prefix=outfiles.prefix, 
              samplename=samplename, 
              subsamples=subsamples, 
              dataset=dataset,
              no.iters=no.iters,
              clustering=clustering, 
              conc_param=conc_param, 
              cluster_conc=cluster_conc)
    
  } else if (analysis_type == "reassign_muts_1d") {
    print("Reassigning mutations...")
    ##############################
    # Reassign mutations to clusters using a previous 1D clustering run
    ##############################
    load(file=file.path(outdir, paste(samplename, "_gsdata.RData", sep="")))
    res = reassign_1D(outdir=outdir, 
                      samplename=samplename, 
                      no.iters=no.iters, 
                      no.iters.burn.in=no.iters.burn.in, 
                      dataset=dataset, 
                      cellularity=cellularity, 
                      GS.data=GS.data,
                      conc_param=conc_param, 
                      cluster_conc=cluster_conc,
                      mut.assignment.type=mut.assignment.type)
    clustering = res$clustering
    outfiles.prefix = res$outfiles.prefix
    
  } else if (analysis_type == "reassign_muts_nd") {
    print("Reassigning mutations...")
    ##############################
    # Reassign mutations to clusters using a previous nD clustering run
    ##############################
    load(file=paste(outfiles.prefix, "_bestConsensusResults.RData", sep=""))
    res = reassign_nd(outdir=outdir, 
                      samplename=samplename, 
                      subsamples=subsamples, 
                      no.iters=no.iters, 
                      no.iters.burn.in=no.iters.burn.in, 
                      dataset=dataset, 
                      GS.data=GS.data,
                      conc_param=conc_param, 
                      cluster_conc=cluster_conc,
                      mut.assignment.type=mut.assignment.type)
    
    clustering = res$clustering
    outfiles.prefix = res$outfiles.prefix
    
  } else if (analysis_type == "phasing_ass") {
    print("Running DPClust with Phasing information")
    ##############################
    ##############################
    # Reassign mutations to clusters using a previous nD clustering run
    ##############################
    
    num.timepoints = dim(dataset$WTCount)[2]
    muts.loc.matrix = data.frame(Chr=dataset$chromosome[,1],
                                 Pos=dataset$position[,1])

    mutphasingfiles = paste(datpath, mutphasingfiles, sep="/")
    attach(dataset)
    if(num.timepoints>1){
      ###############
      
      loc.info = list()
      phasing.info = list()
      phasing.index = c()
      for(i.sample in 1:num.timepoints){  #load phasing file
        clphasing.file.temp = mutphasingfiles[i.sample] #phasing file 
        load( clphasing.file.temp ) #
        cl.phasing.pass = data.frame(cl.phasing.pass)
        
        rm.index = c()
        for(i.cl in 1:dim(cl.phasing.pass)[1]){
          mut1.temp = which(muts.loc.matrix$Chr==cl.phasing.pass[i.cl,1] & muts.loc.matrix$Pos==cl.phasing.pass[i.cl,2])
          mut2.temp = which(muts.loc.matrix$Chr==cl.phasing.pass[i.cl,1] & muts.loc.matrix$Pos==cl.phasing.pass[i.cl,6])
          if(length(mut1.temp)&length(mut2.temp)){
            cl.phasing.pass[i.cl,"Mut1_index"] = mut1.temp
            cl.phasing.pass[i.cl,"Mut2_index"] = mut2.temp
          }else{
            rm.index = c(rm.index,i.cl)
          }
        }
        if(length(rm.index)){
          cl.phasing.pass = cl.phasing.pass[-rm.index,]
        }  #remove pairs that are not included in dataset
        
        phasing.info = c(phasing.info, list(cl.phasing.pass))
        
        phasing.index.temp = paste0(cl.phasing.pass$Mut1_index,"-",cl.phasing.pass$Mut2_index)
        phasing.index = c(phasing.index,list(phasing.index.temp))
      }
      
      
      phasing.index.unique = unique(unlist(phasing.index))
      phasing.pair.num = length(phasing.index.unique)
      all.phasing.reads = array(0,c(phasing.pair.num,4,num.timepoints))#Fistr 4: Read on both sites; Last 2: index
      all.phasing.index = all.phasing.info = matrix(NA,nrow=phasing.pair.num,ncol=2)
      for(i.sample in 1:num.timepoints){
        phasing.loc = sapply(phasing.index[[i.sample]],function(t){which(phasing.index.unique==t)})
        all.phasing.reads[phasing.loc,1:4,i.sample] = as.matrix(phasing.info[[i.sample]][,10:13])
        all.phasing.index[phasing.loc,1] = phasing.info[[i.sample]][,16]
        all.phasing.index[phasing.loc,2] = phasing.info[[i.sample]][,17]
        all.phasing.info[phasing.loc,1] = phasing.info[[i.sample]][,14]
        all.phasing.info[phasing.loc,2] = phasing.info[[i.sample]][,15]
      }
      
      mut_single_1 = mutCount[all.phasing.index[,1],] - all.phasing.reads[,2,] - all.phasing.reads[,3,]
      
      mut_single_2 = mutCount[all.phasing.index[,2],] - all.phasing.reads[,2,] - all.phasing.reads[,4,]
      
      
      wt_single_1 = WTCount[all.phasing.index[,1],] -  all.phasing.reads[,1,] - all.phasing.reads[,4,]
      wt_single_2 = WTCount[all.phasing.index[,2],] -  all.phasing.reads[,2,] - all.phasing.reads[,3,]
      
      nonneg.indi = apply(cbind(mut_single_1,mut_single_2,wt_single_1,wt_single_2),1,function(t){length(which(t<0))})
      nonneg.index = which(nonneg.indi == 0)
      if(length(nonneg.index>0)){
        all.phasing.index = all.phasing.index[nonneg.index,]
        all.phasing.reads = all.phasing.reads[nonneg.index,,]
      }else{
        print("No Phasing information available")
        q(save="no", status=1)
      }
      phasing.pair.num = dim(all.phasing.index)[1]
      
      mutCount[all.phasing.index[,1],] =  mut_single_1[nonneg.index,]
      mutCount[all.phasing.index[,2],] =  mut_single_2[nonneg.index,]
      WTCount[all.phasing.index[,1],] = wt_single_1[nonneg.index,]
      WTCount[all.phasing.index[,2],] = wt_single_2[nonneg.index,]
      
    }else{
      load(mutphasingfiles) #
      
      rm.index = c()  #phasing pairs that are not listed in dataset
      for(i.cl in 1:dim(cl.phasing.pass)[1]){
        mut1.temp = which(muts.loc.matrix$Chr==cl.phasing.pass[i.cl,1] & muts.loc.matrix$Pos==cl.phasing.pass[i.cl,2])
        mut2.temp = which(muts.loc.matrix$Chr==cl.phasing.pass[i.cl,1] & muts.loc.matrix$Pos==cl.phasing.pass[i.cl,6])
        
        if(length(mut1.temp)&length(mut2.temp)){
          cl.phasing.pass[i.cl,"Mut1_index"] = mut1.temp
          cl.phasing.pass[i.cl,"Mut2_index"] = mut2.temp
        }else{
          rm.index = c(rm.index,i.cl)
        }
      }
      
      
      if(length(rm.index)){
        cl.phasing.pass = cl.phasing.pass[-rm.index,]
      }
      
      cl.phasing.pass = data.frame(cl.phasing.pass)
      mutCount[cl.phasing.pass$Mut1_index] = cl.phasing.pass$mut.count1 #need to make sure the type is numeric
      mutCount[cl.phasing.pass$Mut2_index] = cl.phasing.pass$mut.count2
      WTCount[cl.phasing.pass$Mut1_index] = cl.phasing.pass$WT.count1
      WTCount[cl.phasing.pass$Mut2_index] = cl.phasing.pass$WT.count2
      
      all.phasing.index = cl.phasing.pass[,c("Mut1_index","Mut2_index")]
      all.phasing.reads = array(NA,dim=c(dim(cl.phasing.pass)[1],4,num.timepoints))
      all.phasing.reads[,,1] = unlist(cl.phasing.pass[,c("Num_WT_WT","Num_MUT_MUT","Num_WT_MUT","Num_MUT_WT")])
      
      all.phasing.info = cl.phasing.pass[,c("Copy_Major","Copy_Minor")]
      
      
    }
    
    
    colnames(all.phasing.index) = c("Mut1_index","Mut2_index")
    colnames(all.phasing.info) = c("Copy_Major","Copy_Minor")
    all.phasing.index = as.data.frame(all.phasing.index)
    all.phasing.info = as.data.frame(all.phasing.info)
    
    
    assign_result =   PhasingAssignment(mutCount=mutCount, 
                                        WTCount=WTCount, 
                                        no.iters=no.iters, #input
                                        num.timepoints = num.timepoints,
                                        error.rate = error.rate,
                                        killclu = kill.clu.frac,
                                        cellularity=cellularity, 
                                        conc_param = 1,
                                        C = num_clu_samplers,
                                        all.phasing.index = all.phasing.index ,
                                        all.phasing.reads = all.phasing.reads,
                                        all.phasing.info = all.phasing.info ,
                                        totalCopyNumber=dataset$totalCopyNumber, 
                                        mutation.copy.number=dataset$mutation.copy.number,
                                        copyNumberAdjustment=dataset$copyNumberAdjustment, 
                                        mutationTypes=dataset$mutationType,
                                        samplename=samplename, 
                                        subsamplesrun=subsamples,
                                        output_folder=outdir, 
                                        cluster_conc=cluster_conc,
                                        mut.assignment.type=mut.assignment.type,
                                        most.similar.mut=most.similar.mut,
                                        max.considered.clusters=max.considered.clusters,
                                        normalCopyNumber=array(2,dim(mutCount)))
    
    
    
    GS.data = list(S.i=assign_result$S.i, V.h=assign_result$V.h, pi.h=assign_result$pi.h,mutBurdens=assign_result$mutBurdens, 
                   alpha=assign_result$alpha, y1=assign_result$mutCount, N1=assign_result$N1,
                   whole.tree = assign_result$whole.tree,
                   tree.selected = assign_result$tree.selected)
    
    save(file=file.path(outdir, paste(samplename, "_gsdata.RData", sep="")), GS.data)
    
    tree.output = assign_result$tree.selected
    whole.tree = assign_result$whole.tree
    
    #find the most common tree after burn-in
    
    tree.size.iter = sapply(2:no.iters,function(t){dim(whole.tree[[tree.output[t]]])[1]+1})
    tree.size.iter = c(NA,tree.size.iter)
    tree.size.final = tree.size.iter[no.iters]
    converged.iter = which(tree.size.iter==tree.size.final)
    used.iters = max(min(converged.iter),no.iters.burn.in+1)
    most.fre.tree.index = sort(table(tree.output[used.iters:no.iters]),decreasing=TRUE)[1]
    most.fre.tree.index = as.numeric(names(most.fre.tree.index))
    most.fre.tree = whole.tree[[most.fre.tree.index]]
    
    write.table(paste0(outdir,samplename,"_most_frequent_tree.txt"),sep="\t",row.names=F,quote=F)
    
    #print
    # nD dataset, plot sample versus sample
    if (ncol(mutCount) > 1) {
      ########################
      # Plot density and Assign mutations to clusters - nD
      ########################
      print("Estimating density between pairs of samples...")
      subsamplesrun = subsamples
      for (i in 1:(length(subsamplesrun)-1)) {
        for (j in (i+1):length(subsamplesrun)) {
          print(paste("Samples", subsamplesrun[i], "and", subsamplesrun[j], sep=" "))
          imageFile = file.path(outdir, paste(samplename,subsamplesrun[i],subsamplesrun[j],"_iters",no.iters,"_concParam",conc_param,"_clusterWidth",1/cluster_conc,"_2D_binomial.png",sep=""))
          density = Gibbs.subclone.density.est(mutation.copy.number[,c(i,j)]/copyNumberAdjustment[,c(i,j)],
                                               GS.data,
                                               imageFile,
                                               post.burn.in.start = no.iters.burn.in,
                                               post.burn.in.stop = no.iters,
                                               samplenames = paste(samplename,subsamplesrun[c(i,j)],sep=""),
                                               indices=c(i,j))
          save(file=file.path(outdir, paste(samplename, subsamplesrun[i], subsamplesrun[j], "_densityoutput.RData", sep="")), GS.data, density)
        }
      }
      
    } else {
      ########################
      # Plot density and Assign mutations to clusters - 1D
      ########################
      # 1D dataset, plot just the single density
      wd = getwd()
      setwd(outdir)
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
      opts = list(samplename=samplename, subsamplenames=subsamplesrun, no.iters=no.iters, no.iters.burn.in=no.iters.burn.in, no.iters.post.burn.in=no.iters-no.iters.burn.in, outdir)
      
      
      # Make a second set of figures with the mutation assignments showing
      # Replot the data with cluster locations
      plot1D(density=density, 
             polygon.data=polygon.data, 
             pngFile=paste(outdir, "/", samplename, "_DirichletProcessplot_with_cluster_locations.png", sep=""), 
             density.from=0, 
             x.max=1.5, 
             mutationCopyNumber=mutation.copy.number, 
             no.chrs.bearing.mut=copyNumberAdjustment,
             samplename=samplename,
             cluster.locations=consClustering$cluster.locations,
             mutation.assignments=consClustering$best.node.assignments)
      
      plot1D_2(density=density, 
               polygon.data=polygon.data, 
               pngFile=paste(outdir, "/", samplename, "_DirichletProcessplot_with_cluster_locations_2.png", sep=""), 
               density.from=0, 
               x.max=1.5, 
               mutationCopyNumber=mutation.copy.number, 
               no.chrs.bearing.mut=copyNumberAdjustment,
               samplename=samplename,
               cluster.locations=consClustering$cluster.locations,
               mutation.assignments=consClustering$best.node.assignments,
               mutationTypes=mutationTypes)
      
    }
    
  } else {
    print(paste("Unknown type of analysis",analysis_type))
    q(save="no", status=1)
  }
  
  ####################################################################################################################
  # Write the final output
  ####################################################################################################################
  if (all(!analysis_type %in% c('replot_1d', 'replot_nd'))) {
    print("Writing out final output...")
    
    # Load the MCMC output as its needed to get cluster confidence intervals
    density_file = file.path(outdir, paste(samplename, "_DirichletProcessplotdensity.txt", sep=""))
    if (file.exists(density_file)) {
      density = read.table(density_file, header=T)
      colnames(density)[1] = "fraction.of.tumour.cells"
    } else {
      density = NA
    }
    
    polygon_file = file.path(outdir, paste(samplename, "_DirichletProcessplotpolygonData.txt", sep=""))
    if (file.exists(polygon_file)) {
      polygon.data = read.table(polygon_file, header=T)
    } else {
      polygon.data = NA
    }
    
    load(file.path(outdir, paste(samplename, "_gsdata.RData", sep="")))
    write_tree = analysis_type != 'nd_dp' & analysis_type != 'reassign_muts_1d' & analysis_type != 'reassign_muts_nd'
    writeStandardFinalOutput(clustering=clustering, 
                             dataset=dataset,
                             most.similar.mut=most.similar.mut,
                             outfiles.prefix=outfiles.prefix,
                             outdir=outdir,
                             samplename=samplename,
                             subsamplenames=subsamples,
                             GS.data=GS.data,
                             density=density,
                             polygon.data=polygon.data,
                             no.iters=no.iters, 
                             no.iters.burn.in=no.iters.burn.in,
                             assign_sampled_muts=assign_sampled_muts,
                             write_tree=write_tree,
                             generate_cluster_ordering=generate_cluster_ordering,
                             min_muts_cluster=min_muts_cluster,
                             min_frac_muts_cluster=min_frac_muts_cluster)
  }
  
  ####################################################################################################################
  # Remove intermediate files
  ####################################################################################################################
  if (!keep_temp_files) {
    .remove_file = function(filename) {
      if (file.exists(filename)) file.remove(filename)
    }
    
    # 1D method and general files
    .remove_file(paste(outfiles.prefix, "_removedMutationsIndex.txt", sep=""))
    .remove_file(file.path(outdir, paste(samplename, "_DP_and_cluster_info.txt", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_DirichletProcessplot.png", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_DirichletProcessplot_with_cluster_locations.png", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_DirichletProcessplotdensity.txt", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_DirichletProcessplotpolygonData.txt", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_localOptima.txt", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_optimaInfo.txt", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_gsdata.RData", sep="")))
    .remove_file(file.path(outdir, "dataset.RData"))
    
    # nD method files
    .remove_file(file.path(outdir, paste(samplename, "_DP_and cluster_info_0.01.txt", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_confInts_0.01.txt", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_localHighConfidenceMultidimensionalOptima_0.01.txt", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_localMultidimensionalOptima_0.01.txt", sep="")))
    .remove_file(file.path(outdir, paste(samplename, "_optimaInfo_0.01.txt", sep="")))
    
    for (i in 1:(length(subsamples)-1)) {
      for (j in (i+1):length(subsamples)) {
        .remove_file(file.path(outdir, paste(samplename, subsamples[i], subsamples[j], "_densityoutput.RData", sep="")))
        .remove_file(file.path(outdir, paste(samplename, subsamples[i], subsamples[j], "_densityoutput.csv", sep="")))
        .remove_file(file.path(outdir, pattern=glob2rx(paste(samplename, subsamples[i], subsamples[j], "*densityData1.csv", sep="")), full.names=T))
        density_csv_files = list.files(outdir, pattern=glob2rx(paste(samplename, subsamples[i], subsamples[j], "*vals.csv", sep="")), full.names=T)
        for (infile in density_csv_files) { .remove_file(infile) }
      }
    }
    
    nd_density_files = list.files(outdir, pattern="_2D_binomial_")
    if (length(nd_density_files) > 0) { file.remove(nd_density_files) }
  }
  print("Done.")
}
