#
# Plotting functions
#

#' Original plotting function: Creates the 1D density plot in either allele fraction, mutation copy number or fraction of tumour cells space.
#' 
#' @param density Is a two column data frame. First column x-axis density line, second column y-axis density line.
#' @param polygon.data Contains the coordinates of the 95\% confidence interval around the density line.
#' @param pngFile Where to save the figure (Default: NA)
#' @param density.from Startpoint from where density is drawn (Default: 0)
#' @param x.max Maximum values on x axis in the plot (Default: NA, determine from data)
#' @param y.max Maximum values on y axis in the plot (Default: NA, determine from data)
#' @param y MutCount data, only used when mutationCopyNumber is not supplied (Default: NULL)
#' @param N Total count (mutCount+WTCount), only used when mutationCopyNumber is not supplied (Default: NULL)
#' @param mutationCopyNumber Mutation copy number data, when not supplied y and N are required (Default: NULL)
#' @param no.chrs.bearing.mut Copynumber (multiplicity) adjustment per mutation, when not supplied the plot will either contain mutation copy number or allele frequencies (Default: NULL)
#' @param samplename Name of the sample under analysis (Default: empty string)
#' @param CARL Legacy annotation parameter, no longer used (Default: empty numeric)
#' @param cluster.locations Locations in the space to be plotted where clusters reside (Default: NULL)
#' @param mutation.assignments Mutation assignments (Default: NULL)
#' @param mutationTypes Type of each mutation, different mutation types can be plot with different colours (Default: NULL)
#' @author dw9, sd11
plot1D = function(density, polygon.data, pngFile=NA, density.from=0, x.max=NA, y.max=NA, y=NULL, N=NULL, mutationCopyNumber=NULL, no.chrs.bearing.mut=NULL,samplename="",CALR=numeric(0), cluster.locations=NULL, mutation.assignments=NULL, mutationTypes=NULL) {
  if (!is.na(pngFile)) { png(filename=pngFile,width=1500,height=1000) }
  
  # Convert data into the space that was used for the clustering. This is done dynamically through the given data.
  xlabel = "Mutation Copy Number"
  if(is.null(mutationCopyNumber)){
    print("No mutationCopyNumber. Using mutation burden")
    if (is.null(y) | is.null(N)) { 
      print("When not supplying mutationCopyNumber, y (mutCount) and N (totalCount) are required")
      q(save="no", status=1)
    }
    mutationCopyNumber = y/N
    xlabel = "Mutation Burden"
  }
  
  if(!is.null(no.chrs.bearing.mut)){
    mutationCopyNumber = mutationCopyNumber / no.chrs.bearing.mut
    xlabel = "Fraction of Tumour Cells"
  }
  
  # X and Y coordinates of the density
  xx = density[,1]
  yy = density[,2]

  if(is.na(y.max)) { y.max=ceiling(max(polygon.data)) } #/10)*10
  
  # Plot the histogram, the density line and add the plot title
  par(mar = c(5,6,4,1)+0.1)
  if (is.null(mutationTypes)) {
    # Plot all the data as a lightgrey histogram
    hist(mutationCopyNumber[mutationCopyNumber<=x.max], breaks=seq(-0.1, x.max, 0.025), col="lightgrey",freq=FALSE, xlab=xlabel,main="", ylim=c(0,y.max),cex.axis=2,cex.lab=2)
  } else {
    # Plot SNVs and CNAs with different colours
    hist(mutationCopyNumber[mutationCopyNumber<=x.max & mutationTypes=="SNV"], breaks=seq(-0.1, x.max, 0.025), col=rgb(211/255,211/255,211/255,0.8),freq=FALSE, xlab=xlabel,main="", ylim=c(0,y.max),cex.axis=2,cex.lab=2)
    hist(mutationCopyNumber[mutationCopyNumber<=x.max & mutationTypes=="CNA"], breaks=seq(-0.1, x.max, 0.025), col=rgb(255/255,0,0,0.8), freq=FALSE, cex.axis=2, cex.lab=2, add=T)
  }
  polygon(c(xx, rev(xx)), polygon.data, border="plum4", col=cm.colors(1,alpha=0.3))
  lines(xx, yy, col="plum4", lwd=3)
  title(samplename, cex.main=3)
  
  # If cluster locations are provided, add them as a vertical line with nr of mutations mentioned
  if(!is.null(cluster.locations) & !is.null(mutation.assignments)) {
    assign.counts = table(mutation.assignments)
    
    for (cluster in unique(mutation.assignments)) {
      x = cluster.locations[cluster.locations[,1]==cluster, 2]
      lines(x=c(x, x), y=c(0, y.max), col="black", lwd=3)
      text(paste("Cluster",cluster, sep=" "), x=x+0.01, y=(9/10)*y.max, adj=c(0,0), cex=2)
      text(paste(as.numeric(assign.counts[names(assign.counts)==as.character(cluster)]), "mutations", sep=" "), x=x+0.01, y=(9/10)*y.max-0.35, adj=c(0,0), cex=2)
    }
  }
  
  # if(length(CALR)>0){
  #   x.index = sapply(1:length(CALR),function(i){which.min(abs(CALR[i]-xx))})
  #   CALR.yvals = (polygon.data[x.index]+polygon.data[2*length(xx)-x.index])/2
  #   points(CALR,CALR.yvals,pch=20,col="red",cex=3)
  # }
  
  if (!is.na(pngFile)) { dev.off() }
}

#' ggplot2 based density figure that can show CNA-pseudo-SNVs separately. Currently this figure only can plot the CCF space.
#' 
#' @param density Density as returned by Gibbs.subclone.density.est.1d.
#' @param polygon.data Polygon data as returned by Gibbs.subclone.density.est.1d.
#' @param mutationCopyNumber Mutation copy number matrix to be used to construct the CCF space.
#' @param no.chrs.bearing.mut Multiplicity matrix to be used to construct the CCF space.
#' @param pngFile Filename of the PNG file to be created. If NA is supplied the figure is returned (Default: NA).
#' @param density.from Legacy parameter that is no longer used (Default: 0).
#' @param x.max Max value on the x-axis to be plotted (Default: NA, then determined based on the data).
#' @param y.max Max value on the y-axis to be plotted (Default: NA, then determined based on the data).
#' @param y Legacy parameter that is used in the original figure to be able to plot different spaces. This value should be the number of supporting variant reads. It has no effect here. (Default: NULL).
#' @param N Legacy parameter that is used in the original figure to be able to plot different spaces. This value should be the total number of reads. It has no effect here. (Default: NULL).
#' @param samplename Name of the sample to be put on top of the figure.
#' @param CALR Legacy parameter that is no longer used.
#' @param cluster.locations Locations of where clusters were found. A vectical line is plotted for each cluster.
#' @param mutation.assignments Vector with hard assignment of mutations to clusters - optional, is used for plotting (Default: NULL)
#' @param mutationTypes Vector with mutation types. The background histogram can be coloured by mutation type (Default: SNV)
#' @param font_sizes A list that defines the size of the fonts in the title, axis and legends (Default: Best settings for 1500x1000 plot).
#' @param cbPalette Colours for mutation types: Grey for first mutation type (SNVs), orange for second (CNAs), as defined in LoadData
#' @author sd11
plot1D_2 = function(density, polygon.data, mutationCopyNumber, no.chrs.bearing.mut, pngFile=NA, density.from=0, x.max=NA, y.max=NA, y=NULL, N=NULL, samplename="", CALR=numeric(0), cluster.locations=NULL, mutation.assignments=NULL, mutationTypes=rep("SNV", length(mutationCopyNumber)), font_sizes = list(plot.title = 50, axis.text = 25, axis.title = 35, legend.text = 25, legend.title = 25, legend.position = "bottom"), cbPalette = c("lightgray", "red", "blue")) {
  colnames(density)[1] = "fraction.of.tumour.cells"
  
  conf.interval = data.frame(x=density[,1], ymax=(polygon.data[1:512] / sum(density$median.density)), ymin=(rev(polygon.data[513:1024]) / sum(density$median.density)))
  density$median.density = density$median.density / sum(density$median.density)
  ccf.df = as.data.frame(mutationCopyNumber / no.chrs.bearing.mut)
  ccf.df$mutationType = mutationTypes
  
  if (is.na(y.max)) { y.max=max(conf.interval$ymax) }
  if (is.na(x.max)) { x.max=ceiling(max(ccf.df)) }
  
  p = ggplot() +
    geom_histogram(data=ccf.df, mapping=aes(x=V1, y=(..count..)/sum(..count..), fill=mutationType, alpha=0.3), binwidth=0.025, position="stack", alpha=0.8, colour="black") +
    geom_ribbon(data=conf.interval, mapping=aes(x=x,ymin=ymin,ymax=ymax), fill=cm.colors(1, alpha=0.6)) +
    geom_line(data=density, mapping=aes(x=fraction.of.tumour.cells, y=median.density), colour="plum4", size=2) +
    xlab("Fraction of Tumour Cells") +
    ylab("Density") +
    ggtitle(samplename) +
    theme_bw() +
    xlim(0, x.max) +
    theme(axis.text=element_text(size=font_sizes$axis.text),
          axis.title=element_text(size=font_sizes$axis.title),
          plot.title=element_text(size=font_sizes$plot.title),
          legend.text=element_text(size=font_sizes$legend.text),
          legend.title=element_text(size=font_sizes$legend.title),
          legend.position=font_sizes$legend.position) +
    scale_fill_manual(values=cbPalette) +
    scale_colour_discrete(drop=F, limits=levels(ccf.df$mutationTypes))
  
  # If cluster locations are provided, add them as a vertical line with nr of mutations mentioned
  if(!is.null(cluster.locations) & !is.null(mutation.assignments)) {
    # Get non empty clusters and their ids
    clusters = unique(mutation.assignments)
    clusters = sort(clusters[!is.na(clusters)])
    assignment_counts = array(NA, length(clusters))
    for (c in 1:length(clusters)) {
      assignment_counts[c] = sum(mutation.assignments==clusters[c], na.rm=T)
    }
    dat = data.frame(cluster.locations[cluster.locations[,1] %in% clusters, c(1,2), drop=FALSE])
    colnames(dat) = c("non_empty_cluster_ids", "non_empty_cluster_locations")
    dat$assignment_counts = assignment_counts
    dat$y.max = y.max

    # Only attempt to plot when the lines will be within the bounds of the figure, otherwise this will crash
    if (any(dat$non_empty_cluster_locations < x.max)) {
      # Plot a line for each cluster, the cluster id and the number of mutations assigned to it
      p = p + geom_segment(data=dat, mapping=aes(x=non_empty_cluster_locations, xend=non_empty_cluster_locations, y=0, yend=y.max)) +
        geom_text(data=dat, mapping=aes(x=(non_empty_cluster_locations+0.01), y=(9/10)*y.max, label=paste("Cluster", dat$non_empty_cluster_ids, sep=" "), hjust=0), size=8) +
        geom_text(data=dat, mapping=aes(x=(non_empty_cluster_locations+0.01), y=(9/10)*y.max-((1/20)*y.max), label=paste(assignment_counts, "mutations", sep=" "), hjust=0), size=8)
    }
  }
  
  if (!is.na(pngFile)) {
    png(filename=pngFile,width=1500,height=1000)
    print(p)  
    dev.off()
  } else {
    return(p)
  }
}

#' Plot a table with the assignment counts
#' @param cluster_locations Cluster table with cluster number, cluster location and number of mutations as columns.
#' @param pngFile Output file to save the image.
#' @param cndata Optional CNA data table. This must be the table from after assigning CNA events to clusters (Default: NULL).
#' @param num_samples Optional parameter representing the number of samples that have been clustered (Default: 1)
#' @param indeldata Optional indel data table. This must be the table from after assigning indel events to clusters (Default: NULL).
#' @author sd11
plotAssignmentTable = function(cluster_locations, pngFile, cndata=NULL, num_samples=1, indeldata=NULL) {
  # Set the naming for the figure
  cluster_locations = as.data.frame(cluster_locations)
  colnames(cluster_locations) = c("Cluster", "Location", rep("", num_samples-1), "Num SNVs")
  # Order by ascending cluster number
  cluster_locations = cluster_locations[with(cluster_locations, order(-Cluster)),]
  
  for (i in 1:num_samples) {
    # First column is the cluster number, so take i+1
    cluster_locations[,i+1] = round(cluster_locations[,i+1], 2)
  }
  
  # Add in the indels if they are available
  if (!is.null(indeldata)) {
    cluster_locations$no.of.indels = 0
    for (cluster.no in unique(cluster_locations$Cluster)) {
      cluster_locations[cluster_locations$Cluster==cluster.no, "no.of.indels"] = sum(indeldata$cluster==cluster.no, na.rm=T)
    }
    colnames(cluster_locations)[ncol(cluster_locations)] = "Num indels"
  }
  
  # Add in the CNAs if they are available
  if (!is.null(cndata)) {
    cluster_locations$no.of.cnas = 0
    for (cluster.no in unique(cluster_locations$Cluster)) {
      cndata_cluster = cndata[cndata$cluster_assignment==cluster.no,]
      cndata_cluster = cndata_cluster[unique(cndata_cluster$startpos),]
      cluster_locations[cluster_locations$Cluster==cluster.no, "no.of.cnas"] = nrow(cndata_cluster)
    }
    colnames(cluster_locations)[ncol(cluster_locations)] = "Num CNAs"
  }
  
  png(filename=pngFile,width=500,height=500)
  grid.table(cluster_locations, rows=NULL)
  dev.off()
}

#' Create a 2D plot for multi-dimensional clustering
#' 
#' @param xvals Density coordinates of the x-axis
#' @param yvals Density coordinates of the y-axis
#' @param zvals The density itself
#' @param subclonal.fraction_x Fraction of cells estimate for each mutation in sample on x-axis
#' @param subclonal.fraction_y Fraction of cells estimate for each mutation in sample on y-axis
#' @param pngFile Filename to which the figure should be pushed
#' @param samplename_x Samplename to assign as x-axis label
#' @param samplename_y Samplename to assign as y-axis label
#' @param max.plotted.value The maximum value to plot on both x- and y-axis (Default: NA)
#' @param cluster.locations Locations of clusters to be overlayed (Default: NULL)
#' @param plot.mutations Boolean, supply TRUE when mutations should be overlayed on top of the density (Default: FALSE)
#' @author dw9
plotnD = function(xvals, yvals, zvals, subclonal.fraction_x, subclonal.fraction_y, pngFile, samplename_x, samplename_y, max.plotted.value=NA, cluster.locations=NULL, plot_mutations=F) {
  colours=colorRampPalette(c("white","red"))
  
  # Determine the minimum and maximum plotted values
  if (!is.na(max.plotted.value)) {
    range=list(c(-0.1,max.plotted.value), c(-0.1,max.plotted.value))
    # Subset x,y,z to remove data outside the plot
    zvals = data.matrix(zvals[xvals<=max.plotted.value,yvals<=max.plotted.value])
    xvals = data.matrix(xvals[xvals<=max.plotted.value])
    yvals = data.matrix(yvals[yvals<=max.plotted.value])
  } else {
    # Set dynamic range based on the subclonal fraction data
    range=list(c(floor(min(subclonal.fraction_x)*10)-1,ceiling(max(subclonal.fraction_x)*10)+1)/10, 
               c(floor(min(subclonal.fraction_y)*10)-1,ceiling(max(subclonal.fraction_y)*10)+1)/10)
  }

  #plot.data = cbind(data[[i]]$subclonal.fraction,data[[j]]$subclonal.fraction)
  plot.data = cbind(subclonal.fraction_x, subclonal.fraction_y)
  if(!is.na(max.plotted.value)) {
    plot.data = plot.data[plot.data[,1]<=max.plotted.value & plot.data[,2]<=max.plotted.value,]
  }
  
  if (!is.null(cluster.locations) & !plot_mutations) {
    # Plot the cluster locations on top of the density
    panel_function = function(...) { 
      lattice::panel.levelplot(...)
      lattice::panel.abline(h = 0:floor(max(plot.data[,2])))
      lattice::panel.abline(v = 0:floor(max(plot.data[,1])))
      lattice::lpoints(cluster.locations, pch=".", cex=6, col="black")
    }
  } else if (plot_mutations) {
    # Plot the mutations overlayed on top of the density
    panel_function = function(...) { 
      lattice::panel.levelplot(...)
      lattice::panel.abline(h = 0:floor(max(plot.data[,2])))
      lattice::panel.abline(v = 0:floor(max(plot.data[,1])))                   
      lattice::lpoints(plot.data, pch=".", cex=4, col="black") 
    }
  } else {
    # No overlay to be plotted
    panel_function = function(...) { 
      lattice::panel.levelplot(...)
      lattice::panel.abline(h = 0:floor(max(plot.data[,2])))
      lattice::panel.abline(v = 0:floor(max(plot.data[,1])))
    }
  }
  
  png(filename=pngFile, width=1500, height=1000) 
  image.wid = 500 * (range[[1]][2] - range[[1]][1])
  image.ht = 500 * (range[[2]][2] - range[[2]][1])
  fig = lattice::levelplot(zvals,
                row.values=xvals,
                column.values=yvals,
                xlim=range[[1]],
                ylim=range[[2]],
                xlab=list(label=samplename_x,cex=2),
                ylab=list(label=samplename_y,cex=2),
                scales=list(x=list(cex=1.5),y=list(cex=1.5)),
                col.regions=colours,
                colorkey=F,
                panel=panel_function
  )    
  print(fig)
  dev.off()
}
