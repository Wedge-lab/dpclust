library(lattice)

plot1D = function(density, polygon.data, pngFile=NA, density.from=0, x.max=NA, y.max=NA, y=NULL, N=NULL, mutationCopyNumber=NULL, no.chrs.bearing.mut=NULL,samplename="",CALR=numeric(0), cluster.locations=NULL, mutation.assignments=NULL) {
  #
  # Creates the 1D density plot in either allele fraction, mutation copy number or fraction of tumour cells space.
  #   density:      Is a two column data frame. First column x-axis density line, second column y-axis density line.
  #   polygon.data: Contains the coordinates of the 95% confidence interval around the density line.
  #   pngFile:      Where to save the figure.
  #   density.from: Startpoint from where density is drawn.
  #   x.max/y.max:  Maximum values on x and y axis in the plot.
  #   y:            MutCount data.
  #   N:            Total count (mutCount+WTCount).
  #   mutationCopyNumber:   Mutation copy number data.
  #   no.chrs.bearing.mut:  Copynumber adjustment per mutation.
  #   samplename:   Name of the sample under analysis.
  #   CALR:         ?? some kind of annotation
  #   cluster.locations Locations in the space to be plotted where clusters reside.
  #
  # ggplot equivalent: 
  # conf.interval = data.frame(x=c(density[,1], rev(density[,1])), y=as.vector(polygon.data[,1]))
  # mutationCopyNumber.df = as.data.frame(mutationCopyNumber)
  # p = ggplot() + geom_histogram(data=mutationCopyNumber.df, mapping=aes(x=V1, y=..density..), binwidth=0.1) + geom_polygon(data=conf.interval, mapping=aes(x=x, y=y), fill='lightblue', alpha=0.7) + geom_line(data=density, mapping=aes(x=fraction.of.tumour.cells, y=median.density), colour="purple")
  #
  if (!is.na(pngFile)) { png(filename=pngFile,,width=1500,height=1000) }
  
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
  hist(mutationCopyNumber[mutationCopyNumber<=x.max], breaks=seq(-0.1, x.max, 0.025), col="lightgrey",freq=FALSE, xlab=xlabel,main="", ylim=c(0,y.max),cex.axis=2,cex.lab=2)
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
      text(paste(as.numeric(assign.counts[cluster]), "mutations", sep=" "), x=x+0.01, y=(9/10)*y.max-0.35, adj=c(0,0), cex=2)
    }
  }
  
  if(length(CALR)>0){
    x.index = sapply(1:length(CALR),function(i){which.min(abs(CALR[i]-xx))})
    CALR.yvals = (polygon.data[x.index]+polygon.data[2*length(xx)-x.index])/2
    points(CALR,CALR.yvals,pch=20,col="red",cex=3)
  }
  
  if (!is.na(pngFile)) { dev.off() }
}

plotnD = function(xvals, yvals, zvals, subclonal.fraction_x, subclonal.fraction_y, pngFile, samplename_x, samplename_y, max.plotted.value=NA) {
  #  
  # Create a 2D density plot for the nD clustering.
  #   xvals:                Density coordinates of the x-axis
  #   yvals:                Density coordinates of the y-axis
  #   zvals:                The density itself
  #   subclonal.fraction_x: Fraction of cells estimate for each mutation in sample on x-axis
  #   subclonal.fraction_y: Fraction of cells estimate for each mutation in sample on y-axis
  #   pngFile:              Filename to which the figure should be pushed
  #   samplename_x:         Samplename to assign as x-axis label
  #   samplename_y:         Samplename to assign as y-axis label
  #   max.plotted.value:    The maximum value to plot on both x- and y-axis
  #
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
  
  # First plot without mutations
  png(filename=gsub(".png","_withoutMutations.png",pngFile),width=1500,height=1000)       
  image.wid = 500 * (range[[1]][2] - range[[1]][1])
  image.ht = 500 * (range[[2]][2] - range[[2]][1])
  fig=levelplot(zvals,row.values=xvals,column.values=yvals,xlim=range[[1]],ylim=range[[2]],xlab=list(label=samplename_x,cex=2),ylab=list(label=samplename_y,cex=2),scales=list(x=list(cex=1.5),y=list(cex=1.5)),col.regions=colours,colorkey=F,
                panel = function(...) { 
                  panel.levelplot(...)
                  panel.abline(h = 0:floor(max(plot.data[,2])))
                  panel.abline(v = 0:floor(max(plot.data[,1])))
                }
  )    
  print(fig)
  dev.off()
  
  # Second plot with mutations
  png(filename=gsub(".png","_withMutations.png",pngFile),width=1500,height=1000)
  image.wid = 500 * (range[[1]][2] - range[[1]][1])
  image.ht = 500 * (range[[2]][2] - range[[2]][1])
  fig=levelplot(zvals,row.values=xvals,column.values=yvals,xlim=range[[1]],ylim=range[[2]],xlab=list(label=samplename_x,cex=2),ylab=list(label=samplename_y,cex=2),scales=list(x=list(cex=1.5),y=list(cex=1.5)),col.regions=colours,colorkey=F,
                panel = function(...) { 
                  panel.levelplot(...)
                  panel.abline(h = 0:floor(max(plot.data[,2])))
                  panel.abline(v = 0:floor(max(plot.data[,1])))                   
                  lpoints(plot.data,pch=".",cex=6,col="black") 
                  # Left over from previous code: differentiate between how well mutations are covered
                  #if(nrow(burden>=500)){
                  #  lpoints(burden,pch=".",cex=1,col="black")
                  #}else if(nrow(burden>=100)){
                  #  lpoints(burden,pch=".",cex=2,col="black")
                  #}else{
                  #  lpoints(burden,pch=".",cex=4,col="black")
                  #}
                }
  )
  print(fig)
  dev.off()
      
}
