
oneDimensionalClustering <- function(samplename, subclonal.fraction, GS.data, density, no.iters, no.iters.burn.in) {
  no.muts = length(subclonal.fraction)
  normal.copy.number = rep(2,no.muts)
  
  S.i = GS.data$S.i
  V.h = GS.data$V.h
  pi.h = GS.data$pi.h[,,1]
  
  hypercube.size = 20
  localOptima = NULL
  peak.indices = NULL
  for(i in (1+hypercube.size):(nrow(density)-hypercube.size)){
  	if(density$median.density[i] == max(density$median.density[(i-hypercube.size):(i+hypercube.size)])){
  	  localOptima = c(localOptima,density$fraction.of.tumour.cells[i])
  		peak.indices = c(peak.indices,i)
  	}
  }
  
  print("localOptima")
  print(localOptima)
  
  write.table(localOptima,paste(samplename,"_localOptima.txt",sep=""),quote=F,sep="\t")
  
  #ASSIGN mutations to clusters
  no.optima = length(localOptima)
  if(no.optima>1){	
  	boundary = array(NA,no.optima-1)
  	mutation.preferences = array(0,c(no.muts,no.optima))
  	for(i in 1:(no.optima-1)){
  		min.density = min(density$median.density[(peak.indices[i]+1):(peak.indices[i+1]-1)])
  		min.indices = intersect(which(density$median.density == min.density),(peak.indices[i]+1):(peak.indices[i+1]-1))
  		
  		#what distance along the line between a pair of optima do we have to go to reach the minimum density
      boundary[i] = (density$fraction.of.tumour.cells[max(min.indices)] + density$fraction.of.tumour.cells[min(min.indices)])/2
  	}
  
    sampledIters = (no.iters.burn.in + 1) : no.iters
  	#don't use the intitial state
  	sampledIters = sampledIters[sampledIters!=1]	
  	if(length(sampledIters) > 1000){
  	  sampledIters=floor(post.burn.in.start + (1:1000) * (no.iters - no.iters.burn.in)/1000)  
  	}
  
  	S.i = data.matrix(S.i)
  	for(s in sampledIters){
  		temp.preferences = array(0,c(no.muts,no.optima))
  		for(c in unique(S.i[s,])){
  			bestOptimum = sum(pi.h[s,c]>boundary)+1
  			temp.preferences[S.i[s,]==c,bestOptimum] = temp.preferences[S.i[s,]==c,bestOptimum] + 1		
  		}
      iter.preferences = t(apply(temp.preferences, 1, function(p) { as.integer(p==max(p)) / sum(p==max(p)) }))
  		mutation.preferences = mutation.preferences + iter.preferences
  	}
  	mutation.preferences = mutation.preferences / length(sampledIters)
    most.likely.cluster = max.col(mutation.preferences)
  	
    out = cbind(mutation.preferences,most.likely.cluster)
  	colnames(out)[(ncol(out)-no.optima):ncol(out)] = c(paste("prob.cluster",1:no.optima,sep=""),"most.likely.cluster")
  	write.table(cbind(1:no.optima,localOptima,colSums(mutation.preferences)),paste(samplename,"_optimaInfo.txt",sep=""),col.names = c("cluster.no","location","no.of.mutations"),row.names=F,sep="\t",quote=F)		
  	write.table(out,paste(samplename,"_DP_and_cluster_info.txt",sep=""),sep="\t",row.names=F,quote=F)
    
  	most.likely.cluster.likelihood = apply(mutation.preferences, 1, max)
  
  }else{
  	most.likely.cluster = rep(1,no.muts)
  	most.likely.cluster.likelihood = rep(1,no.muts)
  }
  
  return(list(best.node.assignments=most.likely.cluster, best.assignment.likelihoods=most.likely.cluster.likelihood))
}


# twoDimensionalClustering <- function(samplename, subclonal.fraction, GS.data, density, no.iters, no.iters.burn.in) {
#   
#   hypercube.size = 20
#   localOptima.x = NULL
#   localOptima.y = NULL
#   peak.indices.x = NULL
#   peak.indices.y = NULL
#   # Obtain density for each point in the data.frame
#   for(i in (1+hypercube.size):(nrow(density$median.density)-hypercube.size)){
#     for(j in (1+hypercube.size):(ncol(density$median.density)-hypercube.size)){
#       if(density$median.density[i,j] == max(density$median.density[(i-hypercube.size):(i+hypercube.size), (j-hypercube.size):(j+hypercube.size)])){
#         # Only take those clusters that have a non-trivial density
#         if(density$median.density[i,j] > 0.1) { 
#           localOptima.x = c(localOptima.x,density$xvals[i])
#           localOptima.y = c(localOptima.y,density$yvals[j])
#           peak.indices.x = c(peak.indices.x,i)
#           peak.indices.y = c(peak.indices.y,j)
#         }
#       }
#     }
#   }
# 
#   localOptima = data.frame(x=localOptima.x, y=localOptima.y)
#   print("localOptima")
#   print(localOptima)
#   
#   write.table(localOptima,paste(samplename,"_localOptima.txt",sep=""),quote=F,sep="\t")
#   
#   # Below is not adapted to 2D case
# #   no.optima = nrow(localOptima)
# #   if(no.optima>1){  
# #     boundary = array(NA,no.optima-1)
# #     mutation.preferences = array(0,c(no.muts,no.optima))
# #     for(i in 1:(no.optima-1)){
# #       min.density = min(density$median.density[(peak.indices.x[i]+1):(peak.indices.x[i+1]-1), (peak.indices.y[i]+1):(peak.indices.y[i+1]-1)])
# #       min.indices = intersect(which(density$median.density == min.density),(peak.indices[i]+1):(peak.indices[i+1]-1))
# #       
# #       #what distance along the line between a pair of optima do we have to go to reach the minimum density
# #       boundary[i] = (density$fraction.of.tumour.cells[max(min.indices)] + density$fraction.of.tumour.cells[min(min.indices)])/2
# #     }
# #     
# #     sampledIters = (no.iters.burn.in + 1) : no.iters
# #     #don't use the intitial state
# #     sampledIters = sampledIters[sampledIters!=1]	
# #     if(length(sampledIters) > 1000){
# #       sampledIters=floor(post.burn.in.start + (1:1000) * (no.iters - no.iters.burn.in)/1000)  
# #     }
# #     
# #     S.i = data.matrix(S.i)
# #     for(s in sampledIters){
# #       temp.preferences = array(0,c(no.muts,no.optima))
# #       for(c in unique(S.i[s,])){
# #         bestOptimum = sum(pi.h[s,c]>boundary)+1
# #         temp.preferences[S.i[s,]==c,bestOptimum] = temp.preferences[S.i[s,]==c,bestOptimum] + 1		
# #       }
# #       iter.preferences = t(apply(temp.preferences, 1, function(p) { as.integer(p==max(p)) / sum(p==max(p)) }))
# #       mutation.preferences = mutation.preferences + iter.preferences
# #     }
# #     mutation.preferences = mutation.preferences / length(sampledIters)
# #     most.likely.cluster = max.col(mutation.preferences)
# #     
# #     out = cbind(mutation.preferences,most.likely.cluster)
# #     colnames(out)[(ncol(out)-no.optima):ncol(out)] = c(paste("prob.cluster",1:no.optima,sep=""),"most.likely.cluster")
# #     write.table(cbind(1:no.optima,localOptima,colSums(mutation.preferences)),paste(samplename,"_optimaInfo.txt",sep=""),col.names = c("cluster.no","MCN","no.of.mutations"),row.names=F,sep="\t",quote=F)		
# #     write.table(out,paste(samplename,"_DP_and_cluster_info.txt",sep=""),sep="\t",row.names=F,quote=F)
# #     
# #     most.likely.cluster.likelihood = apply(mutation.preferences, 1, max)
# #     
# #   }else{
# #     most.likely.cluster = rep(1,no.muts)
# #     most.likely.cluster.likelihood = rep(1,no.muts)
# #   }
# }


# q(save="no")