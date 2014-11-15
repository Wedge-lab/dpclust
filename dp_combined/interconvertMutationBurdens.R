mutationBurdenToMutationCopyNumber<-function(burden,totalCopyNumber,cellularity,normalCopyNumber = rep(2,length(burden))){
	#mutCopyNumber = burden/cellularity*(cellularity*(totalCopyNumber-burden)+2*(1-cellularity))
	mutCopyNumber = burden/cellularity*(cellularity*totalCopyNumber+normalCopyNumber*(1-cellularity))
	mutCopyNumber[is.nan(mutCopyNumber)]=0
	return(mutCopyNumber)
}

#useful for chrX
mutationCopyNumberToMutationBurden<-function(copyNumber,totalCopyNumber,cellularity,normalCopyNumber = rep(2,length(copyNumber))){
	#burden = copyNumber*cellularity/(cellularity*(totalCopyNumber-burden)+2*(1-cellularity))
	burden = copyNumber*cellularity/(cellularity*totalCopyNumber+normalCopyNumber*(1-cellularity))
	#burden[is.nan(burden)]=0
	burden[is.nan(burden)|(burden<0.000001)]=0.000001
	#190412
	#burden[burden>1]=1
	burden[burden>0.999999]=0.999999
	return(burden)	
}
