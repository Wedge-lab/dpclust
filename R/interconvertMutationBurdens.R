#
# Core conversion functions
#

mutationBurdenToMutationCopyNumber<-function(burden,totalCopyNumber,cellularity,normalCopyNumber = rep(2,length(burden))){
	mutCopyNumber = burden/cellularity*(cellularity*totalCopyNumber+normalCopyNumber*(1-cellularity))
	mutCopyNumber[is.nan(mutCopyNumber)]=0
	return(mutCopyNumber)
}

mutationCopyNumberToMutationBurden<-function(copyNumber,totalCopyNumber,cellularity,normalCopyNumber = rep(2,length(copyNumber))){
	burden = copyNumber*cellularity/(cellularity*totalCopyNumber+normalCopyNumber*(1-cellularity))
	burden[is.nan(burden)|(burden<0.000001)]=0.000001
	burden[burden>0.999999]=0.999999
	return(burden)	
}
