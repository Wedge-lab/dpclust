mutationBurdenToMutationCopyNumber<-function(burden,totalCopyNumber,cellularity){
	mutCopyNumber = burden/cellularity*(cellularity*totalCopyNumber+2*(1-cellularity))
	mutCopyNumber[is.nan(mutCopyNumber)]=0
	return(mutCopyNumber)
}

mutationCopyNumberToMutationBurden<-function(copyNumber,totalCopyNumber,cellularity){
	burden = copyNumber*cellularity/(cellularity*totalCopyNumber+2*(1-cellularity))
	burden[is.nan(burden)|(burden<0.000001)]=0.000001
	burden[burden>0.999999]=0.999999
	return(burden)	
}