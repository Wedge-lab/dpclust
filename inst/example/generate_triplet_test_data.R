# simulation of a simple case with 6 samples

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

coverage = 100
n_muts = 100

generate_clonal_cluster = function(coverage, n_muts, ccf) {
  coverage = rpois(n_muts, lambda=coverage)
  mut.count = rbinom(coverage, n_muts, ccf/2)
  WT.count = coverage - mut.count
  vaf = mut.count / (mut.count + WT.count)
  
  # chr	end	WT.count	mut.count	subclonal.CN	mutation.copy.number	subclonal.fraction	no.chrs.bearing.mut	phase
  output = data.frame(WT.count=WT.count, 
                      mut.count=mut.count, 
                      subclonal.CN=2,
                      mutation.copy.number=mutationBurdenToMutationCopyNumber(vaf, 2, 1),
                      subclonal.fraction=mutationBurdenToMutationCopyNumber(vaf, 2, 1),
                      no.chrs.bearing.mut=1,
                      phase="unphased")
  return(output)
}

# cluster 1: clonal, present in 6 samples
sample_01 = generate_clonal_cluster(coverage, n_muts, 1)
sample_02 = generate_clonal_cluster(coverage, n_muts, 1)
sample_03 = generate_clonal_cluster(coverage, n_muts, 1)
sample_04 = generate_clonal_cluster(coverage, n_muts, 1)
sample_05 = generate_clonal_cluster(coverage, n_muts, 1)
sample_06 = generate_clonal_cluster(coverage, n_muts, 1)

# cluster 2: clonal, present in 3 samples
sample_01 = rbind(sample_01, generate_clonal_cluster(coverage, n_muts, 1))
sample_02 = rbind(sample_02, generate_clonal_cluster(coverage, n_muts, 1))
sample_03 = rbind(sample_03, generate_clonal_cluster(coverage, n_muts, 1))
sample_04 = rbind(sample_04, generate_clonal_cluster(coverage, n_muts, 0))
sample_05 = rbind(sample_05, generate_clonal_cluster(coverage, n_muts, 0))
sample_06 = rbind(sample_06, generate_clonal_cluster(coverage, n_muts, 0))

# cluster 3: clonal, present in 3 samples
sample_01 = rbind(sample_01, generate_clonal_cluster(coverage, n_muts, 0))
sample_02 = rbind(sample_02, generate_clonal_cluster(coverage, n_muts, 0))
sample_03 = rbind(sample_03, generate_clonal_cluster(coverage, n_muts, 0))
sample_04 = rbind(sample_04, generate_clonal_cluster(coverage, n_muts, 1))
sample_05 = rbind(sample_05, generate_clonal_cluster(coverage, n_muts, 1))
sample_06 = rbind(sample_06, generate_clonal_cluster(coverage, n_muts, 1))

# cluster 4: clonal, private in 1 sample
sample_01 = rbind(sample_01, generate_clonal_cluster(coverage, n_muts, 1))
sample_02 = rbind(sample_02, generate_clonal_cluster(coverage, n_muts, 0))
sample_03 = rbind(sample_03, generate_clonal_cluster(coverage, n_muts, 0))
sample_04 = rbind(sample_04, generate_clonal_cluster(coverage, n_muts, 0))
sample_05 = rbind(sample_05, generate_clonal_cluster(coverage, n_muts, 0))
sample_06 = rbind(sample_06, generate_clonal_cluster(coverage, n_muts, 0))

# cluster 5: clonal, private in 1 sample
sample_01 = rbind(sample_01, generate_clonal_cluster(coverage, n_muts, 0))
sample_02 = rbind(sample_02, generate_clonal_cluster(coverage, n_muts, 1))
sample_03 = rbind(sample_03, generate_clonal_cluster(coverage, n_muts, 0))
sample_04 = rbind(sample_04, generate_clonal_cluster(coverage, n_muts, 0))
sample_05 = rbind(sample_05, generate_clonal_cluster(coverage, n_muts, 0))
sample_06 = rbind(sample_06, generate_clonal_cluster(coverage, n_muts, 0))

# cluster 6: clonal, private in 1 sample
sample_01 = rbind(sample_01, generate_clonal_cluster(coverage, n_muts, 0))
sample_02 = rbind(sample_02, generate_clonal_cluster(coverage, n_muts, 0))
sample_03 = rbind(sample_03, generate_clonal_cluster(coverage, n_muts, 1))
sample_04 = rbind(sample_04, generate_clonal_cluster(coverage, n_muts, 0))
sample_05 = rbind(sample_05, generate_clonal_cluster(coverage, n_muts, 0))
sample_06 = rbind(sample_06, generate_clonal_cluster(coverage, n_muts, 0))

# cluster 7: clonal, private in 1 sample
sample_01 = rbind(sample_01, generate_clonal_cluster(coverage, n_muts, 0))
sample_02 = rbind(sample_02, generate_clonal_cluster(coverage, n_muts, 0))
sample_03 = rbind(sample_03, generate_clonal_cluster(coverage, n_muts, 0))
sample_04 = rbind(sample_04, generate_clonal_cluster(coverage, n_muts, 1))
sample_05 = rbind(sample_05, generate_clonal_cluster(coverage, n_muts, 0))
sample_06 = rbind(sample_06, generate_clonal_cluster(coverage, n_muts, 0))

# cluster 8: clonal, private in 1 sample
sample_01 = rbind(sample_01, generate_clonal_cluster(coverage, n_muts, 0))
sample_02 = rbind(sample_02, generate_clonal_cluster(coverage, n_muts, 0))
sample_03 = rbind(sample_03, generate_clonal_cluster(coverage, n_muts, 0))
sample_04 = rbind(sample_04, generate_clonal_cluster(coverage, n_muts, 0))
sample_05 = rbind(sample_05, generate_clonal_cluster(coverage, n_muts, 1))
sample_06 = rbind(sample_06, generate_clonal_cluster(coverage, n_muts, 0))

# cluster 9: clonal, private in 1 sample
sample_01 = rbind(sample_01, generate_clonal_cluster(coverage, n_muts, 0))
sample_02 = rbind(sample_02, generate_clonal_cluster(coverage, n_muts, 0))
sample_03 = rbind(sample_03, generate_clonal_cluster(coverage, n_muts, 0))
sample_04 = rbind(sample_04, generate_clonal_cluster(coverage, n_muts, 0))
sample_05 = rbind(sample_05, generate_clonal_cluster(coverage, n_muts, 0))
sample_06 = rbind(sample_06, generate_clonal_cluster(coverage, n_muts, 1))

# add chr/pos
chrpos = data.frame(chr=1, end=1:nrow(sample_01))
sample_01 = cbind(chrpos, sample_01)
sample_02 = cbind(chrpos, sample_02)
sample_03 = cbind(chrpos, sample_03)
sample_04 = cbind(chrpos, sample_04)
sample_05 = cbind(chrpos, sample_05)
sample_06 = cbind(chrpos, sample_06)

# write data files
if (!dir.exists("Data")) { dir.create("Data") }
write.table(sample_01, file="Data/simulated_002_01_dp_input.txt", quote=F, sep="\t", row.names=F)
write.table(sample_02, file="Data/simulated_002_02_dp_input.txt", quote=F, sep="\t", row.names=F)
write.table(sample_03, file="Data/simulated_002_03_dp_input.txt", quote=F, sep="\t", row.names=F)
write.table(sample_04, file="Data/simulated_002_04_dp_input.txt", quote=F, sep="\t", row.names=F)
write.table(sample_05, file="Data/simulated_002_05_dp_input.txt", quote=F, sep="\t", row.names=F)
write.table(sample_06, file="Data/simulated_002_06_dp_input.txt", quote=F, sep="\t", row.names=F)

# create a master file
masterfile = data.frame()
for (i in 1:6) {
  masterfile = rbind(masterfile, 
                     data.frame(sample="simulated_002", subsample=paste0("0", i), datafile=paste0("simulated_002_0", i, "_dp_input.txt"), cellularity=1))
}
write.table(masterfile, file="simulated_triplet.txt", quote=F, sep="\t", row.names=F)



