library(GenomicRanges)

args=commandArgs(TRUE)
#                                                  Battenberg output                                GetAlleleFrequencies input  Battenberg output                               GetAlleleFrequencis output
# R CMD BATCH '--no-restore-data --args PD7404a ../../battenberg/PD7404a/PD7404a_rho_and_psi.txt ../mutation_loci/PD7404a.loci ../../battenberg/PD7404a/PD7404a_subclones.txt PD7404a_alleleFrequencies.txt female' ~/repo/dirichlet/GetDirichletProcessInfo.R PD7404a.Rout


# samplename = "cf583bf4-a5aa-4c7a-9962-1d20ffc42a98"
# cellularity_file = "/lustre/scratch112/sanger/cgppipe/PanCancerDownloads/workspace/sd11/variant_calling_pilot_50/battenberg/cf583bf4-a5aa-4c7a-9962-1d20ffc42a98_rho_and_psi.txt"
# loci_file = "/lustre/scratch112/sanger/cgppipe/PanCancerDownloads/workspace/sd11/variant_calling_pilot_50/haplotype/mutation_loci/cf583bf4-a5aa-4c7a-9962-1d20ffc42a98.loci"
# subclone_file = "/lustre/scratch112/sanger/cgppipe/PanCancerDownloads/workspace/sd11/variant_calling_pilot_50/battenberg/cf583bf4-a5aa-4c7a-9962-1d20ffc42a98_subclones.txt"
# allele_frequencies_file = "/lustre/scratch112/sanger/cgppipe/PanCancerDownloads/workspace/sd11/variant_calling_pilot_50/haplotype/mutation_loci/cf583bf4-a5aa-4c7a-9962-1d20ffc42a98_alleleFrequencies.txt"
# sex = "male"
# outdir = "test/" #/lustre/scratch112/sanger/cgppipe/PanCancerDownloads/workspace/sd11/variant_calling_pilot_50/dirichlet_input/

# samplename = "5c7b3dc8-f64b-4ca9-bd43-2245297450b4"
# cellularity_file = "/lustre/scratch112/sanger/cgppipe/PanCancerDownloads/workspace/sd11/variant_calling_pilot_50/battenberg/5c7b3dc8-f64b-4ca9-bd43-2245297450b4_rho_and_psi.txt"
# loci_file = "/lustre/scratch112/sanger/cgppipe/PanCancerDownloads/workspace/sd11/variant_calling_pilot_50/haplotype/mutation_loci/5c7b3dc8-f64b-4ca9-bd43-2245297450b4.loci" 
# subclone_file = "/lustre/scratch112/sanger/cgppipe/PanCancerDownloads/workspace/sd11/variant_calling_pilot_50/battenberg/5c7b3dc8-f64b-4ca9-bd43-2245297450b4_subclones.txt"
# allele_frequencies_file = "/lustre/scratch112/sanger/cgppipe/PanCancerDownloads/workspace/sd11/variant_calling_pilot_50/haplotype/mutation_loci/5c7b3dc8-f64b-4ca9-bd43-2245297450b4_alleleFrequencies.txt"
# sex = "male"
# outdir = "test/" #/lustre/scratch112/sanger/cgppipe/PanCancerDownloads/workspace/sd11/variant_calling_pilot_50/dirichlet_input/

samplename = toString(args[1])
cellularity_file = toString(args[2])
loci_file = toString(args[3])
subclone_file = toString(args[4])
allele_frequencies_file = toString(args[5])
sex = toString(args[6])
outdir = toString(args[7])
# TODO: implement additional files here

wd = getwd()
setwd("~/repo/dirichlet/dp_combined")
source("interconvertMutationBurdens.R")
setwd(wd)

if(sex == 'male' | sex == 'Male') {
  isMale = T
} else if(sex == 'female' | sex == 'Female') {
  isMale = F
} else {
  stop("Unknown sex supplied, exit.")
}

##################################
# See script below these functions
##################################

GetDirichletProcessInfo<-function(samplename, cellularity, info, subclone.file, is.male = F, out.dir = NULL, phase.dir = NULL, SNP.phase.file = NULL, mut.phase.file = NULL){
	print(samplename)

	subclone.data = read.table(subclone.file,sep="\t",header=T,row.names=1, stringsAsFactors=F)
# 	subclone.data$subclonal.CN = (subclone.data$nMaj1_A + subclone.data$nMin1_A) * subclone.data$frac1_A
  subclone.data.gr = GRanges(subclone.data$chr, IRanges(subclone.data$startpos, subclone.data$endpos), rep('*', nrow(subclone.data)))
  elementMetadata(subclone.data.gr) = subclone.data[,3:ncol(subclone.data)]
# 	info2 = as.data.frame(cbind(as.data.frame(info), array(NA, c(length(info), 7))))
#   colnames(info2) = c('chr','pos','WT.count','mut.count','subclonal.CN','nMaj1','nMin1','frac1','nMaj2','nMin2','frac2')

  info_anno = as.data.frame(cbind(array(NA, c(length(info), 7)))) 
  colnames(info_anno) = c('subclonal.CN','nMaj1','nMin1','frac1','nMaj2','nMin2','frac2')
  inds = findOverlaps(info, subclone.data.gr)  
  info_anno[queryHits(inds),2:7] = subclone.data[subjectHits(inds),][,c("nMaj1_A", "nMin1_A", "frac1_A", "nMaj2_A", "nMin2_A", "frac2_A")]

  CN1 = (info_anno[queryHits(inds),]$nMaj1 + info_anno[queryHits(inds),]$nMin1) * info_anno[queryHits(inds),]$frac1
  # If frac is not one for allele 1 (i.e. not only CN data for allele 1), add the CN contribution of allele 2 as well
  CN2 = (info_anno[queryHits(inds),]$nMaj2 + info_anno[queryHits(inds),]$nMin2) * info_anno[queryHits(inds),]$frac2 * ifelse(info_anno[queryHits(inds),]$frac1 != 1, 1, 0)
  CN2[is.na(CN2)] = 0
  info_anno[queryHits(inds),]$subclonal.CN = CN1 + CN2

  elementMetadata(info) = cbind(elementMetadata(info), info_anno)

#THIS ASSUMES that info contains chr and pos in the first 2 columns
# 	for(r in 1:nrow(subclone.data)){
# 		CN = (subclone.data$nMaj1_A[r] + subclone.data$nMin1_A[r]) * subclone.data$frac1_A[r]
# 		if(subclone.data$frac1_A[r] != 1){
# 			CN = CN + (subclone.data$nMaj2_A[r] + subclone.data$nMin2_A[r]) * subclone.data$frac2_A[r]
# 		}
# 		print(CN)
# 		info2$subclonal.CN[info2[,1]==subclone.data$chr[r] & info2[,2]>=subclone.data$startpos[r] & info2[,2]<=subclone.data$endpos[r]] = CN
# 		info2$nMaj1[info2[,1]==subclone.data$chr[r] & info2[,2]>=subclone.data$startpos[r] & info2[,2]<=subclone.data$endpos[r]] = subclone.data$nMaj1_A[r]
# 		info2$nMin1[info2[,1]==subclone.data$chr[r] & info2[,2]>=subclone.data$startpos[r] & info2[,2]<=subclone.data$endpos[r]] = subclone.data$nMin1_A[r]
# 		info2$frac1[info2[,1]==subclone.data$chr[r] & info2[,2]>=subclone.data$startpos[r] & info2[,2]<=subclone.data$endpos[r]] = subclone.data$frac1_A[r]
# 		info2$nMaj2[info2[,1]==subclone.data$chr[r] & info2[,2]>=subclone.data$startpos[r] & info2[,2]<=subclone.data$endpos[r]] = subclone.data$nMaj2_A[r]
# 		info2$nMin2[info2[,1]==subclone.data$chr[r] & info2[,2]>=subclone.data$startpos[r] & info2[,2]<=subclone.data$endpos[r]] = subclone.data$nMin2_A[r]
# 		info2$frac2[info2[,1]==subclone.data$chr[r] & info2[,2]>=subclone.data$startpos[r] & info2[,2]<=subclone.data$endpos[r]] = subclone.data$frac2_A[r]	
# 	}
#   print(info2[info2[,2]=="81432939",])

	info$phase="unphased"
	if(!is.null(phase.dir)){
		for(chr in unique(info[,1])){
			chr.inds = which(info[,1]==chr)
			ph.file = paste(phase.dir,"/",samplename,"_muts_linkedMuts_segmented_chr",chr,".txt",sep="")
			if(file.exists(ph.file)){
				phase.info = read.table(ph.file,sep="\t",header=F,row.names=NULL,stringsAsFactors=F)
				indices = match(info[chr.inds,2],phase.info[,2])
				info$phase[chr.inds[!is.na(indices)]] = phase.info[indices[!is.na(indices)],14]
			}
		}
	}
	info$phase[is.na(info$phase)]="unphased"

	if(is.male & "chr" %in% names(info)){
		normal.CN = rep(2,nrow(info))
		normal.CN[info$chr=="X"| info$chr=="Y"] = 1
		info$mutation.copy.number = mutationBurdenToMutationCopyNumber(info$mut.count/ (info$mut.count + info$WT.count) , info$subclonal.CN, cellularity, normal.CN)
	}else{
		info$mutation.copy.number = mutationBurdenToMutationCopyNumber(info$mut.count/ (info$mut.count + info$WT.count) , info$subclonal.CN, cellularity)
	}

	#convert MCN to subclonal fraction - tricky for amplified mutations
	info$subclonal.fraction = info$mutation.copy.number
	expected.burden.for.MCN = mutationCopyNumberToMutationBurden(rep(1,length(info)),info$subclonal.CN,cellularity)
	non.zero.indices = which(info$mut.count>0 & !is.na(expected.burden.for.MCN))
	#test for mutations in more than 1 copy
	p.vals = sapply(1:length(non.zero.indices),function(v,e,i){
    prop.test(v$mut.count[i],v$mut.count[i] + v$WT.count[i],e[i],alternative="greater")$p.value 
    },v=info[non.zero.indices,], e=expected.burden.for.MCN[non.zero.indices])
	amplified.muts = non.zero.indices[p.vals<=0.05]

	info$no.chrs.bearing.mut = 1	

	#copy numbers of subclones can only differ by 1 or 0 (as assumed when calling subclones)
	if(length(amplified.muts)>0){		
		for(a in 1:length(amplified.muts)){
			max.CN2=0
			#use phasing info - if on 'deleted' (lower CN chromosome), use the minor copy number
			if(info$phase[amplified.muts[a]]=="MUT_ON_DELETED"){
				print("mut on minor chromosome")
				max.CN1 = info$nMin1[amplified.muts[a]]
				frac1 = info$frac1[amplified.muts[a]]
				frac2=0
				if(!is.na(info$nMin2[amplified.muts[a]])){
					#swap subclones, so that the one with the higher CN is first
					if(info$nMin2[amplified.muts[a]]>max.CN1){
						max.CN2 = max.CN1
						max.CN1 = info$nMin2[amplified.muts[a]]
						frac2 = frac1
						frac1 = info$frac2[amplified.muts[a]]
					}else{
						max.CN2 = info$nMin2[amplified.muts[a]]
						frac2 = info$frac2[amplified.muts[a]]
					}
				}					
			}else{
				max.CN1 = info$nMaj1[amplified.muts[a]]
				frac1 = info$frac1[amplified.muts[a]]
				frac2=0
				if(!is.na(info$nMaj2[amplified.muts[a]])){
					#swap subclones, so that the one with the higher CN is first
					if(info$nMaj2[amplified.muts[a]]>max.CN1){
						max.CN2 = max.CN1
						max.CN1 = info$nMaj2[amplified.muts[a]]
						frac2 = frac1
						frac1 = info$frac2[amplified.muts[a]]						
					}else{
						max.CN2 = info$nMaj2[amplified.muts[a]]
						frac2 = info$frac2[amplified.muts[a]]
					}
				}	
			}
			best.err = info$mutation.copy.number[amplified.muts[a]] - 1
			best.CN=1
			for(j in 1:max.CN1){
				for(k in (j-1):min(j,max.CN2)){
					potential.CN = j * frac1 + k * frac2
					err = abs(info$mutation.copy.number[amplified.muts[a]]/potential.CN-1)
					if(err<best.err){
						info$no.chrs.bearing.mut[amplified.muts[a]] = potential.CN
						best.err=err
						best.CN = potential.CN
					}
				}
			}
			info$subclonal.fraction[amplified.muts[a]] = info$mutation.copy.number[amplified.muts[a]] / best.CN
		}
	}

  ##########################################################################

	#test for subclonal mutations
	
	#test whether mut burden is less than expected value for MCN = 1
	p.vals1 = sapply(1:length(non.zero.indices),function(v,e,i){
    prop.test(v$mut.count[i],v$mut.count[i] + v$WT.count[i],alternative="less")$p.value
    },v=info[non.zero.indices,], e= expected.burden.for.MCN[non.zero.indices])
	#test whether mut burden is above error rate (assumed to be 1 in 200)
	p.vals2 = sapply(1:length(non.zero.indices),function(v,i){
    prop.test(v$mut.count[i],v$mut.count[i] + v$WT.count[i],0.005,alternative="greater")$p.value
    },v=info[non.zero.indices,])

	subclonal.muts = non.zero.indices[p.vals1<=0.05 & p.vals2<=0.05]

	# use subclonal CN that minimises the difference in subclonal fraction from 1
	if(length(subclonal.muts)>0){
		for(a in 1:length(subclonal.muts)){
			#if there are no subclonal CNVs, don't adjust subclonal fraction
			if(is.na(info$frac2[subclonal.muts[a]])){next}
			#assume subclonal muts are on one chromosome copy, therefore mutation copy number must be subclonal fraction of the higher CN subclone (i.e. lost in lower CN subclone) or 1 (i.e. present in both subclones)
			if(info$nMaj1[subclonal.muts[a]]+info$nMin1[subclonal.muts[a]] > info$nMaj2[subclonal.muts[a]]+info$nMin2[subclonal.muts[a]]){	
				possible.subclonal.fractions = c(info$frac1[subclonal.muts[a]],1)
			}else{
				possible.subclonal.fractions = c(info$frac2[subclonal.muts[a]],1)
			}
			best.CN = possible.subclonal.fractions[which.min(abs(info$mutation.copy.number[subclonal.muts[a]]/possible.subclonal.fractions - 1))]
			#extra test 200313 - check whether subclonal CN results in clonal mutation, otherwise subclonal CN doesn't explain subclonal MCN
			if(best.CN != 1 & prop.test(info$mut.count[subclonal.muts[a]],info$mut.count[subclonal.muts[a]]+info$WT.count[subclonal.muts[a]],expected.burden.for.MCN[subclonal.muts[a]] * best.CN)$p.value > 0.05){
				info$subclonal.fraction[subclonal.muts[a]] = info$mutation.copy.number[subclonal.muts[a]] / best.CN
				info$no.chrs.bearing.mut[subclonal.muts[a]] = best.CN
			}
		}
	}	

	possible.zero.muts = intersect((1:length(info))[-non.zero.indices],which(!is.na(info$nMin1)))
	possible.zero.muts = c(possible.zero.muts,non.zero.indices[p.vals2>0.05])
	if(length(possible.zero.muts)>0){
		del.indices = which(info$nMin1[possible.zero.muts]==0 & !info$phase[possible.zero.muts]=="MUT_ON_RETAINED")
		info$subclonal.fraction[possible.zero.muts[del.indices]] = NA
		info$no.chrs.bearing.mut[possible.zero.muts[del.indices]] = 0
	}

  # convert GenomicRanges object to df
#   df = data.frame(chr=seqnames(info),
#                  starts=start(info)-1,
#                  ends=end(info))
  df = data.frame(chr=seqnames(info),
                pos=start(info))
  df = cbind(df, elementMetadata(info))

  out.file = paste(samplename,"_allDirichletProcessInfo.txt",sep="")
	if(!is.null(out.dir)){
		out.file = paste(out.dir,'/',out.file,sep="")
	}
	write.table(df, out.file, sep="\t", row.names=F, quote=F)
}

GetCellularity <- function(rho_and_psi_file) {
  d = read.table(rho_and_psi_file, header=T, stringsAsFactors=F)
  return(d['FRAC_GENOME','rho'])
}

GetWTandMutCount <- function(loci_file, allele_frequencies_file) {
  subs.data = read.table(loci_file, sep='\t', header=F, stringsAsFactors=F)
  subs.data.gr = GRanges(subs.data[,1], IRanges(subs.data[,2], subs.data[,2]), rep('*', nrow(subs.data)))
  elementMetadata(subs.data.gr) = subs.data[,c(3,4)]
  
  alleleFrequencies = read.table(allele_frequencies_file, sep='\t',header=F, stringsAsFactors=F)
  alleleFrequencies.gr = GRanges(alleleFrequencies[,1], IRanges(alleleFrequencies[,2], alleleFrequencies[,2]), rep('*', nrow(alleleFrequencies)))
  elementMetadata(alleleFrequencies.gr) = alleleFrequencies[,3:7]
  
  ref.indices = match(subs.data[,3],nucleotides)
  alt.indices = match(subs.data[,4],nucleotides)
  WT.count = as.numeric(sapply(1:nrow(alleleFrequencies),function(v,a,i){v[i,a[i]+2]},v=alleleFrequencies,a=ref.indices))
  mut.count = as.numeric(sapply(1:nrow(alleleFrequencies),function(v,a,i){v[i,a[i]+2]},v=alleleFrequencies,a=alt.indices))
  
  combined = data.frame(chr=subs.data[,1],pos=subs.data[,2],WTCount=WT.count, mutCount=mut.count)
  colnames(combined) = c("chr","pos","WT.count","mut.count")

#   combined.gr = GRanges(subs.data[,1], IRanges(subs.data[,2], subs.data[,2]+1), rep('*', nrow(subs.data)))
  combined.gr = GRanges(seqnames(subs.data.gr), ranges(subs.data.gr), rep('*', nrow(subs.data)))
  elementMetadata(combined.gr) = data.frame(WT.count=WT.count, mut.count=mut.count)
  return(combined.gr)
}

##############################################
# GetDirichletProcessInfo
##############################################
nucleotides = c("A","C","G","T")
info_counts = GetWTandMutCount(loci_file, allele_frequencies_file)
cellularity = GetCellularity(cellularity_file)
GetDirichletProcessInfo(samplename, cellularity, info_counts, subclone_file, out.dir=outdir, is.male=isMale)
