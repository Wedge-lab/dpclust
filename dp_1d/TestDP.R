setwd("C:/R_code")
source("interconvertMutationBurdens.R")
source("subclone_Dirichlet_Gibbs_sampler_binomial.R")
#get the number of mutant reads into mutReads and the total number of reads at each mutation into totalReads, then run the next line
mutReads <- c(rep(10,25), rep(50,30), rep(90,25))
#mutReads = sample(mutReads)
totalReads <- c(rep(100,80) )
GS.data<-subclone.dirichlet.gibbs(y=mutReads,N=totalReads)
Gibbs.subclone.density.est(GS.data,"DirichletProcessplot.png", post.burn.in.start = 300, post.burn.in.stop = 1000, y.max=50)

mutReads.binomial <- c(rbinom(25,100,0.1), rbinom(30,100,0.5), rbinom(25,100,0.9))
GS.data.binomial<-subclone.dirichlet.gibbs(y=mutReads.binomial,N=totalReads)
Gibbs.subclone.density.est(GS.data.binomial,"DirichletProcessplotBinomial.png", post.burn.in.start = 300, post.burn.in.stop = 1000, y.max=50)
