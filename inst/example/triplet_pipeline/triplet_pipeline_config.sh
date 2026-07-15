#
# Shared configuration for the triplet pipeline example scripts.
#
# no_iters must be identical between steps 1 and 2 - step 2 slices the MCMC
# output saved by step 1 using this value (burn-in is derived from it the
# same way in both steps, so does not need to be repeated here). This is why
# it is passed explicitly by both steps' .sh scripts, unlike other tunable
# parameters (e.g. max.considered.clusters) which are left at the R scripts'
# own defaults since they carry no such cross-step dependency. conc_param,
# cluster_conc and density.smooth are internal tuning parameters not meant
# to be changed by regular users, so are not exposed here at all.
#
# Source this file from each step's .sh script rather than repeating these
# values, so there is exactly one place to change them.
#
no_iters=1000
seed=123

# How many triplets step 2 (density estimation) runs at once. Each triplet
# process can use several GB of memory, so keep this conservative unless 
# you know your machine has the headroom.
n_parallel=1

samplename=simulated_002
datpath=../../extdata/simulated_data/Data/
outputdir=output
input=../../extdata/simulated_data/simulated_triplet.txt
