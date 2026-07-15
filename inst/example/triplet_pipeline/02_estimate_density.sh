#!/bin/bash
#
# Example shell script to run the DPClust triplet pipeline stage 2
# (density estimation) across every possible triplet of subsamples.
# Requires 01_run_core.sh to have been run first.
#
source ./triplet_pipeline_config.sh

echo "Starting up: $(date)"

# Density estimation is done per triplet, so set choose_number to 3
choose_number=3

# Number of subsamples for this donor, straight from the samplesheet
choose_from=$(Rscript -e "sp <- read.table('${input}', header=TRUE, stringsAsFactors=FALSE); cat(sum(sp\$sample=='${samplename}'))")

# Ask the package how many triplets that implies
total_triplets=$(Rscript -e "suppressMessages(library(DPClust)); cat(determine_num_triplets(${choose_from}))")

run_one () {
  choose_index=$1
  echo "Triplet ${choose_index} / ${total_triplets}"
  R --vanilla --slave -q -f 02_estimate_density.R --args \
    -s ${samplename} -d ${datpath} -o ${outputdir} -i ${input} \
    --choose.from ${choose_from} --choose.number ${choose_number} --choose.index ${choose_index} \
    --no.iters ${no_iters}
  echo "Completed triplet ${choose_index} / ${total_triplets}: $(date)"
}
export -f run_one
export samplename datpath outputdir input choose_from choose_number total_triplets no_iters

# Now estimate the density for every combination, n_parallel at a time
seq 1 "${total_triplets}" | xargs -n 1 -P "${n_parallel}" -I{} bash -c 'run_one {}'
