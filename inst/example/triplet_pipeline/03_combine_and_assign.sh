#!/bin/bash
#
# Example shell script to run the DPClust triplet pipeline stage 3
# (combine density estimates across all triplets into a consensus).
# Requires 01_run_core.sh and 02_estimate_density.sh to have completed first.
#
source ./triplet_pipeline_config.sh

R --vanilla --slave -q -f 03_combine_and_assign.R --args \
  -s ${samplename} -d ${datpath} -o ${outputdir} -i ${input}
