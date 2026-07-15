#!/bin/bash
#
# Example shell script to run the DPClust triplet pipeline stage 1
#
source ./triplet_pipeline_config.sh

R --vanilla --slave -q -f 01_run_core.R --args \
  -s ${samplename} -d ${datpath} -o ${outputdir} -i ${input} \
  --no.iters ${no_iters} --seed ${seed}
