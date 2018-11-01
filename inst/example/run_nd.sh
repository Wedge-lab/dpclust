#!/bin/bash
#
# Example shell script to run DPClust on a multi-sample case
#
mkdir -p output
R --vanilla --slave -q -f dpclust_pipeline.R --args -r 1 -d simulated_data/Data/ -o output -i simulated_data/simulated.txt
