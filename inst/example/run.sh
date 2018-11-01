#!/bin/bash
#
# Example shell script to run DPClust on a single sample
#
mkdir -p output
R --vanilla --slave -q -f dpclust_pipeline.R --args -r 1 -d ../extdata/simulated_data/Data/ -o output -i ../extdata/simulated_data/simulated_1d.txt
