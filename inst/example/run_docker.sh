#!/bin/bash
#
# Example shell script to run the Dockerised DPClust on a single sample
#
R --vanilla --slave -q -f /opt/dpclust/inst/example/dpclust_pipeline.R --args -r 1 -d /opt/dpclust/inst/extdata/simulated_data/Data/ -o /mnt/output -i /opt/dpclust/inst/extdata/simulated_data/simulated_1d.txt
