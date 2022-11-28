#!/bin/bash
#
# Example shell script to run DPClust on a single sample
#
mkdir -p output
R --vanilla --slave -q -f Pipeline_phasingversion.R --args -r 1 -d ./05_Data -o ./05_PhaDPClust -i ./${samplename}_phadp.txt -a phasing_ass
