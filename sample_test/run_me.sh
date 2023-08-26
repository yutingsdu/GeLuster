#!/bin/bash

######################## Set the environment variables ########################
#
# export PATH=/path_to_minimap2/:$PATH
# export PATH=/path_to_samtools/:$PATH
#
# e.g.
# export PATH=/home/yuting/yuting/Software/minimap2-2.24/:$PATH
# export PATH=/home/yuting/yuting/Software/samtools-1.10/:$PATH
#

../GeLuster -r Test.fastq -f fq -s cDNA -o geluster_outdir
