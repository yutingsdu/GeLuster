#!/bin/bash

######################## Set the environment variables ########################
#
#export PATH=/path_to_minimap2/:$PATH
#export PATH=/path_to_samtools/:$PATH

../GeLuster -r Test.fastq -f fq -s cDNA -o geluster_outdir
