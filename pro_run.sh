#!/bin/bash
export PATH=/home/yuting/yuting/ThirdG-transcriptome/NewAssembler_v0/src/plugin:$PATH

echo "First stage begin..."
date

Nx=N60

query=../../Notts_Run2_20171108_1D.pass.dedup.fastq

PseudoRef $query fq $Nx #get seudo reference and seudoref-reads.info

ref=$Nx".seudo.fasta"
prefix=$Nx".seudo"
time minimap2 -ax splice -k 11 $ref $query > $prefix.sam # mapping reads to the reference 
samtools sort  -o $prefix.bam $prefix.sam 
get_alignment_info $prefix".bam" $Nx".alignment.info" #from stringtie


query2=$Nx".seudo.fasta"
prefix=$Nx".seudo-self"

time minimap2 -ax splice -k 11 -X -N 1 $ref $query2 > $prefix.sam
samtools sort  -o $prefix.bam $prefix.sam
get_alignment_info $prefix".bam" $Nx".self.alignment.info" #from stringtie


get_cluster $query fq $Nx".seudoref-reads.info" $Nx".alignment.info" $Nx".self.alignment.info" $Nx

date
echo "First stage finish!"


#./pro_run2.sh

