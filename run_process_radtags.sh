#!/bin/bash
#PBS -q aces
#PBS -l walltime=72:00:00,nodes=1:ppn=4
#PBS -N manacus_process_radtags
#PBS -k oe

module load gcc/7.2.0

stacks=/projects/aces/kira2/programs/stacks-2.0Beta9
raw=/projects/aces/kira2/raw_data
processed=/projects/aces/kira2/processed_samples
barcodes=/projects/aces/kira2/sample_info/barcodes_jan_2018.txt

$stacks/process_radtags \
  -i gzfastq \
  -p $raw \
  -P \
  -o $processed \
  -b $barcodes \
  -c \
  -q \
  -r \
  -e sbfI
