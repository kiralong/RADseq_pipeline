#!/bin/bash
#PBS -q aces
#PBS -l walltime=168:00:00,nodes=1:ppn=8
#PBS -N manacus_bwa_alignment
#PBS -j oe

module load gcc/7.2.0

bwa_path=/projects/aces/kira2/programs/bwa-0.7.17
samtools_path=/projects/aces/kira2/programs/samtools-1.7
sample_names_path=/projects/aces/kira2/sample_info/barcodes_jan_2018.tsv
reads_path=/projects/aces/kira2/processed_samples
out_path=/projects/aces/kira2/aligned_samples/ASM171598v3
database_path=/projects/aces/kira2/reference_genome/ASM171598v3/GCF_001715985.3_ASM171598v3_genomic.fna.gz
parallel_path=/projects/aces/kira2/programs/parallel-20180322/bin

cat $sample_names_path | cut -f 2 |
while read sample; do
	echo "$bwa_path/bwa mem -t 2 $database_path $reads_path/${sample}.1.fq.gz $reads_path/${sample}.2.fq.gz | \
		$samtools_path/samtools view -b -h | $samtools_path/samtools sort --threads 2 -o $out_path/${sample}.bam";
done | $parallel_path/parallel -j 4
