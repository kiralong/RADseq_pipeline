# README - RADseq Pipeline
Running Restriction Site associated DNA data instructions

Latest update: 10/27/20

This doc is notes and instructions on a pipeline to process raw RADseq data. The provided reference scripts are the ones I actually used in my own computing cluster, but you will need to tweak them to run on your own data and in your own computing environment.

## Overall Pipeline Summary

raw reads -> `process_radtags` -> `BWA` -> `samtools` -> `gstacks` -> `populations`

## Pipeline Steps

### Step 1: Process Raw Reads

**Requirements:**

- `stacks` version `2.x`, or `2.53` to match exactly with reference code.
- `gcc` version `7.x`, or `7.2.0` to match exactly with reference code.

Take your raw reads (should be a `.fastq` file or `.fastq.gz` if zipped) from the sequencing facility and filter, demultiplex, and trim your raw reads using the Stacks 2 software module `process_radtags`.

See [Stacks - Process RADtags](https://catchenlab.life.illinois.edu/stacks/comp/process_radtags.php) to customize the below example code.
Change the example code to your own paths and select which enzyme you used, adapters, etc.

```bash
#!/bin/bash
/path/to/stacks/process_radtags \
    -p /path/to/raw/fastq \
    -o /path/to/output/directory \
    -b /path/to/barcodes/file \ #You need a file that lists all barcodes used
    --paired \                  #Paired or unpaired data
    --clean \                   #Tells program to clean and remove uncalled bases
    --quality \                 #Discards reads will low quality scores
    --rescue \                  #Rescues barcodes and RADtags
    --renz-1 sbfI \             #Restriction enzyme used, if double digest, add the flag renz_2
    -i gzfastq \                #input file type, default gzfastq if unknown
    --adapter_1 GATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCAGAACAA \  # P2 top from SB's RADseq protocol, seen in read 1
    --adapter_2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \               # P1 bottom from SB's RADseq protocol, seen in read 2 rev comp
    --adapter_mm 2              #Number of mismatches allowed in the adapter sequence
```

For an example of the barcodes file, see [`barcodes_jan_2018.txt`](barcodes_jan_2018.txt) in this repository.
This file should be 2 columns. The first column is the 7 bp barcode and the 2nd column is the sample ID.

For the script to run `process_radtags` in a computing cluster with the shorthand flags, see [`run_process_radtags.sh`](run_process_radtags.sh). This was used for the first *Manacus vitellinus* run.

### Step 2 - Align RAD Data to Reference Genome

**Requirements:**

- `bwa` version `0.7.x`, or `0.7.17` to match exactly reference code.
- A reference genome, usually as a fasta.gz file

> Note: If there is no reference genome, you'll need to use the denovo pipeline, also available in `stacks`.

Download your reference genome. (For my species I used: *M. vitellinus* GenBank Accession Number `GCA_001715985.3`)

Use [BWA](http://bio-bwa.sourceforge.net/) mem (Li & Durbin 2009) to align your RAD reads. You'll first need to make a BWA reference index. 

```bash
#Create a BWA reference index
bwa index /path/to/reference.fasta.gz

#By default, the above command will use the name of your input FASTA file as the basename for your bwa index. You will feed this index into your next script to run bwa. If you want to change the name, use the '-p' flag such as the following: bwa index /path/to/reference.fasta.gz -p my_bwa_database.
```
Now that you have made your index (referred to as your 'database' in the following example), run bwa. You can follow the below sample code:

```bash
cat $sample_names_path | cut -f 2 |
  while read sample
  do
    $bwa_path/bwa mem -t 2 $database_path $reads_path/${sample}.1.fq.gz $reads_path/${sample}2.fq.gz | \
      samtools_path/samtools view -b -h | $samtools_path/samtools sort --threads 2 -o
    $outpath/${sample}.bam
```

See the file [`bwa_alignment.sh`](bwa_alignment.sh) for an example of running `bwa`.
Note that this script and the example code pipes the bwa alignments directly into `samtools`; this prevents having to save the sam files! This script reads the sample names from the barcodes file from step 1: `process_radtags` by cutting just the samples names in column 2 and loops over each of the samples one at a time in the while loop. Every time it  generates the `bwa+samtools` command for that sample and runs them in parallel. At the end of this after using both `bwa` and `samtools` you should have a directory full of .bam files for each of your samples. 

### Step 3 - Process Aligned Reads with Samtools

**Requirements:**

- `samtools` version `1.7`

Sort and process your aligned reads with `samtools`.

Use samtools (Li et al. 2009). May need to install in computing cluster.
(http://www.htslib.org/)

See the file [`bwa_alignment.sh`](bwa_alignment.sh) and look at the last part for adding `samtools`. Note that in this pipeline guide this step was integrated with the previous step 2 with `bwa`.

### Step 4 - RAD Assembly and Genotyping

**Requirements:**

- `stacks` (same version earlier)
- Your sorted bam files
- A ["population map"](https://catchenlab.life.illinois.edu/stacks/manual/#popmap)

Use the Stacks 2 [`gstacks`](https://catchenlab.life.illinois.edu/stacks/comp/gstacks.php) module (Rochette et al. 2019) to assemble and genotype your RAD loci.

Example code:
```bash
#!/bin/bash
/path/to/stacks/gstacks \
  -I /path/to/directory/with/bam/files \
  -O /path/to/output/directory \
  -M /path/to/popmap \                    #You will need a text file with a list of all the samples and the population they are in as two columns seperated by tabs. 
  -t 8                                    #This is the number of threads you want for parallelizing. The default without this flag is 1

#Note that this is the step where you will probably want to remove pcr duplicates if you have paired-end, single-digest RAD data. Then add the flag --rm-pcr-duplicates to the above code. 
```

### Step 5 - Filter Genotypes and Calculate Genome Statistics

Use the Stacks 2 module `populations` to filter your genotypes by removing loci with high missing data, loci that aren't present in all populations, or with low minimum allele frequency.
Populations is also where you can get the sumstats file to look at fst, pi etc.

Example code:
```bash
#!/bin/bash
/path/to/stacks/populations
    --in-path /path/to/gstacks/output
    --out-path /path/to/output/directory
    --popmap /path/to/population/map
    --threads 12                          #for parallelizing in the cluster
    --min-mac $mac                        #minimum allele count, I usually use 3
    --min-population $p                   #Number of populations the allele must be present in
    --min-samples-per-pop $r              #I usually like r .50 - .80
```

You can add many more flags for different file outputs, such as `--vcf` to get a vcf file or '--hzar' for a HZAR file.

See the file [run_populations.sh](run_populations.sh) for an example with a bunch of different flags.
