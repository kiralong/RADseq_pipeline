# README - RADseq Pipeline
Running Restriction Site associated DNA data instructions

Latest update: 10/20/20

This doc is notes and instructions on the pipeline to process raw RADseq data. The provided reference scripts are the ones I actually used in my own computing cluster, but you will need to tweak them to run on your own data and in your own computing environment.

## Overall Pipeline Summary

raw reads -> `process_radtags` -> `BWA` -> `samtools` -> `gstacks` -> `populations`

## Pipeline Steps

### Step 1: Process Raw Reads

**Requirements:**

- `stacks` version `2.x`, or `2.53` to match exactly with reference code.
- `gcc` version `7.x`, or `7.2.0` to match exactly with reference code.
  - TODO QUESTION `gcc` must be on the `PATH` for stacks to run?

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

> Note: If there is no reference genome, you'll need to use the denovo pipeline, also available in `stacks`.

Download your reference genome. (*M. vitellinus* GenBank Accession Number `GCA_001715985.3`)

Use BWA mem (Li & Durbin 2009) to align your RAD reads.

Check if [BWA](http://bio-bwa.sourceforge.net/) is installed in computing cluster.

See the file [`bwa_alignment.sh`](bwa_alignment.sh) for an example of running `bwa`.
Note that this script also pipes the bwa alignments directly into `samtools` in this single script; this prevents having to save the sam files!

### Step 3 - Process Aligned Reads with Samtools

**Requirements:**

- `samtools` version `1.7`

Sort and process your aligned reads with `samtools`.

Use samtools (Li et al. 2009). May need to install in computing cluster.
(http://www.htslib.org/)

See the file [`bwa_alignment.sh`](bwa_alignment.sh) and look at the last part for `samtools` loop.

### Step 4 - RAD Assembly and Genotyping

**Requirements:**

- `stacks` (same version earlier)

Use the Stacks 2 `gstacks` module (Rochette et al. 2019) to assemble and genotype your RAD loci.

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

You can add many more flags for different file outputs, such as `--vcf` to get a vcf file.

See the file [run_populations.sh](run_populations.sh) for an example with a bunch of different flags.
