#!/bin/bash

## Tandem Repeat analysis
## Author: Gavrielatos Marios
## Mayo Clinic, Fl
## Dr. Owen Ross Lab
## Last revised 09/13/2024

#SBATCH --mail-user=surname.name@mail.edu           # CHANGE THIS
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50
#SBATCH --job-name=trgt
#SBATCH --partition=cpu-short
#SBATCH --nodes=1
#SBATCH --ntasks=100

#SBATCH --time=80:00:00
#SBATCH --mem=128G
# #SBATCH --mem-per-cpu=8G
#SBATCH --chdir /path/to/output/directory/tr_analysis   # CHANGE THIS
#SBATCH --output logs/%x.%N.%j.stdout
#SBATCH --error logs/%x.%N.%j.stderr

########################## COMMENTS ##########################
# CHANGE MAIL USER & CHDIR                                   #
# Change directory to output dir: ../tr_analysis             #
##############################################################

# Activate packages 
module load conda
source activate /path/to/mambaforge/envs/trgt


# Directory with your reference data - CHANGE THIS
grch38='/path/to/Reference/data/director/Ref_genome_grch38.fasta'

# Bed file with polymorphic repeats - Download from github - CHANGE THIS
repeats='/path/to/trgt/repeats/polymorphic_repeats.hg38.bed'

# If you have multiple samples add the path to your folders and the extension of those files - CHANGE THIS
## Example: 
## My files are in the /my/folder/bams directory &
## I want to use the alignments on the grch38 ref genome with the pbmm2 alignment. These files have the extension ".grch38.pbmm2.bam"
## Then:
## bams="/my/folder/bams/*.grch38.pbmm2.bam"
bams="/my/folder/bams/*extension.bam"

for bam in $bams; do
filename_bam=$(basename ${bam})
ID="${filename_bam%%.*}"



/home/mayo/m301955/my_data/tools/trgt/trgt_v0.8.0 --genome $grch38 --reads $bam --repeats $repeats --output-prefix ${ID}.grch38 --threads 10

bcftools sort -Ob -o ${ID}.grch38.sorted.vcf.gz ${ID}.grch38.vcf.gz
bcftools index ${ID}.grch38.sorted.vcf.gz

samtools sort -@10 -o ${ID}.grch38.spanning.sorted.bam ${ID}.grch38.spanning.bam
samtools index -@10 ${ID}.grch38.spanning.sorted.bam

 done

