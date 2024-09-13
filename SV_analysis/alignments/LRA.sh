#!/bin/bash

## LRA alignment
## Author: Gavrielatos Marios
## Mayo Clinic, Fl
## Dr. Owen Ross Lab
## Last revised 09/13/2024

#SBATCH --mail-user=surname.name@mail.edu           # CHANGE THIS
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50
#SBATCH --job-name=LRA
#SBATCH --partition=cpu-short
#SBATCH --nodes=1
#SBATCH --ntasks=128

#SBATCH --time=80:00:00
#SBATCH --mem=128G
#SBATCH --chdir /path/to/output/directory/bams
#SBATCH --output logs/%x.%N.%j.stdout
#SBATCH --error logs/%x.%N.%j.stderr

########################## COMMENTS ##########################
# CHANGE MAIL USER & CHDIR                                   #
# Change directory to output dir: ../LRA/bams                #
# If samtools causes issues, comment out the samtools lines  #
# and run them in the terminal after the end of the LRA      #
# alignments                                                 #
# At the end of the script remove the extra files from the   #
# middle steps:                                              #
# *.grch38.lra.unsorted.sam,                                 #   
# *.grch38.lra.unsorted.bam,                                 #
##############################################################

# Activate packages 
module load conda
source activate /path/to/mambaforge/envs/lra

# Directory with your reference data - CHANGE THIS
grch38='/path/to/Reference/data/director/Ref_genome_grch38.fasta'
lra index $grch38

# path to unaligned PacBio bam files - CHANGE THIS
ubams="/path/to/unaligned/bam/files/ubams/*"

for ubam in $ubams; do
filename_ubam=$(basename ${ubam})
ID="${filename_ubam%%.*}"

lra align -CCS $grch38 $ubam -t 128 -p s > ${ID}.grch38.lra.unsorted.sam
samtools view -hb -@10 ${ID}.grch38.lra.unsorted.sam > ${ID}.grch38.lra.unsorted.bam
samtools sort -@8 ${ID}.grch38.lra.unsorted.bam  -o ${ID}.grch38.lra.bam 
samtools index -@10 ${ID}.grch38.lra.bam 
done