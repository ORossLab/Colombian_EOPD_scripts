#!/bin/bash

## Pbmm2 alignment
## Author: Gavrielatos Marios
## Mayo Clinic, Fl
## Dr. Owen Ross Lab
## Last revised 09/13/2024

#SBATCH --mail-user=surname.name@mail.edu           # CHANGE THIS
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50
#SBATCH --job-name=pbmm2_hg38
#SBATCH --partition=cpu-short
#SBATCH --nodes=1
#SBATCH --ntasks=128

#SBATCH --time=80:00:00
#SBATCH --mem=128G
#SBATCH --chdir /path/to/output/directory           # CHANGE THIS
#SBATCH --output logs/%x.%N.%j.stdout
#SBATCH --error logs/%x.%N.%j.stderr

########################## COMMENTS ##########################
# CHANGE MAIL USER & CHDIR                                   #
# Change directory to output dir: ../pbmm2/bams              #
##############################################################

# Activate packages 
module load conda
source activate /path/to/mambaforge/envs/pb-human-wgs

# Directory with your reference data - CHANGE THIS
grch38='/path/to/Reference/data/director/Ref_genome_grch38.fasta'
pbmm2 index $grch38 Ref_genome_grch38.mmi 
mmi='/path/to/Reference/data/director/Ref_genome_grch38.mmi'

# path to unaligned PacBio bam files - CHANGE THIS
ubams="/path/to/unaligned/bam/files/ubams/*"

for ubam in $ubams; do
filename_bam=$(basename ${ubam})
filename="${filename_bam%.*}"

# Run pbmm2
pbmm2 align -j 128 --sort --preset HIFI $mmi $ubam "${filename}.grch38.mapped.sorted.bam" 

done