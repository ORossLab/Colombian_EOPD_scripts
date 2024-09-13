#!/bin/bash

## Sniffles2 annotation
## Author: Gavrielatos Marios
## Mayo Clinic, Fl
## Dr. Owen Ross Lab
## Last revised 09/13/2024

#SBATCH --mail-user=surname.name@mail.edu           # CHANGE THIS
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50
#SBATCH --job-name=Sniffles
#SBATCH --partition=cpu-short
#SBATCH --nodes=1
#SBATCH --ntasks=100

#SBATCH --time=50:00:00
#SBATCH --mem=80G
#SBATCH --chdir /path/to/output/directory/sniffles           # CHANGE THIS
#SBATCH --output logs/%x.%N.%j.stdout
#SBATCH --error logs/%x.%N.%j.stderr

########################## COMMENTS ##########################
# CHANGE MAIL USER & CHDIR                                   #
# Change directory to output dir: ../sniffles                #
# If you want to get ALL SVs, without quality control,       #
# include the --no-qc parameter                              #
# If you only was to keep the SVs with filter=PASS, remove   #
# the --no-qc parameter                                      #
##############################################################

# Activate packages 
module load conda
source activate /path/to/mambaforge/envs/sniffles2


# Directory with your reference data - CHANGE THIS
grch38='/path/to/Reference/data/director/Ref_genome_grch38.fasta'
ref="grch38"

# Directory you will get your alignmnet bam files from - CHANGE THIS
base_dir="/path/to/alignments/dir"



for bam in $(find $base_dir -name \*.bam); do  
filename="${bam##*/}"
ID="${filename%%.*}"

aligner=$(echo "$filename" | rev | cut -d'.' -f2 | rev)

sniffles \
-i $bam \
--reference $ref_path \
--no-qc \
--threads 50 \
--output-rnames \
-v "${ID}.${ref}.${aligner}.sniffles2.vcf"

done


