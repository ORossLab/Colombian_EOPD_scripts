#!/bin/bash

## PhenoSV annotation
## Author: Gavrielatos Marios
## Mayo Clinic, Fl
## Dr. Owen Ross Lab
## Last revised 09/13/2024

#SBATCH --mail-user=surname.name@mail.edu           # CHANGE THIS
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50
#SBATCH --job-name=PhenoSV2
#SBATCH --partition=cpu-short
#SBATCH --nodes=1
#SBATCH --ntasks=128

#SBATCH --time=80:00:00
#SBATCH --mem=100G
#SBATCH --chdir /path/to/output/directory/annotations/phenosv/           # CHANGE THIS
#SBATCH --output logs/%x.%N.%j.stdout
#SBATCH --error logs/%x.%N.%j.stderr

########################## COMMENTS ##########################
# CHANGE MAIL USER & CHDIR                                   #
# Change directory to output dir: ../phenosv                 #
##############################################################

# Activate packages 
module load conda
source activate /path/to/mambaforge/envs/phenosv

# Directory with your reference data - CHANGE THIS
base_dir="/path/to/output/directory/annotations/phenosv"


file="${base_dir}/phenoSV_alignment.bed"

filename=$(basename "$file")

# Output directory - CHANGE THIS
out_dir="/path/to/output/directory/annotations/phenosv"
# mkdir -p $out_dir
# cd $out_dir
# ln -s $file $filename

# find ./  -type f -name 'input_*' -delete
bash /path/to/phenosv.sh $out_dir/${filename} $out_dir 90 'HP:0002180'
echo "Process of $dir_name finished"