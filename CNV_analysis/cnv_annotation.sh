#!/bin/bash

## CNV annotation
## Author: Gavrielatos Marios
## Mayo Clinic, Fl
## Dr. Owen Ross Lab
## Last revised 09/13/2024

#SBATCH --mail-user=surname.name@mail.edu           # CHANGE THIS
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50
#SBATCH --job-name=cnv_phenoSV2
#SBATCH --partition=cpu-short
#SBATCH --nodes=1
#SBATCH --ntasks=128

#SBATCH --time=80:00:00
#SBATCH --mem=100G
#SBATCH --chdir /path/to/output/directory/cnv_analysis/phenosv           # CHANGE THIS
#SBATCH --output logs/%x.%N.%j.stdout
#SBATCH --error logs/%x.%N.%j.stderr

########################## COMMENTS ##########################
# CHANGE MAIL USER & CHDIR                                   #
# Change directory to output dir: ../cnv_analysis/phenosv    #
##############################################################

# Activate packages 
module load conda
source activate /path/to/mambaforge/envs/phenosv


base_dir="/path/to/output/directory/cnv_analysis/phenosv"


file="${base_dir}/phenoSV_CNV_alignment.bed"

filename=$(basename "$file")


out_dir="/path/to/output/directory/cnv_analysis/phenosv"
# mkdir -p $out_dir
# cd $out_dir
# ln -s $file $filename

# find ./  -type f -name 'input_*' -delete
# HPO term: neurodegeneration
bash /home/mayo/m301955/my_data/tools/PhenoSV/phenosv/model/phenosv.sh $file $out_dir 500 'HP:0002180'
