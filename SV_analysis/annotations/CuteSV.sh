#!/bin/bash

## CuteSV annotation
## Author: Gavrielatos Marios
## Mayo Clinic, Fl
## Dr. Owen Ross Lab
## Last revised 09/13/2024

#SBATCH --mail-user=surname.name@mail.edu           # CHANGE THIS
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50
#SBATCH --job-name=cutesv
#SBATCH --partition=cpu-short
#SBATCH --nodes=1
#SBATCH --ntasks=100

#SBATCH --time=50:00:00
#SBATCH --mem=80G
#SBATCH --chdir /path/to/output/directory/cutesv           # CHANGE THIS
#SBATCH --output logs/%x.%N.%j.stdout
#SBATCH --error logs/%x.%N.%j.stderr

########################## COMMENTS ##########################
# CHANGE MAIL USER & CHDIR                                   #
# Change directory to output dir: ../cutesv                  #
##############################################################

# Activate packages 
module load conda
source activate /path/to/mambaforge/envs/cutesv

# Directory with your reference data - CHANGE THIS
grch38='/path/to/Reference/data/director/Ref_genome_grch38.fasta'
ref="grch38"


# Directory you will write your annotations - CHANGE THIS
out_dir="/path/to/output/dir/cutesv"
# Directory you will get your alignmnet bam files from - CHANGE THIS
base_dir="/path/to/alignments/dir"


for bam in $(find $base_dir -name \*.bam); do  
filename="${bam##*/}"
ID="${filename%%.*}"

aligner=$(echo "$filename" | rev | cut -d'.' -f2 | rev)
fi


cuteSV \
--max_size -1 --threads 30 \
--max_cluster_bias_INS 1000 \
--diff_ratio_merging_INS 0.9 \
--max_cluster_bias_DEL 1000 \
--diff_ratio_merging_DEL 0.5 \
--genotype \
$bam $ref_path "${ID}.${ref}.${aligner}.cutesv.vcf" $out_dir


done


