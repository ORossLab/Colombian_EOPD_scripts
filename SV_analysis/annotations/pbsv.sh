#!/bin/bash

## Pbsv annotation
## Author: Gavrielatos Marios
## Mayo Clinic, Fl
## Dr. Owen Ross Lab
## Last revised 09/13/2024

#SBATCH --mail-user=surname.name@mail.edu           # CHANGE THIS
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50
#SBATCH --job-name=pbsv
#SBATCH --partition=cpu-short
#SBATCH --nodes=1
#SBATCH --ntasks=100

#SBATCH --time=50:00:00
#SBATCH --mem=80G
#SBATCH --chdir /path/to/output/directory/pbsv           # CHANGE THIS
#SBATCH --output logs/%x.%N.%j.stdout
#SBATCH --error logs/%x.%N.%j.stderr

########################## COMMENTS ##########################
# CHANGE MAIL USER & CHDIR                                   #
# Change directory to output dir: ../pbsv                    #
##############################################################

# Activate packages 
module load conda
source activate /path/to/mambaforge/envs/pb-human-wgs

# Directory with your reference data - CHANGE THIS
grch38='/path/to/Reference/data/director/Ref_genome_grch38.fasta'
ref="grch38"

# Directory you will write your annotations - CHANGE THIS
out_dir="/path/to/output/dir/pbsv"
# Directory you will get your alignmnet bam files from - CHANGE THIS
base_dir="/path/to/alignments/dir"

for bam in $(find $base_dir -name \*.bam); do  
filename="${bam##*/}"
ID="${filename%%.*}"

aligner=$(echo "$filename" | rev | cut -d'.' -f2 | rev)

pbsv discover --log-level INFO $bam "${out_dir}/${ID}.${ref}.${aligner}.pbsv.svsig.gz"
pbsv call --log-level INFO -j 90 --max-dup-length 20M --max-ins-length 1M $grch38 "${out_dir}/${ID}.${ref}.${aligner}.pbsv.svsig.gz" "${out_dir}/${ID}.${ref}.${aligner}.pbsv.vcf"

done
