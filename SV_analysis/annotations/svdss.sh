#!/bin/bash

## SVDSS pipeline
## Author: Gavrielatos Marios
## Mayo Clinic, Fl
## Dr. Owen Ross Lab
## Last revised 09/13/2024

#SBATCH --mail-user=surname.name@mail.edu           # CHANGE THIS
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50
#SBATCH --job-name=svdss
#SBATCH --partition=cpu-short
#SBATCH --nodes=1
#SBATCH --ntasks=100

#SBATCH --time=80:00:00
#SBATCH --mem=100G
#SBATCH --chdir /path/to/output/directory/svdss           # CHANGE THIS
#SBATCH --output logs/%x.%N.%j.stdout
#SBATCH --error logs/%x.%N.%j.stderr


########################## COMMENTS ##########################
# CHANGE MAIL USER & CHDIR                                   #
# Change directory to output dir: ../svdss                   #
# If you want to get ALL SVs, without quality control,       #
# include the --no-qc parameter                              #
# If you only was to keep the SVs with filter=PASS, remove   #
# the --no-qc parameter                                      #
##############################################################

# Activate packages 
module load conda
source activate /path/to/mambaforge/envs/svdss 

# Directory with your reference data - CHANGE THIS
grch38='/path/to/Reference/data/director/Ref_genome_grch38.fasta'
grch38_fmd='/path/to/Reference/data/director/Ref_genome_grch38.fmd'

# Index Reference
SVDSS index --fasta $grch38 --index $grch38_fmd --threads 30

# Directory you will get your alignmnet bam files from - CHANGE THIS
base_dir="/path/to/alignments/dir/*bam"




for file in $base_dir; do
filename="${file##*/}"
ID="${filename%%.*}"

bam_grch38="/path/to/alignments/dir/${ID}.grch38.pbmm2.bam"

# GRCh38 reference
out_dir="/path/to/output/dir/svdss/${ID}"
rm -rf $out_dir
mkdir -p $out_dir

# Smooth sample
SVDSS smooth --threads 90 --reference ${grch38} --bam $bam_grch38 --workdir ${out_dir}
samtools index -@10 ${out_dir}/smoothed.selective.bam
# Extract SFS from smoothed BAM file
SVDSS search --threads 90 --assemble --index $grch38_fmd --bam ${out_dir}/smoothed.selective.bam --workdir ${out_dir}
# Genotype SVs from the assembled superstrings 
SVDSS call --threads 90 --reference $grch38 --bam ${out_dir}/smoothed.selective.bam --workdir ${out_dir}
mv ${out_dir}/svs_poa.vcf ${out_dir}/../${ID}.svdss.grch38.vcf
rm -rf $out_dir

done