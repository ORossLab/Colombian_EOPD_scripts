#!/bin/bash

## CNV analysis
## Author: Gavrielatos Marios
## Mayo Clinic, Fl
## Dr. Owen Ross Lab
## Last revised 09/13/2024

#SBATCH --mail-user=surname.name@mail.edu           # CHANGE THIS
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50
#SBATCH --job-name=hificnv
#SBATCH --partition=cpu-short
#SBATCH --nodes=1
#SBATCH --ntasks=100

#SBATCH --time=80:00:00sted
#SBATCH --mem=128G
# #SBATCH --mem-per-cpu=8G
#SBATCH --chdir /path/to/output/directory/cnv_analysis           # CHANGE THIS
#SBATCH --output logs/%x.%N.%j.stdout
#SBATCH --error logs/%x.%N.%j.stderr

########################## COMMENTS ##########################
# CHANGE MAIL USER & CHDIR                                   #
# Change directory to output dir: ../cnv_analysis            #
##############################################################

# Activate packages 
module load conda
source activate /path/to/mambaforge/envs/hificnv


# Directory with your reference data - CHANGE THIS
grch38='/path/to/Reference/data/director/Ref_genome_grch38.fasta'

# Replace here the file names of the male samples  - CHANGE THIS
males='sample_1 sample_2 sample_3'
# Excluded regions  - CHANGE THIS
exclude="/path/to/excluded_regions/cnv.excluded_regions.hg38.bed.gz"

# If you have multiple samples add the path to your folders and the extension of those files - CHANGE THIS
## Example: 
## My files are in the /my/folder/bams directory &
## I want to use the alignments on the grch38 ref genome with the pbmm2 alignment. These files have the extension ".grch38.pbmm2.bam"
## Then:
## bams="/my/folder/bams/*.grch38.pbmm2.bam"
bams="/my/folder/bams/*extension.bam"

for bam in $bams; do
filename_bam=$(basename ${bam})
filename="${filename_bam%%.*}"

if [[ " $males " =~ " $filename " ]]; then
    extended_cn="/path/to/expected_cn/male_expected_cn.hg38.bed"
else
    extended_cn="/path/to/expected_cn/female_expected_cn.hg38.bed"
fi

hificnv \
    --bam $bam \
    --ref $grch38 \
    --exclude $exclude \
    --expected-cn $extended_cn \
    --threads 10 \
    --output-prefix "${filename}.grch38"

echo ""
echo "###############################################################"
echo ""

done
