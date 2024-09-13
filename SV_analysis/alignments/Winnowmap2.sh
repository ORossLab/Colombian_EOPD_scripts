#!/bin/bash

## Winnowmap2 alignment
## Author: Gavrielatos Marios
## Mayo Clinic, Fl
## Dr. Owen Ross Lab
## Last revised 09/13/2024

#SBATCH --mail-user=surname.name@mail.edu           # CHANGE THIS
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50
#SBATCH --job-name=Winnowmap
#SBATCH --partition=cpu-short
#SBATCH --nodes=1
#SBATCH --ntasks=100

#SBATCH --time=80:00:00
#SBATCH --mem=128G
#SBATCH --chdir /path/to/output/directory/bams
#SBATCH --output logs/%x.%N.%j.stdout
#SBATCH --error logs/%x.%N.%j.stderr

########################## COMMENTS ##########################
# CHANGE MAIL USER & CHDIR                                   #
# Change directory to output dir: ../winnowmap/bams          #
# At the end of the script remove the extra files from the   #
# middle steps:                                              #
# *.grch38.winnowmap.unsorted.sam,                           #
# *${ID}.grch38.winnowmap.unsorted.bam                       #
##############################################################

# Activate packages 
module load conda
source activate /path/to/mambaforge/envs/Winnowmap

# Directory with your reference data - CHANGE THIS
grch38='/path/to/Reference/data/director/Ref_genome_grch38.fasta'

meryl count threads=8 k=15 output meryl_GRCh38_DB $grch38
meryl print greater-than distinct=0.9998 meryl_GRCh38_DB > repetitive_k15_GRCh38.txt

# Winnowmap2 takes as input fastq.gz files
# You can create these using the bam2fasatq.slrm.sh file provided
fastqs="/path/to/fastq/*fastq.gz"

for fastq in $fastqs; do
filename_fastq=$(basename ${fastq})
ID="${filename_fastq%%.*}"

winnowmap -t 100 -W repetitive_k15_GRCh38.txt -ax map-pb -Y -L --eqx --cs $grch38 $fastq  > ${ID}.grch38.winnowmap.unsorted.sam
samtools view -hb -@10 ${ID}.grch38.winnowmap.unsorted.sam > ${ID}.grch38.winnowmap.unsorted.bam

echo="Start ${ID} GRCh38 shortment."
samtools sort -@8 ${ID}.grch38.winnowmap.unsorted.bam -o ${ID}.grch38.winnowmap.bam

samtools index -@50 ${ID}.grch38.winnowmap.bam


done



