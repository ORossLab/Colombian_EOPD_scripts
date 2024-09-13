#!/bin/bash

## Tandem Repeat analysis
## Author: Gavrielatos Marios
## Mayo Clinic, Fl
## Dr. Owen Ross Lab
## Last revised 09/13/2024

#SBATCH --mail-user=surname.name@mail.edu           # CHANGE THIS
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50
#SBATCH --job-name=tr_truvari
#SBATCH --partition=cpu-short
#SBATCH --nodes=1
#SBATCH --ntasks=100


#SBATCH --time=80:00:00
#SBATCH --mem=100G
#SBATCH --chdir /path/to/output/directory/tr_analysis/   # CHANGE THIS
#SBATCH --output logs/%x.%N.%j.stdout
#SBATCH --error logs/%x.%N.%j.stderr


########################## COMMENTS ##########################
# CHANGE MAIL USER & CHDIR                                   #
# Change directory to output dir: ../tr_analysis             #
##############################################################

# Activate packages 
module load conda
source activate /path/to/mambaforge/envs/truvari

out_dir='/path/to/output/directory/tr_analysis'

out_dir_grch38=${out_dir}/grch38

if [ -d "$out_dir_grch38" ]; then rm -Rf $out_dir_grch38; fi

mkdir -p "${out_dir_chm13}"
mkdir -p "${out_dir_grch38}"


# Directory with your reference data - CHANGE THIS
grch38='/path/to/Reference/data/director/Ref_genome_grch38.fasta'



vcf_path="/path/to/output/directory/tr_analysis"
vcf_path_pathogenic="*.pathogenic.grch38.sorted.vcf.gz"
vcf_path_all="*.grch38.sorted.vcf.gz"



find $vcf_path -type f -name "${vcf_path_all}" > ${out_dir}/tr.list
find $vcf_path -type f -name "${vcf_path_pathogenic}" > ${out_dir}/pathogenic_tr.list

sort tr.list > tr.sorted.list 
sort pathogenic_tr.list > pathogenic_tr.sorted.list 

bcftools merge -m none --force-samples --file-list tr.sorted.list  | bgzip > tr.vcf.gz

bcftools view -c 1 tr.vcf.gz  -O z -o  tr_new.vcf.gz
tabix -p vcf tr.vcf.gz

truvari collapse -i tr.vcf.gz -o tr.merge.vcf -c tr.collapsed.vcf -f $grch38 \
-k common --sizemax 3000000

truvari vcf2df -i -f tr.vcf.gz tr.merge.jl
