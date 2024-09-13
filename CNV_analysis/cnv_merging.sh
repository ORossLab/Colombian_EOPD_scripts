#!/bin/bash

## CNV analysis
## Author: Gavrielatos Marios
## Mayo Clinic, Fl
## Dr. Owen Ross Lab
## Last revised 09/13/2024

#SBATCH --mail-user=surname.name@mail.edu           # CHANGE THIS
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50
#SBATCH --job-name=cnv_truvari
#SBATCH --partition=cpu-short
#SBATCH --nodes=1
#SBATCH --ntasks=100

#SBATCH --time=80:00:00
#SBATCH --mem=100G
#SBATCH --chdir /path/to/output/directory/cnv_analysis/truvari           # CHANGE THIS
#SBATCH --output logs/%x.%N.%j.stdout
#SBATCH --error logs/%x.%N.%j.stderr


########################## COMMENTS ##########################
# CHANGE MAIL USER & CHDIR                                   #
# Remember to create the CHDIR directory before running      #
# Change directory to output dir: ../cnv_analysis/truvari    #
##############################################################

# Activate packages 
module load conda
source activate /path/to/mambaforge/envs/truvari

out_dir='/path/to/output/directory/cnv_analysis/truvari'

out_dir_grch38=${out_dir}/grch38

if [ -d "$out_dir_grch38" ]; then rm -Rf $out_dir_grch38; fi

mkdir -p "${out_dir_grch38}"

# Directory with your reference data - CHANGE THIS
grch38='/path/to/Reference/data/director/Ref_genome_grch38.fasta'

vcf_path="/path/to/output/directory/cnv_analysis/"
vcf_path_chm13="*chm13*.vcf.gz"
vcf_path_grch38="*grch38*.vcf.gz"



find $vcf_path -type f -name "*grch38*.vcf.gz" > ${out_dir}/grch38_tmp.list
find $vcf_path -type f -name "*chm13*.vcf.gz" > ${out_dir}/chm13_tmp.list

sort  ${out_dir}/grch38_tmp.list >  ${out_dir}/grch38.list 
sort  ${out_dir}/chm13_tmp.list > ${out_dir}/chm13.list

bcftools merge -m none --force-samples --file-list ${out_dir}/grch38.list | bgzip > ${out_dir_grch38}/grch38.vcf.gz
tabix -p vcf ${out_dir_grch38}/grch38.vcf.gz

truvari collapse -i ${out_dir_grch38}/grch38.vcf.gz -o ${out_dir_grch38}/grch38.merge.vcf -c ${out_dir_grch38}/grch38.collapsed.vcf -f $grch38 \
--pctseq 0.0 --pctsize 0.5 --pctovl 0.5 --refdist 1000000000 --sizemax 1000000000 -k common

truvari vcf2df -i -f ${out_dir_grch38}/grch38.merge.vcf ${out_dir_grch38}/grch38.merge.jl

done