#!/bin/bash

## Truvari analysis
## Author: Gavrielatos Marios
## Mayo Clinic, Fl
## Dr. Owen Ross Lab
## Last revised 09/13/2024

#SBATCH --mail-user=surname.name@mail.edu           # CHANGE THIS
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50
#SBATCH --job-name=truvari
#SBATCH --partition=cpu-short
#SBATCH --nodes=1
#SBATCH --ntasks=100

#SBATCH --time=80:00:00
#SBATCH --mem=100G
#SBATCH --chdir /path/to/output/directory/annotations/truvari
#SBATCH --output logs/%x.%N.%j.stdout
#SBATCH --error logs/%x.%N.%j.stderr

########################## COMMENTS ##########################
# CHANGE MAIL USER & CHDIR                                   #
# Change directory to output dir: ../truvari                 #
# MOVE ALL VCF FILES TO ONE DIRECTORY:                       #
# .../alignments/annotations/VCFs                            #
##############################################################

# Activate packages 
module load conda
source activate /path/to/mambaforge/envs/truvari

out_dir='/home/mayo/m301955/my_data/my_projects/Genome_Assembly/colombia/alignments/truvari'

out_dir_grch38=${out_dir}/grch38

if [ -d "$out_dir_grch38" ]; then rm -Rf $out_dir_grch38; fi
mkdir -p "${out_dir_grch38}"


callers="sniffles2 pbsv cutesv svdss"

# Run for each caller individually
for caller in $callers; do
echo $caller

# Reference genomes
# Directory with your reference data - CHANGE THIS
grch38='/path/to/Reference/data/director/Ref_genome_grch38.fasta'



vcf_path="/path/to/output/directory/annotations/VCFs/"
vcf_path_grch38="grch38/*${caller}*.vcf.gz"



find $vcf_path -type f -name "*grch38*${caller}*.vcf.gz" > ${out_dir}/grch38_tmp.list

sort grch38_tmp.list > grch38.${caller}.list 

bcftools merge -m none --force-samples --file-list ${out_dir}/grch38.${caller}.list | bgzip > ${out_dir_grch38}/grch38.${caller}.vcf.gz
tabix -p vcf ${out_dir_grch38}/grch38.${caller}.vcf.gz

truvari collapse -i ${out_dir_grch38}/grch38.${caller}.vcf.gz -o ${out_dir_grch38}/grch38.${caller}.merge.vcf -c ${out_dir_grch38}/grch38.${caller}.collapsed.vcf -f $grch38 \
-k common --sizemax 3000000

truvari vcf2df -i -f ${out_dir_grch38}/grch38.${caller}.merge.vcf ${out_dir_grch38}/grch38.${caller}.merge.jl

python truvari_df.py ${out_dir_grch38}/grch38_merge.jl $out_dir_truvari


done

