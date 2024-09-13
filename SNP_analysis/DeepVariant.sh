#!/bin/bash

## DeepVariant analysis script
## Author: Gavrielatos Marios
## Mayo Clinic, Fl
## Dr. Owen Ross Lab
## Last revised 09/13/2024

#SBATCH --mail-user=surname.name@mail.edu           # CHANGE THIS
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --job-name=deepvariant                      # Job name
#SBATCH --partition=cpu-short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128                           # Number of CPUs required per task

#SBATCH --time=80:00:00
#SBATCH --mem=128G
# #SBATCH --mem-per-cpu=8G
#SBATCH --chdir /path/to/work/directory
#SBATCH --output logs/%x.%N.%j.stdout
#SBATCH --error logs/%x.%N.%j.stderr

########################## COMMENTS ##########################
# CHANGE MAIL USER & CHDIR                                   #
# ALL YOUR PATHS MUST START FROM /research/                  #
# NO SOFT LINKS ALLOWED                                      #
# The scripts will run the analysis of ONE SAMPLE at a time  #
##############################################################

set -e
set -x

#i=$SLURM_ARRAY_TASK_ID

# Directory with your reference data - CHANGE THIS
REF_PATH_DIR=/path/to/Reference/data/director

# Path to reference genome - CHANGE THIS
REF_PATH_grch38=/path/to/Reference/data/director/Ref_genome_grch38.fasta

# Work directory (where temporary files will be stored) - CHANGE THIS
WD_PATH=/path/to/snp_analysis/directory
# Data directory (where your alignment files are stored) - CHANGE THIS
DATA_PATH=/path/to/alignment/directory/bams

#i=$SLURM_ARRAY_TASK_ID

# If you have multiple samples add the path to your folders and the extension of those files - CHANGE THIS
## Example: 
## My files are in the /my/folder/bams directory &
## I want to use the alignments on the grch38 ref genome with the pbmm2 alignment. These files have the extension ".grch38.pbmm2.bam"
## Then:
## bams="/my/folder/bams/*.grch38.pbmm2.bam"
bams="/my/folder/bams/*extension.bam"

for BAM_FILE in $bams; do
# BAM file we want to call SNPs from
filename_bam=$(basename ${BAM_FILE})
# Extract sample name
sample_name="${filename_bam%%.*}"


# Output GVCF file
GVCF_FILE="${WD_PATH}/${sample_name}/${sample_name}.gvcf"
# Output VCF file
VCF_FILE="${WD_PATH}/${sample_name}/${sample_name}.vcf"
mkdir -p "${WD_PATH}/${sample_name}"


# Log file
LOG_FILE="${WD_PATH}/LOG_OUT"
# Where to store temporary files
TEMP_FILE="${WD_PATH}/TEMP_OUT"



module load apptainer
# Run DeepVariant
apptainer exec --unsquash --bind ${REF_PATH_DIR},${WD_PATH},${DATA_PATH} \
/research/bsi/tools/biotools/deepvariant/1.6.0/deepvariant_1.6.0.sif \
/opt/deepvariant/bin/run_deepvariant \
--intermediate_results_dir ${TEMP_FILE} \
--logging_dir ${LOG_FILE} \
--model_type PACBIO \
--num_shards=32 \
--output_gvcf ${GVCF_FILE} \
--output_vcf ${VCF_FILE} \
--reads ${BAM_FILE} \
--ref ${REF_PATH_grch38} \
--sample_name $sample_name

rm -rf ${TEMP_FILE}
rm -rf ${LOG_FILE}

done


# Merge file - 
# If you want to skip mergin comment out the rest of the script
rm -p glnexus.list
for BAM_FILE in $bams; do
# BAM file we want to call SNPs from
filename_bam=$(basename ${BAM_FILE})
sample_name="${filename_bam%%.*}"

echo "${WD_PATH}/${sample_name}/${sample_name}.gvcf" >> glnexus.list

done

# Merged file name - CHANGE THIS
glnexus_name="sample_name"

glnexus_cli \
--config DeepVariantWGS \
--list glnexus.list \
--threads 10 > ${glnexus_name}.bcf

module load bcftools
bcftools view ${glnexus_name}.bcf | bgzip -@ 4 -c > ${glnexus_name}.vcf.gz
tabix ${glnexus_name}.bcf.vcf.gz
rm -rf GLnexus.DB

# Keep only the main chromosomes
bcftools view ${glnexus_name}.vcf.gz --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY > ${glnexus_name}.reduced.vcf.gz