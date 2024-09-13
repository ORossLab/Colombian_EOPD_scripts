## Set up Oakvar for SNP analysis
## Author: Gavrielatos Marios
## Mayo Clinic, Fl
## Dr. Owen Ross Lab
## Last revised 09/13/2024

########################## COMMENTS #############
# Run this in the terminal                      #
#################################################

export OV_ROOT_DIR=/path/to/your/dir/oakvar
export OV_LOGS_DIR=/path/to/your/dir/oakvar/oakvar_logs
export OV_LIFTOVER_DIR=/path/to/your/dir/oakvar/conf/liftover
export OV_MODULES_DIR=/path/to/your/dir/oakvar/modules
export OV_JOBS_DIR=/path/to/your/dir/oakvar/jobs
export OV_LOG_DIR=/path/to/your/dir/oakvar/logs
ov system setup

ov module install cadd variantreport hpo clinvar go dbsnp gnomad3 gencode ncbigene varity_r spliceai sift revel regulomedb polyphen2 phi grasp dann_coding fathmm_xf_coding mutation_assessor phdsnpg regulomedb dann fathmm_xf

# If you want to run ONE vcf file:
# Run an annotation job
# -a option controls annotation sources and -t option report formats.
out_dir=/path/to/output/directory/
vcf="/path/to/target/vcf/file/target.vcf"
ov run $vcf -a cadd variantreport hpo clinvar go dbsnp gnomad3 gencode ncbigene varity_r spliceai sift revel regulomedb polyphen2 phi grasp dann_coding fathmm_xf_coding mutation_assessor phdsnpg regulomedb dann fathmm_xf -t vcf -d ${out_dir} --mp 5
sqlite=/path/to/output/directory/target.vcf.sqlite
output=/path/to/output/directory/oakvar_target.reduced
ov report $sqlite -t tsv -s $out_dir



# If you want to run MULTIPLE vcf files:

bams="/path/to/bam/files/*bam"
for bam in $bams; do
filename_bam=$(basename ${bam})
sample_id="${filename_bam%%.*}"
{
# Run an annotation job
# -a option controls annotation sources and -t option report formats.
out_dir=/path/to/output/directory/
vcf="/path/to/target/vcf/file/target.vcf"
ov run $vcf -a cadd variantreport hpo clinvar go dbsnp gnomad3 gencode ncbigene varity_r spliceai sift revel regulomedb polyphen2 phi grasp dann_coding fathmm_xf_coding mutation_assessor phdsnpg regulomedb dann fathmm_xf -t vcf -d ${out_dir} --mp 5
sqlite=/path/to/output/directory/target.vcf.sqlite
output=/path/to/output/directory/oakvar_target.reduced
ov report $sqlite -t tsv -s $out_dir
} &
done
sleep 0.1 # For sequential output
echo "Waiting for processes to finish"
echo "All processes finished"