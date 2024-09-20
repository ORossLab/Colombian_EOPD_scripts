#!/bin/bash

####################################################
# Run with:
# chmod +x homozygosity_test.sh
# ./homozygosity_test.sh file1.txt
####################################################


# Check if one VCF file is provided as input
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 input_vcf"
    exit 1
fi

# Input VCF file
input_vcf="$1"

# Get the directory of the input VCF file
input_dir=$(dirname "$input_vcf")

# Compress and index the VCF file
bgzip -@15 -f "$input_vcf"
tabix "${input_vcf}.gz"

# Set output file names in the same directory as the input VCF file
check_left_vcf="${input_dir}/check.prkn.left.vcf"
check_right_vcf="${input_dir}/check.prkn.right.vcf"

# Extract specific regions from the VCF file
tabix "${input_vcf}.gz" chr6:161785908-161920242 > "$check_left_vcf"
tabix "${input_vcf}.gz" chr6:162102085-162201252 > "$check_right_vcf"

# Set file names for processing
file1=$(realpath "$check_left_vcf")
file2=$(realpath "$check_right_vcf")

# Combine the files into one temporary file
combined_file=$(mktemp)
cat "$file1" "$file2" > "$combined_file"

# Total number of lines in the combined file
total_lines=$(wc -l < "$combined_file")

# Count occurrences of each genotype in the combined file
count_1_1=$(grep -o "1/1" "$combined_file" | wc -l)
count_0_0=$(grep -o "0/0" "$combined_file" | wc -l)
count_0_1=$(grep -o "0/1" "$combined_file" | wc -l)
count_1_2=$(grep -o "1/2" "$combined_file" | wc -l)
count_2_2=$(grep -o "2/2" "$combined_file" | wc -l)
count_dot_dot=$(grep -o "\./\." "$combined_file" | wc -l)

# Calculate percentages with more decimal places
percentage_1_1=$(echo "scale=4; $count_1_1 / $total_lines * 100" | bc)
percentage_0_0=$(echo "scale=4; $count_0_0 / $total_lines * 100" | bc)
percentage_0_1=$(echo "scale=4; $count_0_1 / $total_lines * 100" | bc)
percentage_1_2=$(echo "scale=4; $count_1_2 / $total_lines * 100" | bc)
percentage_2_2=$(echo "scale=4; $count_2_2 / $total_lines * 100" | bc)
percentage_dot_dot=$(echo "scale=4; $count_dot_dot / $total_lines * 100" | bc)

# Calculate combined counts and percentages for homozygosity
combined_count_homozygosity=$((count_1_1 + count_0_0 + count_2_2))
combined_percentage_homozygosity=$(echo "scale=4; $combined_count_homozygosity / $total_lines * 100" | bc)

combined_count_heterozygosity=$((count_0_1 + count_1_2))
combined_percentage_heterozygosity=$(echo "scale=4; $combined_count_heterozygosity / $total_lines * 100" | bc)



# Output the results
echo "Total lines: $total_lines"
echo "Occurrences of 1/1: $count_1_1 ($percentage_1_1%)"
echo "Occurrences of 0/0: $count_0_0 ($percentage_0_0%)"
echo "Occurrences of 0/1: $count_0_1 ($percentage_0_1%)"
echo "Occurrences of 1/2: $count_1_2 ($percentage_1_2%)"
echo "Occurrences of 2/2: $count_2_2 ($percentage_2_2%)"
echo "Occurrences of ./.: $count_dot_dot ($percentage_dot_dot%)"

echo "Occurrences of 1/1 or 0/0 or 2/2 (homozygosity): $combined_count_homozygosity ($combined_percentage_homozygosity%)"
echo "Occurrences of 0/1 or 1/2 (heterozygosity): $combined_count_heterozygosity ($combined_percentage_heterozygosity%)"

# Clean up temporary files
rm "$combined_file"

