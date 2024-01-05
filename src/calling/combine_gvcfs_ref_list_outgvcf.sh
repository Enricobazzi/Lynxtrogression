#!/bin/bash
#SBATCH --time=00-4:00:00
#SBATCH --mem=10G
#SBATCH --output=logs/calling/combine_gvcf-%j.out
#SBATCH --error=logs/calling/combine_gvcf-%j.err
#SBATCH --cpus-per-task=1

# Purpose: Combine gvcfs using a list of gvcfs as input and a single gvcf as output
# Input:
#   - $1: reference genome
#   - $2: file with list of gvcfs (one per line)
#   - $3: output gvcf file
# Output:
#   - $3: output gvcf file
# Usage:
#   - bash combine_gvcfs_ref_list_outgvcf.sh ref.fa gvcf.list out.g.vcf

module load gatk

ref=$1
gvcf_list=$2
gvcf_array=($(cat $gvcf_list))
outgvcf=$3

gatk CombineGVCFs \
   -R $ref \
   $(for i in ${gvcf_array[@]}; do echo "--variant ${i}";done) \
   -O $outgvcf
