#!/bin/bash

### Before running verify if:
### GATK 
### BEDTOOLS
### BCFTOOLS
### are installed and how they should be called

# Purpose: Filter out repetitive regions, indels, non-biallelic sites, invariant sites, and sites with low quality.
# Input:
#   - $1: reference genome
#   - $2: input vcf
#   - $3: mask bed file (repetitive and low mappability regions)
# Output:
#   - $2.filter5.vcf: filtered vcf
# Usage:
#   - bash filter_1to5_ref_invcf_mask.sh ref.fa in.vcf mask.bed

# fasta reference
ref=$1
# input vcf
invcf=$2
# mask bed file (repetitive and low mappability regions)
mask=$3
# basename of input vcf
vcf_basename=$(basename "${invcf}" .vcf)
# path of input vcf
vcf_dir=$(dirname "${invcf}")
# prefix of vcfs
vcf_pre="${vcf_path}/${vcf_basename}"

echo "starting step 1: filtering out repetitive regions"

# Apply the filter with BedTools subtract
bedtools subtract -a ${invcf} -b ${mask} -header | uniq > ${vcf_pre}.filter1.vcf

echo "starting step 2: filtering out indels and non-biallelic sites"

# Apply the filter with GATK SelectVariants
/opt/gatk-4.1.0.0/gatk SelectVariants \
  -select-type SNP \
  --restrict-alleles-to BIALLELIC \
  -R ${ref} \
  -V ${vcf_pre}.filter1.vcf \
  -O ${vcf_pre}.filter2.vcf

echo "starting step 3: filtering out invariant sites"

# Apply the filter with BCFtools view
/opt/bcftools-1.6/bcftools view \
  -e 'INFO/AF=1.00'  \
  ${vcf_pre}.filter2.vcf \
  > ${vcf_pre}.filter3.vcf

echo "starting step 4 and 5: filtering out sites with low quality"

# Filter all except for the RanksSums:
/opt/gatk-4.1.0.0/gatk SelectVariants \
  --selectExpressions "QUAL >= 30 && QD >= 2.0 && FS <= 60.0 && MQ >= 40.0" \
  -R ${ref} \
  -V ${vcf_pre}.filter3.vcf \
  -O ${vcf_pre}.filter4.vcf

# Filter RankSums with bcftools view:
/opt/bcftools-1.6/bcftools view -e 'INFO/MQRankSum<-12.5 | INFO/ReadPosRankSum<-8.0' \
  ${vcf_pre}.filter4.vcf \
  > ${vcf_pre}.filter5.vcf
