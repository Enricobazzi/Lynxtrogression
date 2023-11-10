#!/bin/bash

# With this script I want to apply all of the General filters to my VCF dataset.
# These include:

# (1) Repetitive/Low mappability regions
# (2) Indels + Non-biallelic sites
# (3) Non-variant SNPs
# (4) and (5) Hard quality filters, as GATK standard practices

# It will generate an output VCF file at each step.

# As these filters are independent of population and can be
# applied to the complete VCF file directly.

# Usage:
# ./lp_ll_introgression_vcf_filters_1-5.sh

# This will be run on genomics-a cluster of the EBD

##################################
### Before running verify if : ###
##################################

# GATK
# BEDTOOLS
# BCFTOOLS
# are installed and how they should be used (specific to path and installation)

# and files/directories are correct

###################################
## VARIABLE and PATHS definition ##
###################################

# Reference Genome:
REF=/GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.fa

# Output prefix
OUTpre=lp_ll_introgression_LyCa_ref.sorted

# Input VCF File:
INvcf=lp_ll_introgression_LyCa_ref.sorted.vcf

# BED File of Masked regions:
MASKbed=/GRUPOS/grupolince/reference_genomes/lynx_canadensis/repetitive_regions/lc_rep_ALL_scaffold_coord.bed


# Go to where VCF is
cd /GRUPOS/grupolince/LyCaRef_vcfs/

############################################
## (1) Repetitive/Low mappability regions ##
############################################

echo "starting step 1"

# Apply the filter with BedTools subtract
bedtools subtract -a ${INvcf} -b ${MASKbed} -header | uniq > ${OUTpre}.filter1.vcf

######################################
## (2) Indels + Non-biallelic sites ##
######################################

echo "starting step 2"

# Apply the filter with GATK SelectVariants
/opt/gatk-4.1.0.0/gatk SelectVariants \
  -select-type SNP \
  --restrict-alleles-to BIALLELIC \
  -R ${REF} \
  -V ${OUTpre}.filter1.vcf \
  -O ${OUTpre}.filter2.vcf


#########################
## (3) Invariant sites ##
#########################

echo "starting step 3"

# Apply the filter with BCFtools view
bcftools view -e 'INFO/AF=1.00' ${OUTpre}.filter2.vcf > ${OUTpre}.filter3.vcf

######################################
## (4) and (5) Hard quality filters ##
######################################

echo "starting step 4"

# Filter all except for the RanksSums:
/opt/gatk-4.1.0.0/gatk SelectVariants \
  --selectExpressions "QUAL >= 30 && QD >= 2.0 && FS <= 60.0 && MQ >= 40.0" \
  -R ${REF} \
  -V ${OUTpre}.filter3.vcf \
  -O ${OUTpre}.filter4.vcf

echo "starting step 5"

# Filter RankSums with bcftools view:
bcftools view -e 'INFO/MQRankSum<-12.5 | INFO/ReadPosRankSum<-8.0' ${OUTpre}.filter4.vcf \
 > ${OUTpre}.filter5.vcf
