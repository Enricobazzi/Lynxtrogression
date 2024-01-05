#!/bin/bash
#SBATCH --time=00-6:00:00
#SBATCH --mem=10G
#SBATCH --output=logs/calling/genotype_gvcf-%j.out
#SBATCH --error=logs/calling/genotype_gvcf-%j.err
#SBATCH --cpus-per-task=1

# Purpose: GenotypeGVCFs on a single gvcf file
# Input:
#   - $1: reference genome
#   - $2: input gvcf file
#   - $3: output gvcf file
# Output:
#   - $3: output gvcf file
# Usage:
#   - bash genotype_gvcf_ref_ingvcf_outgvcf.sh ref.fa in.g.vcf out.g.vcf

module load gatk

ref=$1
ingvcf=$2
outgvcf=$3

gatk  GenotypeGVCFs \
  -R $ref \
  -V $ingvcf \
  -O $outgvcf
