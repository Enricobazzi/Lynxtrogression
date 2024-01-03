#!/bin/bash
#SBATCH --time=00-6:00:00
#SBATCH --mem=20G
#SBATCH --output=logs/calling/gvcf-%j.out
#SBATCH --error=logs/calling/gvcf-%j.err
#SBATCH --cpus-per-task=4

# Purpose: Call gvcf from bam file using a bed file
# Input:
#   - $1: reference genome
#   - $2: input bam file
#   - $3: output gvcf file
#   - $4: bed file
# Output:
#   - $3: output gvcf file
# Usage:
#   - bash call_gvcf_ref_bam_outgvcf_bed.sh ref.fa in.bam out.g.vcf bed.bed

module load gatk

ref=$1
inbam=$2
outgvcf=$3
bed=$4

gatk HaplotypeCaller  \
   -R $ref \
   -I $inbam \
   -O $outgvcf \
   -L $bed \
   --native-pair-hmm-threads 4 \
   -ERC GVCF
