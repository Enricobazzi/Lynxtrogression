#!/bin/bash
#SBATCH -t 4-00:00
#SBATCH -p thinnodes
#SBATCH -c 4
#SBATCH --mail-type=TIME_LIMIT_80
#SBATCH --mail-user=enricobazzical@gmail.com

# setup environment
module load cesga/2020
module load gatk/4.2.0.0


# reference genome
ref=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/Lynx_canadensis_Ref/lc4.fa

# bed of region to haplotype call
bed=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/Lynx_canadensis_Ref/CHR_BEDs/${1}_CHR_coordinates.bed

# folder containing bam files and bamlist
bamfolder=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/LyCaRef_bams

# name of bamlist
bamlist=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/LyCaRef_bams/lp_ll_introgression.bamlist

# output vcf path and name
out=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/LyCaRef_vcfs/lp_ll_introgression_${1}.vcf

gatk HaplotypeCaller \
   -R ${ref} \
   $(for bam in $(cat ${bamlist}); do echo "-I ${bamfolder}/${bam}";done) \
   -L ${bed} \
   -O ${out} \
   --native-pair-hmm-threads 4
