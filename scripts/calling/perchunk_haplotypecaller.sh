#!/bin/bash
#SBATCH -t 4-00:00
#SBATCH -p thinnodes
#SBATCH -c 24
#SBATCH -n 1

module load cesga/2020
module load gatk/4.2.0.0

###################################
## VARIABLE and PATHS definition ##
###################################

scaffold=($(echo ${1}))

i=($(echo ${2}))

# reference genome
ref=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/Lynx_canadensis_Ref/lc4.fa

# folder containing bam files and bamlist
bamfolder=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/LyCaRef_bams

# name of bamlist
bamlist=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/LyCaRef_bams/lp_ll_introgression.bamlist


###############################
## GATK 4.1.1.0 CombineGVCFs ##
###############################

for k in {1..8}
 do
  echo "HaplotypeCaller of ${scaffold}_w${i}_w${k}.bed"
  
  out=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/LyCaRef_vcfs/lp_ll_introgression_${scaffold}_w${i}_w${k}_LyCa_ref.vcf
  bed=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/Lynx_canadensis_Ref/Scaffold_Partitions/${scaffold}_${i}/${scaffold}_w${i}_w${k}.bed
  
  gatk HaplotypeCaller \
   -R ${ref} \
   $(for bam in $(cat ${bamlist}); do echo "-I ${bamfolder}/${bam}";done) \
   -L ${bed} \
   -O ${out} \
   --native-pair-hmm-threads 3 &
done

wait