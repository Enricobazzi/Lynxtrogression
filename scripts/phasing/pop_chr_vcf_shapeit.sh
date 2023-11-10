#!/bin/bash
#SBATCH -t 0-02:00
#SBATCH -p thinnodes
#SBATCH -c 1
#SBATCH -n 1
#SBATCH --mem=80GB

# This script will run SHAPEIT4 for phasing the zipped duplicated phase set VCFs of 
# a specific chromosome and population pair

# Usage:
# pop_chr_vcf_shapeit.sh <population> <chromosome>

module load cesga/2020
module load gcccore/system shapeit4/4.2.1

# input directory
INdir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/LyCaRef_vcfs/lp_ll_introgression

# population
pop=($(echo "${1}"))

# chromosome
chr=($(echo "${2}"))

# prefix
PREfix=lp_ll_introgression_filtered_${pop}_${chr}_ps_duplicate

# gmap
gmap=${INdir}/gmaps/${chr}.gmap

# SHAPEIT4 command
shapeit4.2 \
 --input ${INdir}/${PREfix}.vcf.gz \
 --map ${gmap} \
 --region ${chr} \
 --use-PS 0.0001 \
 --output ${INdir}/${PREfix}_phased.vcf \
 --mcmc-iterations 10b,1p,1b,1p,1b,1p,1b,1p,10m
