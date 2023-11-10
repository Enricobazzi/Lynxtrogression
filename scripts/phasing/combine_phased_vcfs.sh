#!/bin/bash

# this script will be used to combine the phased VCFs of each chromosome of each population
# into a single VCF file

#module load gcccore/6.4.0 vcftools/0.1.16
module load gcccore/6.4.0 bcftools/1.10.2

# go to folder
cd /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/LyCaRef_vcfs/lp_ll_introgression

# prefix
lp_ll_introgression_filtered_${chr}_ps_duplicate_phased_masked_merged.vcf

# bgzip and tabix all VCFs for bcftools to work
for vcf in $(ls *_masked.vcf)
 do
  echo "${vcf}"
  bgzip ${vcf}
  tabix -p vcf ${vcf}.gz
done

# get list of chromosomes
chr_list=($(cat /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/Ref_Genome_LyCa/big_scaffolds.bed | cut -f1))

# merge each chromosome's vcf
for chr in ${chr_list[@]}
 do
 echo "merging ${chr} VCFs..."
 # merge chromosome vcf using bcftools
 bcftools merge -O "z" $(ls *_masked.vcf.gz | grep "_${chr}_") \
  > lp_ll_introgression_filtered_${chr}_ps_duplicate_phased_masked_merged.vcf.gz
 tabix -p vcf lp_ll_introgression_filtered_${chr}_ps_duplicate_phased_masked_merged.vcf.gz
done

# concatenate the merged chromosome vcfs
module load gcccore/6.4.0 vcftools/0.1.16

vcf-concat $(ls *_merged.vcf.gz) > ../lp_ll_introgression_LyCa_ref.sorted.filter5.phased.vcf
