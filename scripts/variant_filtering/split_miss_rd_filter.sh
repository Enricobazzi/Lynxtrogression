#!/bin/bash

# this script will be used to split the phased VCF into the VCFs of the different population
# pairs (lpa-wel, lpa-eel, lpa-sel) and filter missingness and read depth accordingly

cd /GRUPOS/grupolince/LyCaRef_vcfs

# prefix
prefix=lp_ll_introgression_LyCa_ref.sorted.filter5.phased.fixed

# for each population of eurasian lynx
for pop in wel eel sel
 do
  
  # define population pair
  pop_pair=($(echo "lpa-${pop}"))
  # define population names in population pair
  pop_names=($(grep "${pop}" lp_ll_introgression/lp_ll_introgression_populations.txt | cut -f1 | paste -s -d '|'))
  # define samples in population pair
  pop_pair_samples=($(grep -m1 "#CHR" ${prefix}.vcf | tr '\t' '\n' | grep -E "sm|${pop_names}"))
  
  # extract VCF of population pair
  echo "extracting VCF of ${pop_pair} population pair..."
  /opt/gatk-4.1.0.0/gatk SelectVariants \
   -R /GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.fa \
   -V ${prefix}.vcf \
   $(for j in ${pop_pair_samples[@]}; do echo "-sn ${j}";done) \
   -O ${prefix}.${pop_pair}.vcf

  # apply the filter based on missingness
  echo "filtering ${pop_pair} VCF removing SNPs with >=15% missing data"
  bedtools subtract -header \
   -a ${prefix}.${pop_pair}.vcf \
   -b lp_ll_introgression/filter_beds/lpa_miss_filter.bed |
   bedtools subtract -header \
   -a stdin \
   -b lp_ll_introgression/filter_beds/${pop}_miss_filter.bed \
   > ${prefix}.${pop_pair}.miss.vcf
 
  # apply the filter based on read depth
  echo "filtering ${pop_pair} VCF removing windows with an excess of read depth"
  bedtools subtract -header \
   -a ${prefix}.${pop_pair}.miss.vcf \
   -b lp_ll_introgression/filter_beds/lpa_rd_filter.bed |
   bedtools subtract -header \
   -a stdin \
   -b lp_ll_introgression/filter_beds/${pop}_rd_filter.bed \
   > ${prefix}.${pop_pair}.miss.rd_fil.vcf

done
