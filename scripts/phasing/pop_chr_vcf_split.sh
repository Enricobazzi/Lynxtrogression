#!/bin/bash

# output dir
OUTdir=/GRUPOS/grupolince/LyCaRef_vcfs/lp_ll_introgression

# input vcf
INvcf=/GRUPOS/grupolince/LyCaRef_vcfs/lp_ll_introgression_LyCa_ref.sorted.filter5.vcf

# list of chromosomes
chr_list=($(cat /GRUPOS/grupolince/reference_genomes/lynx_canadensis/big_scaffolds.bed | cut -f1))

# table of pop information
pop_table=/GRUPOS/grupolince/LyCaRef_vcfs/lp_ll_introgression/lp_ll_introgression_populations.txt

# list of sample populations
pops=($(cat ${pop_table} | cut -f2 | sort -u))

for pop in ${pops[@]}
 do
  # get list of sample population names for the study population
  pop_list=($(grep ${pop} ${pop_table} | cut -f1))
  # get list of samples from the population
  sample_list=($(grep -f <(echo ${pop_list[@]} | tr ' ' '\n') <(grep -m1 "#CHR" ${INvcf} | tr '\t' '\n')))
  
  echo "extracting ${pop_list[@]} samples of pop ${pop}:"
  echo "${sample_list[@]}"
  
  # get VCF including only selected samples
  /opt/gatk-4.1.0.0/gatk SelectVariants \
   -R /GRUPOS/grupolince/reference_genomes/lynx_canadensis/lc4.fa \
   -V ${INvcf} \
   $(for j in ${sample_list[@]}; do echo "-sn ${j}";done) \
   -O ${OUTdir}/lp_ll_introgression_filtered_${pop}.vcf
   
  # split VCF into single chrs
  for chr in ${chr_list[@]}
   do
    echo "extracting ${chr} from ${pop} VCF"
    grep -w ${chr} /GRUPOS/grupolince/reference_genomes/lynx_canadensis/big_scaffolds.bed
    bedtools intersect -header -a ${OUTdir}/lp_ll_introgression_filtered_${pop}.vcf \
     -b <(grep -w ${chr} /GRUPOS/grupolince/reference_genomes/lynx_canadensis/big_scaffolds.bed) \
     > ${OUTdir}/lp_ll_introgression_filtered_${pop}_${chr}.vcf
  done
done
