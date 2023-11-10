#!/bin/bash

# script to duplicate the samples in the vcf obtained by whatshap phase set generation 

# Usage:
# ./duplicate_pop_chr_vcf.sh <population> <chromosome>

# input directory
INdir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/LyCaRef_vcfs/lp_ll_introgression

# population
pop=($(echo "${1}"))

# chromosome
chr=($(echo "${2}"))

# prefix
PREfix=${INdir}/lp_ll_introgression_filtered_${pop}_${chr}_ps

# number of samples
nsamples=($(grep -m1 "#CHR" ${PREfix}.vcf | tr '\t' '\n' | grep "_" | wc -l))


# copy the header except for the table line
echo "creating duplicated VCF header"
grep "##" ${PREfix}.vcf > ${PREfix}_duplicate.vcf

# create a new VCF table header (#CHR line)
# with double the samples (new ones are called the same + "_2")
# sed command removes trailing \t 
echo "adding samples names"
paste <(grep "#CHR" ${PREfix}.vcf) \
 <(paste -d '_' <(grep "#CHR" ${PREfix}.vcf | tr '\t' '\n' | grep "_") <(yes "2" | head -n ${nsamples}) | 
   tr '\n' '\t') | sed 's/^[ \t]*//;s/[ \t]*$//' \
 >> ${PREfix}_duplicate.vcf

# extract extract sample's genotypes
echo "extracting genotypes of ${nsamples} samples"
grep -v "#" ${PREfix}.vcf | rev | cut -f1-${nsamples} | rev > tmp.gts
  
# paste the vcf without the header with the genotype and add it to the duplicated vcf
echo "pasting the GTs and completing duplicated VCF"
paste <(grep -v "#" ${PREfix}.vcf) <(cat tmp.gts) >> ${PREfix}_duplicate.vcf

# remove temporary file
echo "removing tmp file"
rm tmp.gts
