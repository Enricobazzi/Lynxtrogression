#!/bin/bash

# script to take output of shapeit phasing and remove imputation by reverting missing data

# Usage:
# ./gt_masker_pop_chr_vcf.sh <population> <chromosome>

# input directory
INdir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/LyCaRef_vcfs/lp_ll_introgression

# population
pop=($(echo "${1}"))

# chromosome
chr=($(echo "${2}"))

# prefix
PREfix=${INdir}/lp_ll_introgression_filtered_${pop}_${chr}

# original VCF with missing data
orig_VCF=${PREfix}_ps.vcf

# imputed VCF without missing data
imput_VCF=${PREfix}_ps_duplicate_phased.vcf

# number of samples
nsamples=($(grep -m1 "#CHR" ${orig_VCF} | tr '\t' '\n' | grep "_" | wc -l))

echo "masking missing genotypes from ${pop}'s phased VCF of chromosome ${chr}"

# extract each sample's GT in orig and imput
START=1
END=${nsamples}

for (( c=$START; c<=$END; c++ ))
 do
  echo "extracting sample_${pop}_${chr} $c orig GTs"
  grep -v "#" ${orig_VCF} | rev | cut -f1-${nsamples} | rev | cut -f$c | cut -d':' -f1 \
    > ${INdir}/sample_${pop}_${chr}_${c}.orig.tmp.gts

  echo "extracting sample_${pop}_${chr} $c imput GTs"
  grep -v "#" ${imput_VCF} | rev | cut -f1-${nsamples} | rev | cut -f$c \
    > ${INdir}/sample_${pop}_${chr}_${c}.imput.tmp.gts

  echo "masking sample_${pop}_${chr} $c missing GTs"
  paste ${INdir}/sample_${pop}_${chr}_${c}.orig.tmp.gts ${INdir}/sample_${pop}_${chr}_${c}.imput.tmp.gts |
    awk '{FS="\t"; OFS="\t"; if ($1 != "./.") print $2; else print $1; }' \
    > ${INdir}/sample_${pop}_${chr}_${c}.miss.tmp.gts
  
  if [ $c == 1 ]
   then
    echo "creating file for all samples starting with sample_${pop}_${chr} $c"
    cp ${INdir}/sample_${pop}_${chr}_${c}.miss.tmp.gts ${INdir}/${pop}_${chr}.miss.tmp.gts
   else
    echo "joining sample_${pop}_${chr} $c with rest"
    paste ${INdir}/${pop}_${chr}.miss.tmp.gts ${INdir}/sample_${pop}_${chr}_${c}.miss.tmp.gts \
     > ${INdir}/${pop}_${chr}.tmp && mv ${INdir}/${pop}_${chr}.tmp ${INdir}/${pop}_${chr}.miss.tmp.gts
  fi
done

# build header 
echo "building header for ${pop}'s phased and masked VCF of chromosome ${chr}"
grep "##" ${imput_VCF} > ${PREfix}_ps_duplicate_phased_masked.vcf

grep -m1 "#CHR" ${imput_VCF} | tr '\t' '\n' | grep -v "_2$" | paste -s -d '\t' \
  >> ${PREfix}_ps_duplicate_phased_masked.vcf
  
# build rest of vcf
echo "building rest of ${pop}'s phased and masked VCF of chromosome ${chr}"

ncols=($(grep -m1 "#CHR" ${PREfix}_ps_duplicate_phased_masked.vcf | tr '\t' '\n' | grep -v "_" | wc -l))

paste <(grep -v "#" ${imput_VCF} | cut -f1-${ncols}) <(cat ${INdir}/${pop}_${chr}.miss.tmp.gts) \
  >> ${PREfix}_ps_duplicate_phased_masked.vcf

# remove tmp files
echo "removing temporary files"
rm ${INdir}/*${pop}_${chr}*.tmp.gts
