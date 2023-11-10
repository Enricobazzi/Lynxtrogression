#!/bin/bash

# generate a CSV table with the number of missing genotypes at each SNP for each population
echo "chromsosome, position, lpa, wel, eel, sel" > lp_ll_introgression_perpop_missing_gts.csv

INdir=/GRUPOS/grupolince/LyCaRef_vcfs/lp_ll_introgression

# chromosome and position from any of the VCFs
grep -v "#" ${INdir}/lp_ll_introgression_filtered_wel.vcf | cut -f1-2 > chr_pos.tmp

# adding population's missing gts to table
for pop in lpa wel eel sel
 do
  echo "calculate missing gts in ${pop}"
  nsamples=($(grep -m1 "#CHR" ${INdir}/lp_ll_introgression_filtered_${pop}.vcf | tr '\t' '\n' | grep "_" | wc -l))
  grep -v "#" ${INdir}/lp_ll_introgression_filtered_${pop}.vcf | cut -f8 | cut -d';' -f3 | 
   cut -d'=' -f2 | awk -v nsam="${nsamples}" '{print ((nsam*2)-$1)/2}' > ${pop}_miss.tmp
  
  echo "adding $pop missing gts to table" 
  paste chr_pos.tmp ${pop}_miss.tmp > tmp && mv tmp chr_pos.tmp
done

# combine table
cat chr_pos.tmp | tr '\t' ',' >> lp_ll_introgression_perpop_missing_gts.csv

# remove tmp files
rm *.tmp
