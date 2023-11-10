#!/bin/bash

# output dir
OUTdir=/GRUPOS/grupolince/LyCaRef_vcfs/lp_ll_introgression/gmaps

# input vcf
INvcf=/GRUPOS/grupolince/LyCaRef_vcfs/lp_ll_introgression_LyCa_ref.sorted.filter5.vcf

# list of chromosomes
chr_list=($(cat /GRUPOS/grupolince/reference_genomes/lynx_canadensis/big_scaffolds.bed | cut -f1))

# calculating genetic map of each chromosome from VCF 
for chr in ${chr_list[@]}
 do
  echo "calculating genetic map of ${chr} from ${INvcf}"
  echo "${OUTdir}/${chr}.gmap"
  grep -v "#" ${INvcf} | grep -w ${chr} | cut -f1,2 |
  awk '{ print $2, $1 }' |
  awk {'if ( NR==1 ) print $1, $2, 0; else print $1, $2, $1-p, ($1-p)*0.0000019; p=$1'} |
  awk 'BEGIN{print "pos", "chr", "cM"} {sum+=$4} {print $1, $2, sum}' |
  tr ' ' '\t' > ${OUTdir}/${chr}.gmap
done
