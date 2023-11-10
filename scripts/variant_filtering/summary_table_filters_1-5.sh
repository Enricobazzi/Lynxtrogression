#!/bin/bash

echo "Variant Filtering 1 - 5 Summary!"

# Output prefix
OUTpre=/GRUPOS/grupolince/LyCaRef_vcfs/lp_ll_introgression_LyCa_ref.sorted
 
# Output file
OUTlog=lp_ll_introgression_LyCa_ref.filters_1-5.log.csv

echo "step,name,e_vars,f_vars" > ${OUTlog}

# Print number of variants before filtering to log:
ST1start=($(grep -v "#" ${OUTpre}.vcf | wc -l))
echo "initial number of variants : ${ST1start}"
echo "0,start,${ST1start},0" >> ${OUTlog}

############################################
## (1) Repetitive/Low mappability regions ##
############################################

# Print number of variants after filtering to log:
ST1end=$(grep -v "#" ${OUTpre}.filter1.vcf | wc -l)
ST1filtered="$(echo "$ST1start - $ST1end" | bc)"

echo "1,low_map,${ST1end},${ST1filtered}" >> ${OUTlog}
echo "Repetitive/Low mappability regions - number of variants filtered : $ST1filtered"
echo "Repetitive/Low mappability regions - final number of variants : $ST1end"

######################################
## (2) Indels + Non-biallelic sites ##
######################################

# Print number of variants after filtering to log:
ST2end=$(grep -v "#" ${OUTpre}.filter2.vcf | wc -l)
ST2filtered="$(echo "$ST1end - $ST2end" | bc)"

echo "2,indel_bial,${ST2end},${ST2filtered}" >> ${OUTlog}
echo "Indels + Non-biallelic sites - number of variants filtered : ${ST2filtered}"
echo "Indels + Non-biallelic sites - final number of variants : ${ST2end}"

#########################
## (3) invariant sites ##
#########################

# Print number of variants after filtering to log:
ST3end=$(grep -v "#" ${OUTpre}.filter3.vcf | wc -l)
ST3filtered="$(echo "$ST2end - $ST3end" | bc)"

echo "3,invar,${ST3end},${ST3filtered}" >> ${OUTlog}
echo "Lynx genus wide exclusive substitutions - final number of variants : ${ST3end}"
echo "Lynx genus wide exclusive substitutions - number of variants filtered : ${ST3filtered}"

######################################
## (4) and (5) Hard quality filters ##
######################################

# Print number of variants after filtering to log:
ST4end=$(grep -v "#" ${OUTpre}.filter4.vcf | wc -l)
ST4filtered="$(echo "$ST3end - $ST4end" | bc)"

echo "4,gatk_qual1,${ST4end},${ST4filtered}" >> ${OUTlog}
echo "Hard quality except for the RanksSums - final number of variants : ${ST4end}"
echo "Hard quality except for the RanksSums - number of variants filtered : ${ST4filtered}"

# Print number of variants after filtering to log:
ST5end=$(grep -v "#" ${OUTpre}.filter5.vcf | wc -l)
ST5filtered="$(echo "$ST4end - $ST5end" | bc)"

echo "5,gatk_qual2,${ST5end},${ST5filtered}" >> ${OUTlog}
echo "RanksSums excluded by Hard quality - final number of variants : ${ST5end}"
echo "RanksSums excluded by Hard quality - number of variants filtered : ${ST5filtered}"
