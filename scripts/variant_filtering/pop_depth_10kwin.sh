#!/bin/bash

module load samtools

cd /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/LyCaRef_bams

pop=($(echo "${1}"))

bamlist=($(echo ${pop}.bamlist))

chr_list=($(cat /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/Lynx_canadensis_Ref/big_scaffolds.bed | cut -f1))

for chr in ${chr_list[@]}
 do
  length=($(grep "${chr}" /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/Lynx_canadensis_Ref/big_scaffolds.bed | cut -f3))
  n_win=($(echo "($length/10000) +1" | bc))
  echo "chromosome ${chr} has a length of ${length} and a total of ${n_win} 10Kbp windows"
 
  # bed name
  bed=lp_ll_introgression_depth/${pop}_${chr}_10kb_depth.bed
  # rm bed if it already exists
  if [ -f "$bed" ]
   then
    echo "$bed already exists! removing it..."
    rm ${bed}
  fi
  
  START=1
  END=$n_win

  for (( c=$START; c<=$END; c++ ))
   do
    echo "window $c :"
    if [ $c = 1 ]
     then
      ws=1
      we=10000
     else
      ws=($(echo "$we + 1" | bc))
      we=($(echo "$ws + 10000 - 1" | bc))
    fi
    # build bed
    echo "$chr:$ws-$we"
    samtools depth -a -r ${chr}:${ws}-${we} -q 0 -Q 0 $(cat ${bamlist}) | cut -f3- | tr '\t' '\n' | 
     awk 'BEGIN{SUM=0}{ SUM+=$0 }END{print SUM/10000}' | 
     awk -v chr="${chr}" -v ws="${ws}" -v we="${we}" '{print chr, ws-1, we, $1}' |
     tr ' ' '\t' >> ${bed}
    
  done
done
