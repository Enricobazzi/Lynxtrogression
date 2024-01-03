#!/bin/bash
#SBATCH --output=logs/alignment/qualimap-%j.out
#SBATCH --error=logs/alignment/qualimap-%j.err
#SBATCH --time=6:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=10

module load qualimap

bam=$1
out_dir=$2
sample=$(echo $bam | rev | cut -d'/' -f1 | rev | cut -d'_' -f1,2,3,4)

qualimap bamqc \
  -bam ${bam} \
  --java-mem-size=19G \
  -outfile ${sample}_qualimap.html \
  -outformat html \
  -outdir ${out_dir}/${sample}
