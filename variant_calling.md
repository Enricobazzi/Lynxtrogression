### Variant Calling using GATK

Variant calling was performed on the ft3 CESGA server.

```
ref_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/Lynx_canadensis_Ref
bam_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/LyCaRef_bams
gvcf_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/LyCaRef_gvcfs

scaffolds=$(cat ${ref_dir}/autosomic_scaffolds_list.txt)
ref=${ref_dir}/lc4.fa

# samples=(c_lp_sm_0220 c_lp_sm_0474 c_lp_sm_0614 c_ll_ki_0090 c_ll_ur_0202)
# samples=(c_ll_ya_0141 c_ll_ya_0146 c_ll_ca_0240 c_ll_ca_0243)

for sample in ${samples[@]}
 do
  inbam=$(ls ${bam_dir}/*er.bam | grep ${sample})
  
  for chr in ${scaffolds[@]}
   do
    
    bed=${ref_dir}/CHR_BEDs/${chr}_CHR_coordinates.bed
    outgvcf=${gvcf_dir}/${sample}.${chr}.g.vcf

    echo "sbatch call gvcf of $sample $chr"
    sbatch src/calling/call_gvcf_ref_bam_outgvcf_bed.sh \
    $ref \
    $inbam \
    $outgvcf \
    $bed
    
  done
done
```
Combine individual sample files of each chromosome into a unified file for each chromosome.
```
ref_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/Lynx_canadensis_Ref
gvcf_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/LyCaRef_gvcfs

scaffolds=$(cat ${ref_dir}/autosomic_scaffolds_list.txt)
ref=${ref_dir}/lc4.fa

for chr in ${scaffolds[@]}
 do
  echo "combining gvcfs of ${chr}"
  
  ls ${gvcf_dir}/*.${chr}.g.vcf > tmp_${chr}_gvcf.list
  
  sbatch src/calling/combine_gvcfs_ref_list_outgvcf.sh \
    $ref \
    tmp_${chr}_gvcf.list \
    ${gvcf_dir}/demo_inference.${chr}.g.vcf

done

rm tmp_*
```
Genotype the gvcf of each chromosome
```
ref_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/Lynx_canadensis_Ref
gvcf_dir=/mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/LyCaRef_gvcfs

scaffolds=$(cat ${ref_dir}/autosomic_scaffolds_list.txt)
ref=${ref_dir}/lc4.fa

for chr in ${scaffolds[@]}
 do
  
  ingvcf=${gvcf_dir}/demo_inference.${chr}.g.vcf
  outgvcf=${gvcf_dir}/demo_inference.${chr}.genotyped.vcf
  
  echo "genotype gvcf of ${chr}"
  
  sbatch src/calling/genotype_gvcf_ref_ingvcf_outgvcf.sh \
    $ref \
    $ingvcf \
    $outgvcf

done
```
Concatenate the single chromosome VCFs into one
```
module load picard

gvcf_dir=/GRUPOS/grupolince/LyCaRef_gvcfs
vcf_dir=/GRUPOS/grupolince/LyCaRef_vcfs

# Create a VCF file list
ls ${gvcf_dir}/demo_inference.*.genotyped.vcf > data/calling/demo_inference.genotyped.vcf.list

# Concatenate chromosome VCFs into whole genome VCF
/opt/bcftools-1.6/bcftools concat \
  -f data/calling/demo_inference.genotyped.vcf.list \
  --output-type v \
  --output ${vcf_dir}/demo_inference.LyCa_ref.vcf \
  --threads 22

# Sort with bcftools sort:
/opt/bcftools-1.6/bcftools sort -O v \
  -o ${vcf_dir}/demo_inference.LyCa_ref.sorted.vcf \
  ${vcf_dir}/demo_inference.LyCa_ref.vcf

```

I made a list of the bams to be included in the calling:

```
# On CESGA FT3 server go to Canada lynx reference genome BAM files
cd /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/LyCaRef_bams

# Create a bamlist of samples from populations we want to include in our analysis
ls *er.bam | grep -E "ll|lp" | grep -E "ca|sm|ya|vl|ki|ur" | 
  grep -vE "ca_0249|ca_0253|sm_0138|sm_0140|sm_0185|sm_0186|sm_0221|sm_0298|sm_0359" \
  > lp_ll_introgression.bamlist
```

To generate a VCF file from these BAMs, we performed variant calling on each scaffold in parallel, while additionally dividing all of the bigger scaffolds (18 autosomic + X) into 80 chunks each. This way more jobs have to be sbatched, but we have no memory or time problems.

We called scaffolds that could be called without being divided into chunks using the script [perchr_haplotypecaller.sh](scripts/calling/perchr_haplotypecaller.sh)

```
# short non-autosomic (??) scaffolds 
chr_list=($(cat /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/Lynx_canadensis_Ref/non-autosomic_scaffolds_list.txt | grep -vE "mLynCan4_MT|Super_Scaffold_10"))
for chr in ${chr_list[@]}
 do
  echo "sbatching perchr_haplotypecaller for $chr"
  sbatch --mem=120GB perchr_haplotypecaller.sh $chr
done
```

Chunk divided scaffolds were then called using the script [perchunk_haplotypecaller.sh](./perchunk_haplotypecaller.sh)

```
# big autosomic scaffolds (only 10 at the time for cluster job submission limitations):
chr_list=($(cat /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/reference_genomes/Lynx_canadensis_Ref/autosomic_scaffolds_list.txt))
for chr in ${chr_list[@]:0:10}
 do
  for i in {1..10}
   do
    echo "sbatching perchunk_haplotypecaller for $chr w $i"
    sbatch --mem=120GB perchunk_haplotypecaller.sh $chr $i
  done
done

# X chromosome (Super_Scaffold_10)
for i in {1..10}
 do
  echo "sbatching perchunk_haplotypecaller for Super_Scaffold_10 w $i"
  sbatch --mem=120GB perchunk_haplotypecaller.sh Super_Scaffold_10 $i
done
```

Resulting chunk and scaffold VCFs were concatenated using bcf-tools and then sorted

```
cd /mnt/lustre/hsm/nlsas/notape/home/csic/ebd/jgl/lynx_genome/lynx_data/LyCaRef_vcfs

# Create list of vcf files to concatenate
ls lp_ll_introgression_*caffold*.vcf | tr ' ' '\n' > lp_ll_introgression_vcfs.list

# Concatenate chromosome VCFs into whole genome VCF
bcftools concat -f lp_ll_introgression_vcfs.list \
 --output-type v \
 --output lp_ll_introgression_LyCa_ref.unsorted.vcf

# Sort with bcftools sort:
bcftools sort -O v -o lp_ll_introgression_LyCa_ref.sorted.vcf lp_ll_introgression_LyCa_ref.unsorted.vcf
```
