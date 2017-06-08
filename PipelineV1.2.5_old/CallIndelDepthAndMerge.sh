##V1.2.4
##Run this script as: ./CallIndelDepthAndMerge.sh <SampleName/SampleName> 

#!/bin/sh
 


awk 'BEGIN {print "#CHROM\tINDEL_POS-2\tINDEL_POS\tREF\tALT\tTYPE\tLENGTH"  ;} {if ($5 == "DEL" || $5 == "INS") {OFS="\t"; print "chr"$1,$2-2,$2,$3,$4,$5,length($4)-1 ;} }' $1.combined.filtered.HRR.canon_snpeff.txt | uniq > $1.indels.BED

samtools mpileup $1.novo.sorted.fixed.ontarget.bam -A -D -B -Q 15 -d 1000000 -f /Drive_D/Databases/genomes/hg19_genome.fa -l $1.indels.BED > $1.novo.sorted.fixed.ontarget.bam.mpileup


samtools mpileup $1.novo.sorted.fixed.ontarget.afterFilter.clipped.bam -A -D -B -Q 15 -d 1000000 -f /Drive_D/Databases/genomes/hg19_genome.fa -l $1.indels.BED > $1.novo.sorted.fixed.ontarget.afterFilter.clipped.bam.mpileup


python CallIndelDepth.py $1.indels.BED $1.novo.sorted.fixed.ontarget.bam.mpileup $1.novo.sorted.fixed.ontarget.afterFilter.clipped.bam.mpileup $1.IndelDepthForAnnotation.vcf

vcf-sort -c $1.IndelDepthForAnnotation.vcf > $1.IndelDepthForAnnotation.sorted.vcf

bgzip -f $1.IndelDepthForAnnotation.sorted.vcf
bgzip -f $1.DepthForAnnotation.sorted.vcf
tabix -p vcf -f $1.IndelDepthForAnnotation.sorted.vcf.gz
tabix -p vcf -f $1.DepthForAnnotation.sorted.vcf.gz

vcf-concat $1.DepthForAnnotation.sorted.vcf.gz $1.IndelDepthForAnnotation.sorted.vcf.gz | vcf-sort -c > $1.FinalDepthForAnnotation.vcf


