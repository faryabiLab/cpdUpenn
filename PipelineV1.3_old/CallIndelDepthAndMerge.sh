##V1.3
##Run this script as: ./CallIndelDepthAndMerge.sh <SampleName/SampleName> 

#!/bin/sh
 


awk 'BEGIN {print "#CHROM\tINDEL_POS-2\tINDEL_POS\tREF\tALT\tTYPE\tLENGTH"  ;} {if ($5 == "DEL" || $5 == "INS") {OFS="\t"; print "chr"$1,$2-2,$2,$3,$4,$5,length($4)-1 ;} }' $1.combined.filtered.HRR.canon_snpeff.txt | uniq > $1.indels.BED

samtools mpileup $1.novo.sorted.fixed.ontarget.bam -A -D -B -Q 15 -d 1000000 -f /Drive_D/Databases/genomes/hg19_genome.fa -l $1.indels.BED > $1.novo.sorted.fixed.ontarget.bam.mpileup


samtools mpileup $1.novo.sorted.fixed.ontarget.afterFilter.clipped.bam -A -D -B -Q 15 -d 1000000 -f /Drive_D/Databases/genomes/hg19_genome.fa -l $1.indels.BED > $1.novo.sorted.fixed.ontarget.afterFilter.clipped.bam.mpileup


python CallIndelDepth.py $1.indels.BED $1.novo.sorted.fixed.ontarget.bam.mpileup $1.novo.sorted.fixed.ontarget.afterFilter.clipped.bam.mpileup $1.IndelDepthForAnnotation.vcf

bgzip $1.IndelDepthForAnnotation.vcf
bgzip $1.DepthForAnnotation.sorted.vcf
tabix $1.IndelDepthForAnnotation.vcf.gz
tabix $1.DepthForAnnotation.sorted.vcf.gz

vcf-concat $1.DepthForAnnotation.sorted.vcf.gz $1.IndelDepthForAnnotation.vcf.gz | vcf-sort -c > $1.FinalDepthForAnnotation.vcf







#while read line
#do
#IndelType=`echo $line|cut -f6`
#Pos=`echo $line|cut -f3`
#IndelString=`echo $line|cut -f5`
#IndelLength=`echo $line|cut -f7`
#if [ $IndelType == "DEL" ];
#then
#	IndelPos=$Pos
#elif [ $IndelType == "INS" ];
#	IndelPos=`expr $Pos - 1`
#fi


