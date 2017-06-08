##V1.2.4
##Run this script as: ./MakingOfGenotypeGivenAlleles.sh <SampleName> <hg19.ref."Description or Panel Name">
##"Description or Panel Name" could be FFPE_V1.1 or HEME_V1.1 or HEME_V1.2 respective of the sample.
##Change in V1.2.4 : The path to vcf-sort command has been removed as now the path is included in the .zshrc file located in /home/bdaber01 on biolinux. 

#!/bin/sh
## Extract necessary columns from Depth File using awk command.
sed 1d $1/$1.Depth |awk -F "\t" '{OFS=",";print $1,$2,$5}' |sort -k 1,1 > $1/$1.Depth.awked.sorted

## join the awked Depth File with the hg19 refbase database (FFPE_V1 or HEME_V1) 
join -t "," -1 1 -2 1 $2 $1/$1.Depth.awked.sorted > $1/$1.pileupAndDepth.joined

## Run the RecordsProcessing.py Script to obtain GenotypeGivenAlleles, AmpliconsWithDepthBelow250X and DepthForAnnotation vcf Files.
python RecordsProcessing.py $1/$1.pileupAndDepth.joined $1/$1.GenotypeGivenAlleles.vcf $1/$1.Depth.below250X.vcf $1/$1.DepthForAnnotation.vcf

## sort the GenotypeGivenAlleles and DepthForAnnotation vcf file using vcf-sort from vcf-tools.
vcf-sort -c $1/$1.GenotypeGivenAlleles.vcf > $1/$1.GenotypeGivenAlleles.sorted.vcf
vcf-sort -c $1/$1.DepthForAnnotation.vcf >  $1/$1.DepthForAnnotation.sorted.vcf







