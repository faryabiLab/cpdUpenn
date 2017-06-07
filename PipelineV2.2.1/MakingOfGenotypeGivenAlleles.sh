##V2.2.1
##Run this script as: ./MakingOfGenotypeGivenAlleles.sh <PathForSample> <hg19.ref."Description or Panel Name"> <RecordsProcessing_script_path>
##"Description or Panel Name" could be FFPE_V1.1 or HEME_V1.1 or HEME_V1.2 respective of the sample.
##Change in V1.3 : The path to vcf-sort command has been removed as now the path is included in the .zshrc file located in /home/bdaber01 on biolinux.
##Changes in V2.0 : 1) $1/$1.Depth.below150X.vcf is added to the command line of $RecordsProcessing.
##Change in V2.1.1_hpc: Change the paths of the scripts that are triggered within this script. This is done to match the path on cluster.
##Change in V2.2.1: Removed the vcf_tools_path argument as now vcf-sort could be directly called due to implementation of its module.

#!/bin/sh
## Extract necessary columns from Depth File using awk command.
PathForSample=$1
Sample_name=`echo $PathForSample|cut -d "/" -f6`
refBase_file=$2
RecordsProcessing=$3

sed 1d $PathForSample/$Sample_name.Depth |awk -F "\t" '{OFS=",";print $1,$2,$5}' |sort -k 1,1 > $PathForSample/$Sample_name.Depth.awked.sorted

## join the awked Depth File with the hg19 refbase database (FFPE_V1 or HEME_V1) 
join -t "," -1 1 -2 1 $refBase_file $PathForSample/$Sample_name.Depth.awked.sorted > $PathForSample/$Sample_name.pileupAndDepth.joined

## Run the $RecordsProcessing Script to obtain GenotypeGivenAlleles, AmpliconsWithDepthBelow250X and DepthForAnnotation vcf Files.
python $RecordsProcessing $PathForSample/$Sample_name.pileupAndDepth.joined $PathForSample/$Sample_name.GenotypeGivenAlleles.vcf $PathForSample/$Sample_name.Depth.below250X.vcf $PathForSample/$Sample_name.Depth.below150X.vcf $PathForSample/$Sample_name.DepthForAnnotation.vcf

## sort the GenotypeGivenAlleles and DepthForAnnotation vcf file using vcf-sort from vcf-tools.
vcf-sort -c $PathForSample/$Sample_name.GenotypeGivenAlleles.vcf > $PathForSample/$Sample_name.GenotypeGivenAlleles.sorted.vcf
vcf-sort -c $PathForSample/$Sample_name.DepthForAnnotation.vcf >  $PathForSample/$Sample_name.DepthForAnnotation.sorted.vcf







