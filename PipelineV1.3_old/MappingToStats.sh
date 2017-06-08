##V1.3
##Run this script as : ./MappingToStats.sh
##Changes in V1.3 : 1) The MeanCovByPosition, MeanCovByPositionFinal and NormalizedMeanCovByPosition files are now generated separately for HEME and FFPE panels. 2) Run_MasterVarFinal.txt file now contains the cannon_annotated_flagged content of each sample.

#!/bin/sh

# Select necessary columns from SampleSheet.csv
echo "Processing SampleSheet.csv to extract the SampleName,Barcodes and Description (FFPE_V1 or HEME_V1) ...\n"
sed -n '/Data/,//p' SampleSheet.csv | sed 1,2d | awk -F "," '{print $2,$6,$8,$10}' > Samplesheet_Tobeprocessed.txt


# Create a Sub-Directory for each Sample, copy the fastq files into respective sub-directories.
#echo "Creating Sample Subdirectories ...\n"

for i in `awk '{print $1}' Samplesheet_Tobeprocessed.txt`                                                                        
do
mkdir $i
for j in `ls|grep $i |grep fastq`
do
mv $j $i
done
done




# Create a Pipeline Script for each Sub-Directory
echo "Creating the PipelineScript for each sample ...\n"
python SampleProcessing.py Samplesheet_Tobeprocessed.txt


# Run the Pipe for each sample
for i in `awk '{print $1}' Samplesheet_Tobeprocessed.txt`
do
echo "Running the Pipe for Sample $i ...\n"
nohup $i/$i.PipeScript > $i/$i.out 
cat $i/$i.StatSummary >> RunStats.txt 
cat $i.canon_annotated_flagged >> Run_masterVar.txt

if [ -e $i/$i.HEME_V1.2.MeanCovByPosition ];
then
	cat $i/$i.HEME_V1.2.MeanCovByPosition >> Run_HEME_V1.2_MeanCovByPosition.txt
elif [ -e $i/$i.HEME_V1.1.MeanCovByPosition ];
then
	cat $i/$i.HEME_V1.1.MeanCovByPosition >> Run_HEME_V1.1_MeanCovByPosition.txt
elif [ -e $i/$i.FFPE_V1.1.MeanCovByPosition ];
then
	cat $i/$i.FFPE_V1.1.MeanCovByPosition >> Run_FFPE_V1.1_MeanCovByPosition.txt
fi

cat $i/$i.Depth.below250X.amplicons.summary >> Run_AmpliconSummary.txt
done

## Creating Final Run Stats
echo "Creating Final Run Stats File..\n"
sort -k 1,1 -r RunStats.txt|uniq > RunStatsFinal.txt
sort -k 1,1 -r Run_masterVar.txt|uniq > Run_masterVarFinal.txt
sort -k 1,1 -r Run_AmpliconSummary.txt|uniq > Run_AmpliconSummaryFinal.txt

if [ -e Run_HEME_V1.2_MeanCovByPosition.txt ];
then
	sort -k 1,1 -r Run_HEME_V1.2_MeanCovByPosition.txt|uniq > Run_HEME_V1.2_MeanCovByPositionFinal.txt	
	python CombineStats.py Run_HEME_V1.2_MeanCovByPositionFinal.txt RunStatsFinal.txt Run_HEME_V1.2_NormalizedMeanCovByPosition.txt
fi

if [ -e Run_HEME_V1.1_MeanCovByPosition.txt ];
then
	sort -k 1,1 -r Run_HEME_V1.1_MeanCovByPosition.txt|uniq > Run_HEME_V1.1_MeanCovByPositionFinal.txt
	python CombineStats.py Run_HEME_V1.1_MeanCovByPositionFinal.txt RunStatsFinal.txt Run_HEME_V1.1_NormalizedMeanCovByPosition.txt
fi

if [ -e Run_FFPE_V1.1_MeanCovByPosition.txt ];
then
	sort -k 1,1 -r Run_FFPE_V1.1_MeanCovByPosition.txt|uniq > Run_FFPE_V1.1_MeanCovByPositionFinal.txt
	python CombineStats.py Run_FFPE_V1.1_MeanCovByPositionFinal.txt RunStatsFinal.txt Run_FFPE_V1.1_NormalizedMeanCovByPosition.txt
fi



##rm RunStats.txt
##rm Run_masterVar.txt

##name=`basename $PWD`
#cp *.novo.bam ../../sf_CPP_General_Shrey/$name/.
#cp *.fixed.bam ../../sf_CPP_General_Shrey/$name/.
#cp *.fixed.bam.bai ../../sf_CPP_General_Shrey/$name/.
#cp ReadyToRun ../../sf_CPP_General_Shrey/$name/.
#cp FinalDepthAndStats ../../sf_CPP_General_Shrey/$name/.
#cp nohup.out ../../sf_CPP_General_Shrey/$nam
