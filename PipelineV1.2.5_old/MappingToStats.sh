##V1.2.4
##Run this script as : ./MappingToStats.sh
## Changes in V1.2.4 : 1) The MeanCovByPosition, MeanCovByPositionFinal and NormalizedMeanCovByPosition files are now generated separately for HEME and FFPE panels.
##                     2) If...else conditions are added to make pipe run more smoothly. The sample folders not containing the pipescript are skipped during the pipe-run.
##                     3) The "echo" statement indicating the creation of Sample Subdirectories is now made active in this version.
##                     4) Run_masterVar and hence Run_MasterVarFinal now contains the contents of cannon_annotated_flagged file instead of cannon_annotated file for each sample.

#!/bin/sh

# Select necessary columns from SampleSheet.csv
echo "Processing SampleSheet.csv to extract the SampleName,Barcodes and Description (FFPE_V1 or HEME_V1) ...\n"
sed -n '/Data/,//p' SampleSheet.csv | sed 1,2d | awk -F "," '{print $2,$6,$8,$10}' > Samplesheet_Tobeprocessed.txt


# Create a Sub-Directory for each Sample, copy the fastq files into respective sub-directories.
echo "Creating Sample Subdirectories ...\n"

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
if [ -s $i/$i.PipeScript ];
then
	echo "Running the Pipe for Sample $i ...\n"
	
	nohup $i/$i.PipeScript > $i/$i.out
	
	
	echo "Preparing raw Run Files ... \n"
	
	if [ -e $i/$i.StatSummary ];
	then
		cat $i/$i.StatSummary >> RunStats.txt
	fi
	
	if [ -s $i.canon_annotated_flagged ];
	then
		cat $i.canon_annotated_flagged >> Run_masterVar.txt
	fi
	
	
	
	## Concating the Panel-specific MeanCovByPosition file of each Sample to the final Run_(Panel)_MeanCovByPosition.txt.
	
	if [ -e $i/$i.HEME_V2.0.MeanCovByPosition ];
	then
		cat $i/$i.HEME_V2.0.MeanCovByPosition >> Run_HEME_V2.0_MeanCovByPosition.txt
	elif [ -e $i/$i.HEME_V1.2.MeanCovByPosition ];
	then
		cat $i/$i.HEME_V1.2.MeanCovByPosition >> Run_HEME_V1.2_MeanCovByPosition.txt
	elif [ -e $i/$i.HEME_V1.1.MeanCovByPosition ];
	then
		cat $i/$i.HEME_V1.1.MeanCovByPosition >> Run_HEME_V1.1_MeanCovByPosition.txt
	elif [ -e $i/$i.FFPE_V1.1.MeanCovByPosition ];
	then
		cat $i/$i.FFPE_V1.1.MeanCovByPosition >> Run_FFPE_V1.1_MeanCovByPosition.txt
	fi
	
	
	
	if [ -s $i/$i.Depth.below250X.amplicons.summary ];
	then
		cat $i/$i.Depth.below250X.amplicons.summary >> Run_AmpliconSummary.txt
	fi


else
	echo "No PipeScript found for Sample $i, hence skipping it...\n"


fi


done




## Creating Final Run Files

echo "Creating Final Run Stats File..\n"
sort -k 1,1 -r RunStats.txt|uniq > RunStatsFinal.txt

echo "Creating Final MasterVar file..\n"
sort -k 1,1 -r Run_masterVar.txt|uniq > Run_masterVarFinal.txt

echo "Creating Final Amplicon Summary File..\n"
sort -k 1,1 -r Run_AmpliconSummary.txt|uniq > Run_AmpliconSummaryFinal.txt





echo "Creating Final MeanCovByPosition and NormalizedMeanCovByPosition Files for each Panel..\n"

if [ -e Run_HEME_V2.0_MeanCovByPosition.txt ];
then
	sort -k 1,1 -r Run_HEME_V2.0_MeanCovByPosition.txt|uniq > Run_HEME_V2.0_MeanCovByPositionFinal.txt	
	python CombineStats.py Run_HEME_V2.0_MeanCovByPositionFinal.txt RunStatsFinal.txt Run_HEME_V2.0_NormalizedMeanCovByPosition.txt
fi
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
