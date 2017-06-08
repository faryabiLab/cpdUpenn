##V1.2.3.1
##Run this script as : ./MappingToStats.sh
##Changes in V1.2.3.1 : The script now skips the sample for which the pipescript is not found in its respective sample-folder.

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
if [ -s $i/$i.PipeScript ];
then
	echo "Running the Pipe for Sample $i ...\n"
	nohup $i/$i.PipeScript > $i/$i.out 
	cat $i/$i.StatSummary >> RunStats.txt 
	cat $i.canon_annotated >> Run_masterVar.txt
	cat $i/$i.MeanCovByPosition >> Run_MeanCovByPosition.txt
	cat $i/$i.Depth.below250X.amplicons.summary >> Run_AmpliconSummary.txt
else
	echo "No PipeScript found for Sample $i, hence skipping it...\n"
fi
done

## Creating Final Run Stats
echo "Creating Final Run Stats File..\n"
sort -k 1,1 -r RunStats.txt|uniq > RunStatsFinal.txt
sort -k 1,1 -r Run_masterVar.txt|uniq > Run_masterVarFinal.txt
sort -k 1,1 -r Run_MeanCovByPosition.txt|uniq > Run_MeanCovByPositionFinal.txt
sort -k 1,1 -r Run_AmpliconSummary.txt|uniq > Run_AmpliconSummaryFinal.txt
python CombineStats.py Run_MeanCovByPositionFinal.txt RunStatsFinal.txt Run_NormalizedMeanCovByPosition.txt
##rm RunStats.txt
##rm Run_masterVar.txt

##name=`basename $PWD`
#cp *.novo.bam ../../sf_CPP_General_Shrey/$name/.
#cp *.fixed.bam ../../sf_CPP_General_Shrey/$name/.
#cp *.fixed.bam.bai ../../sf_CPP_General_Shrey/$name/.
#cp ReadyToRun ../../sf_CPP_General_Shrey/$name/.
#cp FinalDepthAndStats ../../sf_CPP_General_Shrey/$name/.
#cp nohup.out ../../sf_CPP_General_Shrey/$nam
