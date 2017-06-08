##V1.4.2
##Run this script as : ./MappingToStats.sh
## Changes in V1.3 : 1) The MeanCovByPosition, MeanCovByPositionFinal and NormalizedMeanCovByPosition files are now generated separately for HEME,FFPE,mPPP,bPPP and CEBPa panels.
##                     2) If...else conditions are added to make pipe run more smoothly. The sample folders not containing the pipescript are skipped during the pipe-run.
##                     3) The "echo" statement indicating the creation of Sample Subdirectories is now made active in this version.
##                     4) Run_masterVar and hence Run_MasterVarFinal now contains the contents of cannon_annotated_flagged file instead of cannon_annotated file for each sample.
## Changes in V1.3.1 : 1) Removed the step that converted SampleSheet.csv to SampleSheet_Tobeprocessed.txt as that is now handled upfornt by 'Pipe.sh' script.
##                     2) Added an extra step to remove the fastqs of samples whose description contains "_skip" word in them.
## Changes in V1.3.2 : 1) Added an extra step to create "tmp" folder within the main RunFolder as the new version of Novosort requires it to be created prior to running any sample through its pipescript.
## Changes in V1.4 :   1) Added steps to process Merge Pipe for Heme samples. Merge pipe merges the data for HEME_V1.2 and HEMEao_V1.3 so as to give Merge Results (i.e. HEME_V2.3 results).
## Changes in V1.4.1 : 1) The raw run files and the final run files now get created separately for HEME_V2.3 (i.e. merge) data.
## Changes in V1.4.2 : 1) Included 'if' conditions at "Merge.py execution step" as well as at "RunStatsFinal_HEME_V2.3.txt, Run_masterVarFinal_HEME_V2.3.txt and Run_AmpliconSummaryFinal.txt generation steps" to eliminate the errors due to absence of 'SampleSheet_forMergeAnalysis.txt' and 'RunStats_HEME_V2.3.txt, Run_masterVarFinal_HEME_V2.3 and Run_AmpliconSummaryFinal.txt' files respectively.
##		       2) Established more control at the "Run pipe for each sample" stage by including the 'egrep -v' command on the 'for loop initialization line' to exclude running the pipe for skipped samples.
##                     3) Included an important step to generate final variant upload files, i.e. Run_masterVarFinal_UP.txt and Run_masterVarFinal_HEME_V2.3_UP.txt files.
#!/bin/sh

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

# Create a temp directory
mkdir tmp

# Run the Pipe for each sample
for i in `egrep -v '_skip|-skip' Samplesheet_Tobeprocessed.txt|awk '{print $1}'`
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
	
	if [ -e $i/$i.HEME_V1.2.MeanCovByPosition ];
	then
		cat $i/$i.HEME_V1.2.MeanCovByPosition >> Run_HEME_V1.2_MeanCovByPosition.txt
	elif [ -e $i/$i.HEME_V1.1.MeanCovByPosition ];
	then
		cat $i/$i.HEME_V1.1.MeanCovByPosition >> Run_HEME_V1.1_MeanCovByPosition.txt
	elif [ -e $i/$i.FFPE_V1.1.MeanCovByPosition ];
	then
		cat $i/$i.FFPE_V1.1.MeanCovByPosition >> Run_FFPE_V1.1_MeanCovByPosition.txt
	elif [ -e $i/$i.mPPP.MeanCovByPosition ];
	then
		cat $i/$i.mPPP.MeanCovByPosition >> Run_mPPP_MeanCovByPosition.txt
	elif [ -e $i/$i.bPPP.MeanCovByPosition ];
	then
		cat $i/$i.bPPP.MeanCovByPosition >> Run_bPPP_MeanCovByPosition.txt
	elif [ -e $i/$i.CEBPa.MeanCovByPosition ];
	then
			cat $i/$i.CEBPa.MeanCovByPosition >> Run_CEBPa_MeanCovByPosition.txt
	fi
	
	
	
	if [ -s $i/$i.Depth.below250X.amplicons.summary ];
	then
		cat $i/$i.Depth.below250X.amplicons.summary >> Run_AmpliconSummary.txt
	fi


else
	echo "No PipeScript found for Sample $i, hence skipping it...\n"


fi


done




# Create a Pipeline Script for each Merged Sub-Directory
if [ -s SampleSheet_forMergeAnalysis.txt ];
then
	echo "Creating the PipelineScript for each sample ...\n"
	python Merge.py SampleSheet_forMergeAnalysis.txt
fi

# Run the Pipe for each sample
for i in `ls|grep merged`
do
if [ -d $i ];
then
	if [ -s $i/$i.PipeScript ];
	then
		echo "Running the Pipe for Sample $i ...\n"
	
		nohup $i/$i.PipeScript > $i/$i.out
	
	
		echo "Preparing raw Run Files ... \n"
	
		if [ -e $i/$i.StatSummary ];
		then
			cat $i/$i.StatSummary >> RunStats_HEME_V2.3.txt
		fi
	
		if [ -s $i.canon_annotated_flagged ];
		then
			cat $i.canon_annotated_flagged >> Run_masterVar_HEME_V2.3.txt
		fi
		
		
		
		## Concating the Panel-specific MeanCovByPosition file of each Sample to the final Run_(Panel)_MeanCovByPosition.txt.
		if [ -e $i/$i.HEME_V1.2_HEMEao_V1.3.MeanCovByPosition ];
		then
			cat $i/$i.HEME_V1.2_HEMEao_V1.3.MeanCovByPosition >> Run_HEME_V2.3_MeanCovByPosition.txt
		fi
	
	
	
		if [ -s $i/$i.Depth.below250X.amplicons_sorted.summary ];
		then
			cat $i/$i.Depth.below250X.amplicons_sorted.summary >> Run_AmpliconSummary_HEME_V2.3.txt
		fi


	else
		echo "No PipeScript found for Sample $i, hence skipping it...\n"


	fi
fi
done






## Creating Final Run Files

echo "Creating Final Run Stats File..\n"
sort -k 1,1 -r RunStats.txt|uniq > RunStatsFinal.txt

echo "Creating Final MasterVar file..\n"
sort -k 1,1 -r Run_masterVar.txt|uniq > Run_masterVarFinal.txt

echo "Creating Final Amplicon Summary File..\n"
sort -k 1,1 -r Run_AmpliconSummary.txt|uniq > Run_AmpliconSummaryFinal.txt

if [ -s RunStats_HEME_V2.3.txt ];
then
	echo "Creating Final Run Stats File for HEME_V2.3..\n"
	sort -k 1,1 -r RunStats_HEME_V2.3.txt|uniq > RunStatsFinal_HEME_V2.3.txt
fi

if [ -s Run_masterVar_HEME_V2.3.txt ];
then
	echo "Creating Final MasterVar file..\n"
	sort -k 1,1 -r Run_masterVar_HEME_V2.3.txt|uniq > Run_masterVarFinal_HEME_V2.3.txt
fi

if [ -s Run_AmpliconSummary_HEME_V2.3.txt ];
then
	echo "Creating Final Amplicon Summary File for HEME_V2.3..\n"
	sort -k 1,1 -r Run_AmpliconSummary_HEME_V2.3.txt|uniq > Run_AmpliconSummaryFinal_HEME_V2.3.txt
fi

if [ -s Run_HEME_SamplePairsInfo.txt ];
then
	echo "Generating the final variant upload files..\n" 
	python PrepareFinalUpload.py Run_masterVarFinal.txt Run_masterVarFinal_HEME_V2.3.txt Run_masterVarFinal_UP.txt Run_masterVarFinal_HEME_V2.3_UP.txt Run_HEME_SamplePairsInfo.txt
else
	echo "Run_HEME_SamplePairsInfo not found...Hence Run_masterVarFinal_UP.txt and Run_masterVarFinal_HEME_V2.3_UP.txt files were not generated"
fi






echo "Creating Final MeanCovByPosition and NormalizedMeanCovByPosition Files for each Panel..\n"

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
if [ -e Run_mPPP_MeanCovByPosition.txt ];
then
	sort -k 1,1 -r Run_mPPP_MeanCovByPosition.txt|uniq > Run_mPPP_MeanCovByPositionFinal.txt
	python CombineStats.py Run_mPPP_MeanCovByPositionFinal.txt RunStatsFinal.txt Run_mPPP_NormalizedMeanCovByPosition.txt
fi
if [ -e Run_bPPP_MeanCovByPosition.txt ];
then
	sort -k 1,1 -r Run_bPPP_MeanCovByPosition.txt|uniq > Run_bPPP_MeanCovByPositionFinal.txt
	python CombineStats.py Run_bPPP_MeanCovByPositionFinal.txt RunStatsFinal.txt Run_bPPP_NormalizedMeanCovByPosition.txt
fi
if [ -e Run_CEBPa_MeanCovByPosition.txt ];
then
	sort -k 1,1 -r Run_CEBPa_MeanCovByPosition.txt|uniq > Run_CEBPa_MeanCovByPositionFinal.txt
	python CombineStats.py Run_CEBPa_MeanCovByPositionFinal.txt RunStatsFinal.txt Run_CEBPa_NormalizedMeanCovByPosition.txt
fi
if [ -e Run_HEME_V2.3_MeanCovByPosition.txt ];
then
	sort -k 1,1 -r Run_HEME_V2.3_MeanCovByPosition.txt|uniq > Run_HEME_V2.3_MeanCovByPositionFinal.txt	
	python CombineStats.py Run_HEME_V2.3_MeanCovByPositionFinal.txt RunStatsFinal_HEME_V2.3.txt Run_HEME_V2.3_NormalizedMeanCovByPosition.txt
fi


echo "Deleting the fastqs from the folders of the skipped Samples.."

for i in `grep '_skip' SampleSheet_Tobeprocessed.txt|awk '{print $1}'`
do
rm $i/*fastq.gz
done
