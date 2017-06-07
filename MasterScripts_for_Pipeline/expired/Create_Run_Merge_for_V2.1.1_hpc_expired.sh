##V2.1.1_hpc
#!/bin/sh
## Execute this script as : ./Create_Run_Merge.sh <RunFolderPath> <PipeScriptsFolderPath> <Min_Memory> <Max_Memory> <Threads>
## This wrapper script executes a pipeline for a RunFolder.
RunFolderPath=$1
PipeScriptsFolderPath=$2
Min_Mem=$3
Max_Mem_init=$4
Max_Mem_init_1=`echo $Max_Mem_init|cut -d 'G' -f 1`
Max_Mem=$((Max_Mem_init_1 * 1024))
Num_Threads=$5
RunFolderName=`echo $RunFolderPath|cut -d "/" -f 5`
## copying fastqs and SampleSheet from the RunFolder on SeqOut
echo "Copying fastqs and SampleSheet from RunFolder on SeqOut.."
rsync -pgt --stats --progress --bwlimit=10000 cpdlab@170.212.141.107:/PathCPD/illumina/SeqOut/$RunFolderName/Data/Intensities/BaseCalls/*fastq.gz /project/cpdlab/RunFolders/$RunFolderName/
rsync -pgt --stats --progress --bwlimit=10000 cpdlab@170.212.141.107:/PathCPD/illumina/SeqOut/$RunFolderName/SampleSheet.csv /project/cpdlab/RunFolders/$RunFolderName/
##Checking SampleSheet
if [ -s $RunFolderPath/SampleSheet.csv ];
then
	for i in `sed -n '/Data/,//p' $RunFolderPath/SampleSheet.csv|sed 1d| awk -F "," '{print NF}'`
	do
	if [ $i != 12 ];
	then
		echo "SampleSheet error!..Check the values in each column."
		exit 71
	fi
	done
else
	echo "SampleSheet.csv not found..Hence, exiting.."
	exit 78
fi
##Processing SampleSheet.
echo "Processing SampleSheet.."
sed 's/_/-/g' $RunFolderPath/SampleSheet.csv|sed 's/-V1/_V1/g'|sed 's/-V2/_V2/g'|sed -n '/Data/,//p'| sed 1,2d | awk -F "," '{print $2,$6,$8,$10}' > $RunFolderPath/SampleSheet_ToBeProcessed.txt
sed 's/_/-/g' $RunFolderPath/SampleSheet.csv|sed 's/-V1/_V1/g'|sed 's/-V2/_V2/g'|sed -n '/Data/,//p'| sed 1,2d | awk -F "," '{print $1,$2,$6,$8,$10}' > $RunFolderPath/SampleSheet_forMergeAnalysis.txt
chmod 775 $RunFolderPath/SampleSheet_ToBeProcessed.txt
chmod 775 $RunFolderPath/SampleSheet_forMergeAnalysis.txt

##Create a Folder for each Sample, copy the fastq files into respective folders.
if [ -s $RunFolderPath/SampleSheet_ToBeProcessed.txt ];
then
	echo "Creating Sample Subdirectories ...\n"
	for i in `egrep -v '_skip|-skip' $RunFolderPath/SampleSheet_ToBeProcessed.txt| awk '{print $1}'`
	do
		mkdir $RunFolderPath/$i
		mv $RunFolderPath/$i*fastq.gz $RunFolderPath/$i
	done
	# Create a Pipeline Script for each Sub-Directory
	echo "Creating the PipelineScript for each sample ...\n"
	python $PipeScriptsFolderPath/SampleProcessing.py $RunFolderPath/SampleSheet_ToBeProcessed.txt $RunFolderPath $Min_Mem $Max_Mem_init $Num_Threads
else
	echo "$RunFolderPath/SampleSheet_ToBeProcessed.txt NOT found, therefore, no samples to process..."
	exit 88
fi

# Create a temp directory
mkdir $RunFolderPath/tmp

## Defining a CountCat function
CountCat ()
{
Panel=$1
cat $PipeScriptsFolderPath/File_Headers/$Panel"_MeanCov_Header" >> $RunFolderPath/Run_$Panel"_MeanCovByPositionFinal".txt
cat $PipeScriptsFolderPath/File_Headers/$Panel"_MeanCov_Header" >> $RunFolderPath/Run_$Panel"_NormalizedMeanCovByPositionFinal".txt
}

## Before Starting the Pipe, writing headers for all the final master files
cat $PipeScriptsFolderPath/File_Headers/RunStatsFinal_Header >> $RunFolderPath/RunStatsFinal.txt
cat $PipeScriptsFolderPath/File_Headers/Run_masterVarFinal_Header >> $RunFolderPath/Run_masterVarFinal.txt

for i in `egrep -v '_skip|-skip|EGFRvIII' $RunFolderPath/SampleSheet_ToBeProcessed.txt|awk '{print $4}'|sort|uniq`
do
CountCat $i
done

# Run the Pipe for each sample
for i in `egrep -v '_skip|-skip' $RunFolderPath/SampleSheet_ToBeProcessed.txt| awk '{print $1","$4}'`
do
SampleName=`echo $i|cut -d "," -f1`
PanelName=`echo $i|cut -d "," -f2`
#PName=`echo $PanelName|sed 's/\./_/g'`
#myString=$PName"_Count"
if [ -s $RunFolderPath/$SampleName/$SampleName.PipeScript ];
then
	chmod 770 $RunFolderPath/$SampleName/$SampleName.PipeScript
	echo "Running the Pipe for Sample $SampleName ...\n"
	bsub -J $SampleName -n $Num_Threads -R "span[hosts=1]" -M $(($Max_Mem + 12288)) -o $RunFolderPath/$SampleName/$SampleName.PipeScript.out -e $RunFolderPath/$SampleName/$SampleName.PipeScript.err $RunFolderPath/$SampleName/$SampleName.PipeScript
	bsub -J CreateFinal_$SampleName -w "done($SampleName)" -o $RunFolderPath/$SampleName/$SampleName.CreateFinal.out -e $RunFolderPath/$SampleName/$SampleName.CreateFinal.err sh CreateFinal.sh $RunFolderPath $SampleName $PanelName $PipeScriptsFolderPath
fi
done

##Create a Pipeline Script for each Merged Sub-Directory
if [ -s $RunFolderPath/SampleSheet_forMergeAnalysis.txt ];
then
	echo "Creating the PipelineScript for each HEME_merged sample ...\n"
	python $PipeScriptsFolderPath/Merge.py $RunFolderPath/SampleSheet_forMergeAnalysis.txt $RunFolderPath $Min_Mem $Max_Mem_init $Num_Threads
else
	echo "$RunFolderPath/SampleSheet_forMergeAnalysis.txt NOT found, therefore no HEME_merge_sample to process..."
fi

if [ -s $RunFolderPath/Run_HEME_SamplePairsInfo.txt ];
then
	cat $PipeScriptsFolderPath/File_Headers/HEME_V1.2_HEMEao_V1.3_MeanCov_Header >> $RunFolderPath/Run_HEME_V1.2_HEMEao_V1.3_MeanCovByPositionFinal.txt
	cat $PipeScriptsFolderPath/File_Headers/HEME_V1.2_HEMEao_V1.3_MeanCov_Header >> $RunFolderPath/Run_HEME_V1.2_HEMEao_V1.3_NormalizedMeanCovByPositionFinal.txt
	#HEME_V2_3_Count=0
        for i in `awk '{print $1}' $RunFolderPath/Run_HEME_SamplePairsInfo.txt`
        do
        HEME_sample=`echo $i|cut -d ',' -f1`
        HEMEao_sample=`echo $i|cut -d ',' -f2`
        HEME_merge_sample=`echo $i|cut -d ',' -f3`
        echo $HEME_sample
        echo $HEMEao_sample
        echo $HEME_merge_sample
	if [ -s $RunFolderPath/$HEME_merge_sample/$HEME_merge_sample.PipeScript ];
	then
		chmod 770 $RunFolderPath/$HEME_merge_sample/$HEME_merge_sample.PipeScript
		bsub -J $HEME_merge_sample -w "done($HEME_sample) && done($HEMEao_sample)" -n $Num_Threads -R "span[hosts=1]" -M $(($Max_Mem + 12288)) -o $RunFolderPath/$HEME_merge_sample/$HEME_merge_sample.PipeScript.out -e $RunFolderPath/$HEME_merge_sample/$HEME_merge_sample.PipeScript.err $RunFolderPath/$HEME_merge_sample/$HEME_merge_sample.PipeScript
		bsub -J CreateFinal_$HEME_merge_sample -w "done($HEME_merge_sample)" -e $RunFolderPath/$HEME_merge_sample/$HEME_merge_sample.CreateFinal.err -o $RunFolderPath/$HEME_merge_sample/$HEME_merge_sample.CreateFinal.out sh CreateFinal.sh $RunFolderPath $HEME_merge_sample,$HEME_sample,$HEMEao_sample HEME_V1.2_HEMEao_V1.3 $PipeScriptsFolderPath
	fi
        done
else
	echo "$RunFolderPath/Run_HEME_SamplePairsInfo.txt NOT found, therefore no HEME_merge_sample to process..."
fi
