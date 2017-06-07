##V2.2.1
#!/bin/sh
## Execute this script as : ./CreateFinal.sh <RunFolderPath> <SampleName or JobName> <Panel> <PipeScriptsFolderPath>
##Changes in V2.2: RunStatsFinal now contains the Heme_merge stats in the name of HEMEao_V1.3 sample.
##Changes in V2.2.1: 1) Introduced PipeScriptsFolder variable. Now all data would be synced to the pipeline folder under the runfolder under 'FromHPC' folder on isilon.
##		     2) For HEME samples, the name of HEME_merge sample is replaced with the name of the respective HEME_ao sample in the respective variant, stats and depth files. These files are then concated to the final master files. NameChange_and_concat function handles this scenario. 
##		     3) Removed 'bwlimit' parameter on rsync command line in order to avoid slow file transfer.
RunFolderPath=$1

JobName=`echo $2|cut -d ',' -f1`
##Extracting HEME_V1.2 and HEMEao_V1.3 sample names from argument 4 (This extraction would only work when argument 4 = HEME_merge_sampleName-HEME_V1.2_sampleName-HEMEao_V1.3_sampleName)
HemeJobName=`echo $2|cut -s -d ',' -f2`
Heme_aoJobName=`echo $2|cut -s -d ',' -f3`
PipeScriptsFolderPath=$4
PipeScriptsFolder=`echo $PipeScriptsFolderPath|cut -d "/" -f5`
PanelName=$3

RunFolderName=`echo $RunFolderPath|cut -d "/" -f5`
sync_to_demux ()
{
file_to_sync=$1
rsync -pgt --stats --progress $file_to_sync cpdlab@170.212.141.107:/PathCPD/FromHPC/$RunFolderName/$PipeScriptsFolder/
}

NameChange_and_concat()
{
Name1=$1
Name2=$2
InputFile=$3
OutputFile=$4
echo "$Name1 is replaced with $Name2 in $InputFile to give $InputFile'_NameChange'"
sed "s/$Name1/$Name2/g" $InputFile > $InputFile'_NameChange'
if [ -s $InputFile'_NameChange' ];
then
	cat $InputFile'_NameChange'|sed 1d >> $OutputFile
	sync_to_demux $OutputFile
else
	echo "$InputFile'_NameChange' is either not found or is empty, hence, it did not get concated to $OutputFile"
fi
}

echo "Concating StatSummary file for sample $JobName to the RunStatsFinal File... \n"
	
if [ -s $RunFolderPath/$JobName/$JobName.StatSummary ];
then
	if [ "$PanelName" != "HEME_V1.2_HEMEao_V1.3" ] && [ "$HemeJobName" == "" ] && [ "$Heme_aoJobName" == "" ];
	then
		cat $RunFolderPath/$JobName/$JobName.StatSummary|sed 1d >> $RunFolderPath/RunStatsFinal.txt
		sync_to_demux $RunFolderPath/RunStatsFinal.txt
	elif [ "$PanelName" == "HEME_V1.2_HEMEao_V1.3" ] && [ -n "$HemeJobName" ] && [ -n "$Heme_aoJobName" ];
        then
		NameChange_and_concat $JobName $Heme_aoJobName $RunFolderPath/$JobName/$JobName.StatSummary $RunFolderPath/RunStatsFinal.txt
	fi
else
	echo "$JobName/$JobName.StatSummary is either not found or is empty. It did not get concated to the RunStatsFinal.txt file"
fi

echo "Concating canon_annotated_flagged file for sample $JobName to the Run_masterVarFinal File... \n"
	
if [ -s $RunFolderPath/$JobName.canon_annotated_flagged ];
then
	if [ "$PanelName" != "HEME_V1.2_HEMEao_V1.3" ] && [ "$HemeJobName" == "" ] && [ "$Heme_aoJobName" == "" ];
	then
		cat $RunFolderPath/$JobName.canon_annotated_flagged|sed 1d >> $RunFolderPath/Run_masterVarFinal.txt
		sync_to_demux $RunFolderPath/Run_masterVarFinal.txt
	elif [ "$PanelName" == "HEME_V1.2_HEMEao_V1.3" ] && [ -n "$HemeJobName" ] && [ -n "$Heme_aoJobName" ];
	then
		NameChange_and_concat $JobName $Heme_aoJobName $RunFolderPath/$JobName.canon_annotated_flagged $RunFolderPath/Run_masterVarFinal.txt
	fi
else
	echo "$JobName.canon_annotated_flagged is either not found or is empty. It did not get concated to the Run_masterVarFinal.txt file"
fi
	
echo "Concating MeanCovByPosition file for sample $JobName to the MeanCovByPositionFinal file, generating its Normalized Mean_Cov file and then concating that to the NormalizedMeanCovByPositionFinal file..."

if [ -s $RunFolderPath/$JobName/$JobName.$PanelName.MeanCovByPosition ];
then
	sed -i 's/_mean_cvg//g' $RunFolderPath/$JobName/$JobName.$PanelName.MeanCovByPosition
	if [ "$PanelName" != "HEME_V1.2_HEMEao_V1.3" ] && [ "$HemeJobName" == "" ] && [ "$Heme_aoJobName" == "" ];
	then
		cat $RunFolderPath/$JobName/$JobName.$PanelName.MeanCovByPosition|sed 1d >> $RunFolderPath/Run_$PanelName"_MeanCovByPositionFinal".txt
		sync_to_demux $RunFolderPath/Run_$PanelName"_MeanCovByPositionFinal".txt
	elif [ "$PanelName" == "HEME_V1.2_HEMEao_V1.3" ] && [ -n "$HemeJobName" ] && [ -n "$Heme_aoJobName" ];
	then
		NameChange_and_concat $JobName $Heme_aoJobName $RunFolderPath/$JobName/$JobName.$PanelName.MeanCovByPosition $RunFolderPath/Run_$PanelName"_MeanCovByPositionFinal".txt
	fi
	
	if [ -s $RunFolderPath/$JobName/$JobName.StatSummary ];
	then
        	python $PipeScriptsFolderPath/CombineStats.py $RunFolderPath/$JobName/$JobName.$PanelName.MeanCovByPosition $RunFolderPath/$JobName/$JobName.StatSummary $RunFolderPath/$JobName/$JobName.$PanelName.NormalizedMeanCovByPosition
		if [ -s $RunFolderPath/$JobName/$JobName.$PanelName.NormalizedMeanCovByPosition ];
		then
			if [ "$PanelName" != "HEME_V1.2_HEMEao_V1.3" ] && [ "$HemeJobName" == "" ] && [ "$Heme_aoJobName" == "" ];
			then
				cat $RunFolderPath/$JobName/$JobName.$PanelName.NormalizedMeanCovByPosition|sed 1d >> $RunFolderPath/Run_$PanelName"_NormalizedMeanCovByPositionFinal".txt
				sync_to_demux $RunFolderPath/Run_$PanelName"_NormalizedMeanCovByPositionFinal".txt
			elif [ "$PanelName" == "HEME_V1.2_HEMEao_V1.3" ] && [ -n "$HemeJobName" ] && [ -n "$Heme_aoJobName" ];
			then
				NameChange_and_concat $JobName $Heme_aoJobName $RunFolderPath/$JobName/$JobName.$PanelName.NormalizedMeanCovByPosition $RunFolderPath/Run_$PanelName"_NormalizedMeanCovByPositionFinal".txt
			fi
		else
			echo "$JobName/$JobName.$PanelName.NormalizedMeanCovByPosition is either not found or is empty. It did not get concated to the Run_$PanelName'_NormalizedMeanCovByPositionFinal'.txt file"
		fi
	else
        	echo "Cannot create $JobName/$JobName.$PanelName.NormalizedMeanCovByPosition.This is because $JobName/$JobName.StatSummary is either not found or is empty"
	fi
else
	echo "$JobName/$JobName.$PanelName.MeanCovByPosition is either not found or is empty. It did not get concated to the Run_$PanelName'_MeanCovByPositionFinal'.txt file. Consequently $JobName/$JobName.$PanelName.NormalizedMeanCovByPosition did not get created."
fi

echo "Concating amplicon summary files for sample $JobName to the final AmpliconSummary file"
if [ -s $RunFolderPath/$JobName/$JobName.Depth.below250X.amplicons*.summary ];
then
	cat $RunFolderPath/$JobName/$JobName.Depth.below250X.amplicons*.summary >> $RunFolderPath/Run_$PanelName"_AmpliconSummary_250X_final".txt
	sync_to_demux $RunFolderPath/Run_$PanelName"_AmpliconSummary_250X_final".txt
else
	echo "$JobName/$JobName.Depth.below250X.amplicons.summary (or amplicons_sorted.summary) is either not found or is empty. It did not get concated to the Run_$PanelName'_AmpliconSummary_250X_final'.txt file"
fi
if [ -s $RunFolderPath/$JobName/$JobName.Depth.below150X.amplicons*.summary ];
then
	cat $RunFolderPath/$JobName/$JobName.Depth.below150X.amplicons*.summary >> $RunFolderPath/Run_$PanelName"_AmpliconSummary_150X_final".txt
	sync_to_demux $RunFolderPath/Run_$PanelName"_AmpliconSummary_150X_final".txt
else
	echo "$JobName/$JobName.Depth.below150X.amplicons.summary (or amplicons_sorted.summary) is either not found or is empty. It did not get concated to the Run_$PanelName'_AmpliconSummary_150X_final'.txt file"
fi

echo "Copying files important for variant review process..."
for i in `ls $RunFolderPath/$JobName*.ba*|grep -v 'SAM_BAM_check'`
do
	sync_to_demux $i
done
sync_to_demux $RunFolderPath/$JobName.canon_annotated
sync_to_demux $RunFolderPath/$JobName.canon_annotated_flagged
sync_to_demux $RunFolderPath/$JobName/$JobName.combined.vcf
sync_to_demux $RunFolderPath/SampleSheet_ToBeProcessed.txt
sync_to_demux $RunFolderPath/SampleSheet_forMergeAnalysis.txt
