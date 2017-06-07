##V2.2.1
#!/bin/sh
## Execute this script as : ./CreateFinal.sh <RunFolderPath> <SampleName or JobName> <Panel> <PipeScriptsFolderPath>
##Changes in V2.2: RunStatsFinal now contains the Heme_merge stats in the name of HEMEao_V1.3 sample.
##Changes in V2.2.1: Introduced PipeScriptsFolder variable. Now all data would be synced to the pipeline folder under the runfolder under 'FromHPC' folder on isilon.
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
rsync -pgt --stats --progress --bwlimit=10000 $file_to_sync cpdlab@170.212.141.107:/PathCPD/FromHPC/$RunFolderName/$PipeScriptsFolder/
}

echo "Concating StatSummary file for sample $JobName to the RunStatsFinal File... \n"
	
if [ -s $RunFolderPath/$JobName/$JobName.StatSummary ];
then
	#if [ "$PanelName" != "HEME_V1.2_HEMEao_V1.3" ] && [ "$HemeJobName" == "" ] && [ "$Heme_aoJobName" == "" ];
	#then
		cat $RunFolderPath/$JobName/$JobName.StatSummary|sed 1d >> $RunFolderPath/RunStatsFinal.txt
		sync_to_demux $RunFolderPath/RunStatsFinal.txt
	#else
		#echo "Since $JobName is an Heme_merge sample, its stats would not get concated to the RunStatsFinal file..."
	#fi
	if [ "$PanelName" == "HEME_V1.2_HEMEao_V1.3" ] && [ -n "$HemeJobName" ] && [ -n "$Heme_aoJobName" ] && [ -s $RunFolderPath/RunStatsFinal.txt ];
        then
		echo "$JobName is replaced with $Heme_aoJobName in the stats entry that contains $JobName..."
                sed -i "s/$JobName/$Heme_aoJobName/g" $RunFolderPath/RunStatsFinal.txt
                sync_to_demux $RunFolderPath/RunStatsFinal.txt
	fi
else
	echo "$JobName/$JobName.StatSummary is either not found or is empty. It did not get concated to the RunStatsFinal.txt file"
fi

echo "Concating canon_annotated_flagged file for sample $JobName to the Run_masterVarFinal File... \n"
	
if [ -s $RunFolderPath/$JobName.canon_annotated_flagged ];
then
	cat $RunFolderPath/$JobName.canon_annotated_flagged|sed 1d >> $RunFolderPath/Run_masterVarFinal.txt
	sync_to_demux $RunFolderPath/Run_masterVarFinal.txt

	if [ "$PanelName" == "HEME_V1.2_HEMEao_V1.3" ] && [ -n "$HemeJobName" ] && [ -n "$Heme_aoJobName" ] && [ -s $RunFolderPath/Run_masterVarFinal.txt ];
	then
		#echo "Since $JobName is an Heme_merge sample, all variant entries corresponding to $Heme_aoJobName are removed from Run_masterVarFinal.txt.."
		#sed -i "/$Heme_aoJobName/d" $RunFolderPath/Run_masterVarFinal.txt
		echo "$JobName is replaced with $Heme_aoJobName in all variant lines that contain $JobName..."
		sed -i "s/$JobName/$Heme_aoJobName/g" $RunFolderPath/Run_masterVarFinal.txt
		sync_to_demux $RunFolderPath/Run_masterVarFinal.txt
	fi
else
	echo "$JobName/$JobName.canon_annotated_flagged is either not found or is empty. It did not get concated to the Run_masterVarFinal.txt file"
fi
	
echo "Concating MeanCovByPosition file for sample $JobName to the MeanCovByPositionFinal file, generating its Normalized Mean_Cov file and then concating that to the NormalizedMeanCovByPositionFinal file..."

if [ -s $RunFolderPath/$JobName/$JobName.$PanelName.MeanCovByPosition ];
then
	sed -i 's/_mean_cvg//g' $RunFolderPath/$JobName/$JobName.$PanelName.MeanCovByPosition
	cat $RunFolderPath/$JobName/$JobName.$PanelName.MeanCovByPosition|sed 1d >> $RunFolderPath/Run_$PanelName"_MeanCovByPositionFinal".txt
	sync_to_demux $RunFolderPath/Run_$PanelName"_MeanCovByPositionFinal".txt

	if [ -s $RunFolderPath/$JobName/$JobName.StatSummary ];
	then
        	python $PipeScriptsFolderPath/CombineStats.py $RunFolderPath/$JobName/$JobName.$PanelName.MeanCovByPosition $RunFolderPath/$JobName/$JobName.StatSummary $RunFolderPath/$JobName/$JobName.$PanelName.NormalizedMeanCovByPosition
		if [ -s $RunFolderPath/$JobName/$JobName.$PanelName.NormalizedMeanCovByPosition ];
		then
			cat $RunFolderPath/$JobName/$JobName.$PanelName.NormalizedMeanCovByPosition|sed 1d >> $RunFolderPath/Run_$PanelName"_NormalizedMeanCovByPositionFinal".txt
			sync_to_demux $RunFolderPath/Run_$PanelName"_NormalizedMeanCovByPositionFinal".txt
		else
			echo "$JobName/$JobName.$PanelName.NormalizedMeanCovByPosition is either not found or is empty. It did not get concated to the Run_$PanelName'_NormalizedMeanCovByPositionFinal'.txt file"
		fi
		
		if [ "$PanelName" == "HEME_V1.2_HEMEao_V1.3" ] && [ -n "$HemeJobName" ] && [ -n "$Heme_aoJobName" ] && [ -s $RunFolderPath/Run_$PanelName"_NormalizedMeanCovByPositionFinal".txt ];
		then
			echo "$JobName is replaced with $Heme_aoJobName in all lines that contain $JobName in the Run_$PanelName'_NormalizedMeanCovByPositionFinal'.txt"
			sed -i "s/$JobName/$Heme_aoJobName/g" $RunFolderPath/Run_$PanelName"_NormalizedMeanCovByPositionFinal".txt
			sync_to_demux $RunFolderPath/Run_$PanelName"_NormalizedMeanCovByPositionFinal".txt
		fi
	else
        	echo "Cannot create $JobName/$JobName.$PanelName.NormalizedMeanCovByPosition.This is because $JobName/$JobName.StatSummary is either not found or is empty"
	fi
else
	echo "$JobName/$JobName.$PanelName.MeanCovByPosition is either not found or is empty. It did not get concated to the Run_$PanelName'_MeanCovByPositionFinal'.txt file"
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
