##V2.2.1
#!/bin/sh
## Run this script as sh New_Wrapper_Final.sh <RunFolderPath> <PipeScriptsFolderPath>
PipeScriptsFolderPath=$2
RunFolderPath=$1
for i in `egrep -v '_skip|-skip|HEME_V1.2|HEMEao_V1.3' $RunFolderPath/SampleSheet_ToBeProcessed.txt| awk '{print $1","$4}'`
do
SampleName=`echo $i|cut -d "," -f1`
PanelName=`echo $i|cut -d "," -f2`
if [ -s $RunFolderPath/$SampleName/$SampleName.PipeScript ];
then
	sh CreateFinal.sh $RunFolderPath $SampleName $PanelName $PipeScriptsFolderPath > $RunFolderPath/$SampleName/$SampleName.CreateFinal.out 2>&1
fi
done


if [ -s $RunFolderPath/Run_HEME_SamplePairsInfo.txt ];
then
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
		sh CreateFinal.sh $RunFolderPath $HEME_merge_sample,$HEME_sample,$HEMEao_sample HEME_V1.2_HEMEao_V1.3 $PipeScriptsFolderPath > $RunFolderPath/$HEME_merge_sample/$HEME_merge_sample.CreateFinal.out 2>&1
	fi
        done
else
	echo "$RunFolderPath/Run_HEME_SamplePairsInfo.txt NOT found, therefore no HEME_merge_sample to process..."
fi
