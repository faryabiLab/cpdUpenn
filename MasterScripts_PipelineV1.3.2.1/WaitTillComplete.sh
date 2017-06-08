##V1.0
##Run this script as ./WaitTillComplete.sh <Path of RunFolder on MiseqAnalysis>

#!/bin/sh


name=$1
while :
do
	if [ -d $name ];
	then
		echo "Run Folder created\nSearching for AnalysisLog.txt..."

		while :
		do
			if [ -r $name/AnalysisLog.txt ];
			then
				echo "AnalysisLog.txt created\nWaiting for analysis to be completed..."
				while :
				do
					grepOut=`grep 'FASTQ generation time' $name/AnalysisLog.txt |wc -l`
	
					if [ $grepOut -eq 0 ];
					then
						echo "Analysis not completed yet, waiting for 10 more minutes..."
						sleep 10m
					elif [ $grepOut -eq 1 ];
					then
						echo "Analysis Completed\nMoving Forward with Copying Sample Fastqs and SampleSheet to the Run Folder located on /MSdata/illumina/MiSeqOutput  ..."
						exit 0 
					fi
				done
		
			else
				echo "AnalysisLog.txt is not created yet, waiting for 10 more minute.."
				sleep 10m
			fi
		done
		

	else
		echo "Run Folder not created yet\nWaiting for 1 hour..."
		sleep 1h
	fi

done




