##V1.3.1
##Run this script as ./WaitForSheet.sh <localRunFolderOnPG11>

#!/bin/sh


name=$1
while :
do
	if [ -d $name ];
	then
		echo "Run Folder created\nSearching for Samplesheet_Tobeprocessed_PG11.txt..."

		while :
		do
			if [ -s $name/SampleSheet_Tobeprocessed_PG11.txt ];
			then
				echo "SampleSheet_Tobeprocessed_PG11.txt found!!!\nStarting the pipe on PG11"
				exit 0
			else
				echo "SampleSheet_Tobeprocessed_PG11.txt is not created yet, waiting for 15 more minute.."
				sleep 15m
			fi
		done
		

	else
		echo "Run Folder not created yet\nWaiting for 10 minutes..."
		sleep 10m
	fi

done




