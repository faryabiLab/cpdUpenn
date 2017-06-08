##V1.3.2.1
##Run this script as ./WaitForSheet.sh <localRunFolderOnComp2> <Comp2>
##Changes in V1.3.2.1: 1) Changed the <localRunFolderOnPG11> to <localRunFolderOnComp2> in the "Run this script comment".
##                     2) Added 'Comp2' as a second argument to the command line.
##		       2) Changed SampleSheet_Tobeprocessed_PG11.txt to SampleSheet_Tobeprocessed_$Comp2.txt in the script.	

#!/bin/sh


name=$1
Comp2=$2
while :
do
	if [ -d $name ];
	then
		echo "Run Folder created\nSearching for Samplesheet_Tobeprocessed_$Comp2.txt..."

		while :
		do
			if [ -s $name/SampleSheet_Tobeprocessed_$Comp2.txt ];
			then
				echo "SampleSheet_Tobeprocessed_$Comp2.txt found!!!\nStarting the pipe on $Comp2"
				exit 0
			else
				echo "SampleSheet_Tobeprocessed_$Comp2.txt is not created yet, waiting for 15 more minute.."
				sleep 15m
			fi
		done
		

	else
		echo "Run Folder not created yet\nWaiting for 10 minutes..."
		sleep 10m
	fi

done




