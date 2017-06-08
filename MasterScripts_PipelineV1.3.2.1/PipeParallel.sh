##V1.3.2.1
##Run this script on Comp2 as ./PipeParallel.sh <RunFolderName> <PipelineVersion> <Comp2>
#!/bin/sh
##Changes in V1.3.2.1: 1)Added the <Comp2> argument to the command line, made related changes to the script below. Replaced PG11 with $Comp2 in the script below.

if [ $# -ne 3 ];
then
	echo "Usage $0 <RunFolderName> <PipelineVersion> <Comp2>"
	exit 2
else
	userInput=$1          ## Run Folder name
	PipeVersion=$2        ## Pipeline Version
	Comp2=$3
	## Defining the Run folder on Drive_D
	localFolder=/Drive_D/illumina/Analysis/$userInput
	
	echo "Running WaitForSheet.sh.."

	## Waiting Till SampleSheet_Tobeprocessed_Comp2.txt is found

	./WaitForSheet.sh $localFolder $Comp2
	
	echo "Entering the Run folder on Drive_D"
	
	cd $localFolder

	##Copying the SampleSheet_Tobeprocessed_Comp2.txt to SampleSheet_Tobeprocessed.txt (which would be used in the pipe-run).
	
	echo "Copying the SampleSheet_Tobeprocessed_$Comp2.txt to SampleSheet_Tobeprocessed.txt"
	
	cp SampleSheet_Tobeprocessed_$Comp2.txt SampleSheet_Tobeprocessed.txt
	
	
	## Running the Pipe

	echo "Running the Pipe..."

	./MappingToStats.sh

	## Copying the Run Folder from Drive D to /MSdata/illumina/Analysis
	
	cd $HOME
	echo "Checking whether the Run Folder and Pipe Folder (within that Run Folder) are already created by Pipe.sh in /MSdata/illumina/Analysis ..."
	
	while :
	do
		if [ -d /MSdata/illumina/Analysis/$userInput/Pipeline$PipeVersion ];
		then
			echo "Yes..the Pipe Folder has been detected in the Run Folder inside /MSdata/illumina/Analysis..."
			echo "Creating the Comp Folder inside the Pipe Folder .."
			mkdir /MSdata/illumina/Analysis/$userInput/Pipeline$PipeVersion/$Comp2
			echo "Copying the Data from Run Folder on Drive D to /MSdata/illumina/Analysis/'Run Folder'/Pipeline$PipeVersion/$Comp2 .."
			for i in `ls $localFolder `
			do
			        if [ -d $localFolder/$i ];
			        then
			        	cp -r $localFolder/$i /MSdata/illumina/Analysis/$userInput/Pipeline$PipeVersion/$Comp2
			        else
			        	cp $localFolder/$i /MSdata/illumina/Analysis/$userInput/Pipeline$PipeVersion/$Comp2
				fi
			done
			exit 0
		else
			echo "Run Folder and Pipe Folder not created yet, hence waiting for 10 minutes.."
			sleep 10m
		fi
	done
fi


