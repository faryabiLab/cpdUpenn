##V1.3.1
##Run this script as ./PipeParallel.sh <RunFolderName> <PipelineVersion>
#!/bin/sh

if [ $# -ne 2 ];
then
	echo "Usage $0 <RunFolderName> <PipelineVersion>"
	exit 2
else
	userInput=$1          ## Run Folder name
	PipeVersion=$2        ## Pipeline Version
	## Defining the Run folder on Drive_D
	localFolder=/Drive_D/illumina/Analysis/$userInput
	
	echo "Running WaitForSheet.sh.."

	## Waiting Till SampleSheet_Tobeprocessed_PG11.txt is found

	./WaitForSheet.sh $localFolder
	
	echo "Entering the Run folder on Drive_D"
	
	cd $localFolder

	##Copying the SampleSheet_Tobeprocessed_PG11.txt to SampleSheet_Tobeprocessed.txt (which would be used in the pipe-run).
	
	echo "Copying the SampleSheet_Tobeprocessed_PG11.txt to SampleSheet_Tobeprocessed.txt"
	
	cp SampleSheet_Tobeprocessed_PG11.txt SampleSheet_Tobeprocessed.txt
	
	
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
			mkdir /MSdata/illumina/Analysis/$userInput/Pipeline$PipeVersion/PG11
			echo "Copying the Data from Run Folder on Drive D to /MSdata/illumina/Analysis/'Run Folder'/Pipeline$PipeVersion/PG11 .."
			for i in `ls $localFolder `
			do
			        if [ -d $localFolder/$i ];
			        then
			        	cp -r $localFolder/$i /MSdata/illumina/Analysis/$userInput/Pipeline$PipeVersion/PG11
			        else
			        	cp $localFolder/$i /MSdata/illumina/Analysis/$userInput/Pipeline$PipeVersion/PG11
			        fi
			done
			exit 0
		else
			echo "Run Folder and Pipe Folder not created yet, hence waiting for 10 minutes.."
			sleep 10m
		fi
	done
fi


