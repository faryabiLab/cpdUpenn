##V1.2.4
#!/bin/sh
##Run this script as ./Pipe-DualRun.sh <RunFolderName> <PipelineVersion>
## Changes in V1.2.4 : a) A sed statement to replace '-V2' to '_V2' is added as one of the step while modifying the SampleSheet.csv so as to handle the 'HEME_V2.0' string in the Description column of the SampleSheet correctly.
##                     b) echo statements are changed to indicate the copying of fastqs directly from the MiseqAnalysisFolder to the localFolder.  

userInput=$1          ## Run Folder name
PipeVersion=$2        ## Pipeline Version , example: V1.2

if [ $# -ne 2 ];
then
	echo "Usage $0 <RunFolderName> <PipelineVersion>"
	exit 2
else
	MachineString=`echo $userInput | cut -d '_' -f2 `   ## detecting the miseq machine type from the folder name [i.e. $userInput ]
	if [ $MachineString = M00351 ];
	then
     		MiSeq=MiSeqAnalysis
	elif [ $MachineString = M01651 ];
	then
     		MiSeq=MiSeqAnalysis2
	else
     		echo "Error in RunFolderName. Please mention the correct name (Example:130123_M00351_0016_000000000-A1T7P) and try again..."
     		echo "Usage $0 <RunFolderName> <PipelineVersion>"
     		exit 2
	fi

	## Defining the Run folders on each drive
	MiSeqAnalysisFolder=/$MiSeq/$userInput
	localFolder=/Drive_D/illumina/Analysis/$userInput
	MSDataFolder=/MSdata/illumina/MiSeqOutput/$userInput

	## Defining the Location of the Pipeline Scripts of the Pipeline Version mentioned on command line

	ScriptsLoc=/MSdata/Informatics/Pipeline/CurrentPipes/Illumina/Pipeline$PipeVersion

	echo "Running WaitTillComplete.sh.."

	## Waiting Till MiSeq Analysis gets completed

	./WaitTillComplete.sh $MiSeqAnalysisFolder

	## Copying MiSeqAnalysis Sample Fastqs and SampleSheet to the run folder on /MSdata/illumina/MiSeqOutput


	#cp $MiSeqAnalysisFolder/Data/Intensities/BaseCalls/*.fastq.gz $MSDataFolder/Data/Intensities/BaseCalls
	#cp $MiSeqAnalysisFolder/SampleSheet.csv $MSDataFolder/Data/Intensities/BaseCalls
 

	echo "Fastqs and SampleSheet were not copied from $MiSeqAnalysisFolder to the run folder on /MSdata/illumina/MiSeqOutput"

	## Copying Data from the Run Folder on $MiSeq to the Run Folder on Local Drive D

	echo "Now Copying the Data from Run Folder on $MiSeq to the Run Folder which would be created on Local Drive D.."
	
 	echo "Making Run Folder inside /Drive_D/illumina/Analysis  .."

	mkdir $localFolder
	cp $MiSeqAnalysisFolder/Data/Intensities/BaseCalls/*.fastq.gz $localFolder
	cp $MiSeqAnalysisFolder/SampleSheet.csv $localFolder

	echo "Copying Complete..."

	## Copying the scripts to the folder on Drive D

	echo "Now Copying the Scripts from /MSdata/Informatics/Pipeline/CurrentPipes/Illumina Folder to the Run Folder on Local Drive D.."

	for i in `ls $ScriptsLoc`
	do
		if [ -d $ScriptsLoc/$i ];
		then
			cp -r $ScriptsLoc/$i $localFolder
		else
			cp $ScriptsLoc/$i $localFolder
		fi
	done

		


	echo "Copying Complete..."

	cd $localFolder

	## Modifying the SampleSheet to replace '_' with '-' in the SampleNames

	echo "Modifying the SampleSheet to replace '_' with '-' in SampleNames..."

	sed -i 's/_/-/g' SampleSheet.csv
	sed -i 's/-V1/_V1/g' SampleSheet.csv
	sed -i 's/-V2/_V2/g' SampleSheet.csv

	## Run the Pipe

	echo "Running the Pipe..."

	./MappingToStats.sh

	## Copying the Run Folder from Drive D to /MSdata/illumina/Analysis
	
	#cd $HOME
	#echo "Creating the Run Folder in /MSdata/illumina/Analysis .."
	#mkdir /MSdata/illumina/Analysis/$userInput
	#echo "Creating the Pipe Folder in the same Run Folder .."
	#mkdir /MSdata/illumina/Analysis/$userInput/Pipeline$PipeVersion
	#echo "Copying the Data from Run Folder on Drive D to /MSdata/illumina/Analysis/'Run Folder'/Pipeline$PipeVersion/ .."
	
	#for i in `ls $localFolder `
	#do
	       # if [ -d $localFolder/$i ];
               # then
               #         cp -r $localFolder/$i /MSdata/illumina/Analysis/$userInput/Pipeline$PipeVersion
               # else
                        #cp $localFolder/$i /MSdata/illumina/Analysis/$userInput/Pipeline$PipeVersion
                #fi
	
	#done
			


fi


