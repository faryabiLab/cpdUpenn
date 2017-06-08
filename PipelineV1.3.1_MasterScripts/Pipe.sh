##V1.3.1
##Run this script as ./Pipe.sh <RunFolderName> <PipelineVersion>
#!/bin/sh
##Changes in V1.3.1: 1) Added Functionality to split the SampleSheet_Tobeprocessed.txt into two halves. First half would be processed on Path-Genomix-21 (PG21) comp and another half on Path-Genomix-11 (PG11). The step of conversion of SampleSheet to
##                      SampleSheet_Tobeprocessed.txt has been moved to this script and would be removed from MappingToStats.sh.
##                   2) The variables 'userInput' and 'PipeVersion' are now initialized inside the 'else' condition as it makes more sense to intialize them in that place.



if [ $# -ne 2 ];
then
	echo "Usage $0 <RunFolderName> <PipelineVersion>"
	exit 2
else
	userInput=$1          ## Run Folder name
	PipeVersion=$2        ## Pipeline Version , example: V1.2
	MachineString=`echo $userInput | cut -d '_' -f2 `   ## detecting the miseq machine type from the folder name [i.e. $userInput ]
	if [ $MachineString = M00351 ];
	then
     		MiSeq=MiSeqAnalysis
	elif [ $MachineString = M01651 ];
	then
     		MiSeq=MiSeqAnalysis2
     	elif [ $MachineString = M03333 ];
     	then
     		MiSeq=MiSeqAnalysis3
	else
     		echo "Error in RunFolderName. Please mention the correct name (Example:130123_M00351_0016_000000000-A1T7P) and try again..."
     		echo "Usage $0 <RunFolderName> <PipelineVersion>"
     		exit 2
	fi

	## Defining the Run folders on each drive
	MiSeqAnalysisFolder=/$MiSeq/$userInput
	localFolder=/Drive_D/illumina/Analysis/$userInput
	localFolder_PG11=/Drive_D_illumina_PG11/Analysis/$userInput
	MSDataFolder=/MSdata/illumina/MiSeqOutput/$userInput

	## Defining the Location of the Pipeline Scripts of the Pipeline Version mentioned on command line

	ScriptsLoc=/MSdata/Informatics/Pipeline/CurrentPipes/Illumina/Pipeline$PipeVersion

	echo "Running WaitTillComplete.sh.."

	## Waiting Till MiSeq Analysis gets completed

	./WaitTillComplete.sh $MiSeqAnalysisFolder

	## Copying MiSeqAnalysis Sample Fastqs and SampleSheet to the run folder on /MSdata/illumina/MiSeqOutput


	cp $MiSeqAnalysisFolder/Data/Intensities/BaseCalls/*.fastq.gz $MSDataFolder/Data/Intensities/BaseCalls
	cp $MiSeqAnalysisFolder/SampleSheet.csv $MSDataFolder/Data/Intensities/BaseCalls
 

	echo "Copying Complete.."

	## Copying Data from the Run Folder on /MSdata/illumina/MiSeqOutput to the Run Folder on Local Drive D

	echo "Now Copying the Data from Run Folder on /MSdata/illumina/MiSeqOutput to the Run Folder which would be created on Local Drive D.."
	
 	echo "Making Run Folder inside /Drive_D/illumina/Analysis  .."

	mkdir $localFolder
	cp $MSDataFolder/Data/Intensities/BaseCalls/*.fastq.gz $localFolder
	cp $MSDataFolder/Data/Intensities/BaseCalls/SampleSheet.csv $localFolder

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
	
	echo "Making Run Folder inside /Drive_D_illumina_PG11/Analysis  .."
	
	mkdir $localFolder_PG11
	
	## Copying the scripts to the folder on Drive_D_illumina_PG11/Analysis
	
	echo "Now Copying the Scripts from /MSdata/Informatics/Pipeline/CurrentPipes/Illumina Folder to the Run Folder on Drive_D_illumina_PG11/Analysis.."
	
	for i in `ls $ScriptsLoc`
	do
		if [ -d $ScriptsLoc/$i ];
		then
			cp -r $ScriptsLoc/$i $localFolder_PG11
		else
			cp $ScriptsLoc/$i $localFolder_PG11
		fi
	done
	
	
	echo "Copying Complete..."
	
	echo "Now Copying the Data(fastqs) from Run Folder on /MSdata/illumina/MiSeqOutput to the Run Folder which would be created on Drive_D_illumina_PG11/Analysis.."
	
	cp $MSDataFolder/Data/Intensities/BaseCalls/*.fastq.gz $localFolder_PG11

	
	echo "Copying Complete..."
	
	
	
	
	echo "Entering the Run Folder located on the local Drive D..."
	
	cd $localFolder
	
	## Modifying the SampleSheet to replace '_' with '-' in the SampleNames

	echo "Modifying the SampleSheet to replace '_' with '-' in SampleNames..."

	sed -i 's/_/-/g' SampleSheet.csv
	sed -i 's/-V1/_V1/g' SampleSheet.csv
	
	echo "Processing SampleSheet.csv to extract the SampleName,Barcodes and Description ...\n"
	sed -n '/Data/,//p' SampleSheet.csv | sed 1,2d | awk -F "," '{print $2,$6,$8,$10}' > Samplesheet_preFinal.txt
	
	
	echo "Modifying the SampleSheet_preFinal.txt to skip first half of the sample-list..."
	
	NumberOfSamples=`wc -l SampleSheet_preFinal.txt|cut -d " " -f1`
	HalfNumberOfSamples=$((NumberOfSamples / 2))
	RemainingHalfNumberOfSamples=$((NumberOfSamples-HalfNumberOfSamples))
	awk -v varHNS=$HalfNumberOfSamples 'BEGIN{count=1} {if(count<=varHNS){print $1,$2,$3,$4"_skip"}else{print}} {count++;}' Samplesheet_preFinal.txt > SampleSheet_Tobeprocessed.txt
	
	echo "Modifying the SampleSheet_preFinal.txt to skip the last half of the sample-list ,i.e. Making SampleSheet_Tobeprocessed_PG11.txt"
	
	awk -v varHNS=$HalfNumberOfSamples 'BEGIN{count=1} {if(count>varHNS){print $1,$2,$3,$4"_skip"}else{print}} {count++;}' Samplesheet_preFinal.txt > SampleSheet_Tobeprocessed_PG11.txt
	
	echo "Copying the SampleSheet_Tobeprocessed_PG11.txt to localfolder on Drive_D_PG11..."
	cp SampleSheet_Tobeprocessed_PG11.txt $localFolder_PG11
	
	
	## Run the Pipe
	echo "Running the Pipe..."
	
	./MappingToStats.sh

	## Copying the Run Folder from Drive D to /MSdata/illumina/Analysis
	
	cd $HOME
	echo "Creating the Run Folder in /MSdata/illumina/Analysis .."
	mkdir /MSdata/illumina/Analysis/$userInput
	echo "Creating the Pipe Folder in the same Run Folder .."
	mkdir /MSdata/illumina/Analysis/$userInput/Pipeline$PipeVersion
	echo "Creating the Comp Folder inside the Pipe Folder .."
	mkdir /MSdata/illumina/Analysis/$userInput/Pipeline$PipeVersion/PG21
	echo "Copying the Data from Run Folder on Drive D to /MSdata/illumina/Analysis/'Run Folder'/Pipeline$PipeVersion/PG21 .."
	
	for i in `ls $localFolder `
	do
	        if [ -d $localFolder/$i ];
                then
                        cp -r $localFolder/$i /MSdata/illumina/Analysis/$userInput/Pipeline$PipeVersion/PG21
                else
                        cp $localFolder/$i /MSdata/illumina/Analysis/$userInput/Pipeline$PipeVersion/PG21
                fi
	
	done
			


fi


