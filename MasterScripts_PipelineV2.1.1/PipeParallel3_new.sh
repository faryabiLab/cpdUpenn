##V2.1.1
##Run this script on Comp2 as ./PipeParallel3.sh <RunFolderName> <PipelineVersion> <Comp2>
#!/bin/sh
##PipeParallel3.sh eliminates NAS from the workflow. It adopts the new storage system ,i.e. Isilon
##Changes in V2.1.1 : 1) Introduced MiSeqAnalysis and MiSeqAnalysisFolder variables. 2) The script now copies the data from MiSeqAnalysis folder instead of (MiSeqOutput folder) to the local drive of Comp2.

if [ $# -ne 3 ];
then
	echo "Usage $0 <RunFolderName> <PipelineVersion> <Comp2>"
	exit 2
else
	userInput=$1          ## Run Folder name
	PipeVersion=$2        ## Pipeline Version
	Comp2=$3
	MachineString=`echo $userInput | cut -d '_' -f2 `   ## detecting the miseq machine type from the folder name [i.e. $userInput ]
	if [ $MachineString = M00351 ];
	then
		MiSeqOutput=MiSeqOutput
		MiSeqAnalysis=MiSeqAnalysis
	elif [ $MachineString = M01651 ];
	then
		MiSeqOutput=MiSeqOutput2
		MiSeqAnalysis=MiSeqAnalysis2
     	elif [ $MachineString = M03333 ];
     	then
		MiSeqOutput=MiSeqOutput3
		MiSeqAnalysis=MiSeqAnalysis3
	else
     		echo "Error in RunFolderName. Please mention the correct name (Example:130123_M00351_0016_000000000-A1T7P) and try again..."
		echo "Usage $0 <RunFolderName> <PipelineVersion> <Comp2>"
     		exit 2
	fi
	## Defining the Run folder on Drive_D
	localFolder=/Drive_D/illumina/Analysis/$userInput
	## Defining the MiSeqOutputFolder
	MiSeqOutputFolder=/$MiSeqOutput/$userInput
	MiSeqAnalysisFolder=/$MiSeqAnalysis/$userInput
	## Defining the Location of the Pipeline Scripts of the Pipeline Version mentioned on command line
	ScriptsLoc=/MSdata/Informatics/Pipeline/CurrentPipes/Illumina/Pipeline$PipeVersion

	echo "Running WaitForSheet.sh.."
	## Waiting Till SampleSheet_Tobeprocessed_Comp2.txt is found

        ./WaitForSheet.sh $localFolder $Comp2

	echo "Now Copying the Data from Run Folder on MiSeqAnalysisFolder to the Run Folder on Local Drive D.."
	cp $MiSeqAnalysisFolder/Data/Intensities/BaseCalls/*.fastq.gz $localFolder
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

	echo "Entering the Run folder on Drive_D"
	
	cd $localFolder

	##Copying the SampleSheet_Tobeprocessed_Comp2.txt to SampleSheet_Tobeprocessed.txt (which would be used in the pipe-run).
	
	echo "Copying the SampleSheet_Tobeprocessed_$Comp2.txt to SampleSheet_Tobeprocessed.txt"
	
	cp SampleSheet_Tobeprocessed_$Comp2.txt SampleSheet_Tobeprocessed.txt
	
	echo "Copying the SampleSheet_forMergeAnalysis_$Comp2.txt to SampleSheet_forMergeAnalysis.txt"
	
	cp SampleSheet_forMergeAnalysis_$Comp2.txt SampleSheet_forMergeAnalysis.txt
	
	
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


