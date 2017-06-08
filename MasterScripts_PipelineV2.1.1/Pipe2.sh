##V1.4
##Run this script as ./Pipe.sh <RunFolderName> <PipelineVersion> <Comp1> <Comp2>
#!/bin/sh
##Changes in V1.3.1: 1) Added Functionality to split the SampleSheet_Tobeprocessed.txt into two halves. First half would be processed on Path-Genomix-21 (PG21) comp and another half on Path-Genomix-11 (PG11). The step of conversion of SampleSheet to
##                      SampleSheet_Tobeprocessed.txt has been moved to this script and would be removed from MappingToStats.sh.
##                   2) The variables 'userInput' and 'PipeVersion' are now initialized inside the 'else' condition as it makes more sense to intialize them in that place.
##Changes in V1.3.2.1: 1) Introduced two new argument <Comp1> and <Comp2> that indicates two computers on which pipe is supposed to be split.
##Changes in V1.4 : 1) The script now makes SampleSheet_forMergeAnalysis.txt file for both the computers to process the HEME_V2.3 data on them.
##                  2) Echo statement indicating the creation of SampleSheet_Tobeprocessed.txt is modified to add the words "i.e. Making SampleSheet_Tobeprocessed.txt".


if [ $# -ne 4 ];
then
	echo "Usage $0 <RunFolderName> <PipelineVersion> <Comp1> <Comp2>"
	exit 2
else
	userInput=$1          ## Run Folder name
	PipeVersion=$2        ## Pipeline Version , example: V1.2
	Comp1=$3
	Comp2=$4
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
     		echo "Usage $0 <RunFolderName> <PipelineVersion> <Comp1> <Comp2>"
     		exit 2
	fi

	## Defining the Run folders on each drive
	MiSeqAnalysisFolder=/$MiSeq/$userInput
	localFolder_Comp1=/Drive_D/illumina/Analysis/$userInput
	localFolder_Comp2=/Drive_D_illumina_$Comp2/Analysis/$userInput
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

	mkdir $localFolder_Comp1
	cp $MSDataFolder/Data/Intensities/BaseCalls/*.fastq.gz $localFolder_Comp1
	cp $MSDataFolder/Data/Intensities/BaseCalls/SampleSheet.csv $localFolder_Comp1

	echo "Copying Complete..."

	## Copying the scripts to the folder on Drive D

	echo "Now Copying the Scripts from /MSdata/Informatics/Pipeline/CurrentPipes/Illumina Folder to the Run Folder on Local Drive D.."

	for i in `ls $ScriptsLoc`
	do
		if [ -d $ScriptsLoc/$i ];
		then
			cp -r $ScriptsLoc/$i $localFolder_Comp1
		else
			cp $ScriptsLoc/$i $localFolder_Comp1
		fi
	done

		


	echo "Copying Complete..."

	echo "Making Run Folder inside /Drive_D_illumina_$Comp2/Analysis  .."
	
	mkdir $localFolder_Comp2
	
	## Copying the scripts to the folder on Drive_D_illumina_Comp2/Analysis
	
	#echo "Now Copying the Scripts from /MSdata/Informatics/Pipeline/CurrentPipes/Illumina Folder to the Run Folder on Drive_D_illumina_$Comp2/Analysis.."
	
	#for i in `ls $ScriptsLoc`
	#do
	#	if [ -d $ScriptsLoc/$i ];
	#	then
	#		cp -r $ScriptsLoc/$i $localFolder_Comp2
	#	else
	#		cp $ScriptsLoc/$i $localFolder_Comp2
	#	fi
	#done
	
	
	#echo "Copying Complete..."
	
	#echo "Now Copying the Data(fastqs) from Run Folder on /MSdata/illumina/MiSeqOutput to the Run Folder which would be created on Drive_D_illumina_$Comp2/Analysis.."
	
	#cp $MSDataFolder/Data/Intensities/BaseCalls/*.fastq.gz $localFolder_Comp2

	
	#echo "Copying Complete..."
	
	
	
	
	echo "Entering the Run Folder located on the local Drive D..."
	
	cd $localFolder_Comp1
	
	## Modifying the SampleSheet to replace '_' with '-' in the SampleNames

	echo "Modifying the SampleSheet to replace '_' with '-' in SampleNames..."

	sed -i 's/_/-/g' SampleSheet.csv
	sed -i 's/-V1/_V1/g' SampleSheet.csv
	
	echo "Processing SampleSheet.csv to extract the SampleName,Barcodes and Description ...\n"
	sed -n '/Data/,//p' SampleSheet.csv | sed 1,2d | awk -F "," '{print $2,$6,$8,$10}' > Samplesheet_preFinal.txt
	echo "Making SampleSheet_forMergeAnalysis_preFinal.txt from SampleSheet.csv...\n"
	sed -n '/Data/,//p' SampleSheet.csv | sed 1,2d | awk -F "," '{print $1,$2,$6,$8,$10}' > SampleSheet_forMergeAnalysis_preFinal.txt
	
	
	echo "Modifying the SampleSheet_preFinal.txt to skip first half of the sample-list, i.e. Making SampleSheet_Tobeprocessed.txt"
	
	NumberOfSamples=`wc -l SampleSheet_preFinal.txt|cut -d " " -f1`
	HalfNumberOfSamples=$((NumberOfSamples / 2))
	RemainingHalfNumberOfSamples=$((NumberOfSamples-HalfNumberOfSamples))
	awk -v varHNS=$HalfNumberOfSamples 'BEGIN{count=1} {if(count<=varHNS){print $1,$2,$3,$4"_skip"}else{print}} {count++;}' Samplesheet_preFinal.txt > SampleSheet_Tobeprocessed.txt
	
	echo "Modifying the SampleSheet_forMergeAnalysis_preFinal.txt to completely remove the first half of the sample-list, i.e. Making SampleSheet_forMergeAnalysis.txt"
	
	awk -v varHNS=$HalfNumberOfSamples 'BEGIN{count=1} {if(count>varHNS){print}} {count++;}' SampleSheet_forMergeAnalysis_preFinal.txt > SampleSheet_forMergeAnalysis.txt
	
	echo "Modifying the SampleSheet_preFinal.txt to skip the last half of the sample-list ,i.e. Making SampleSheet_Tobeprocessed_$Comp2.txt"
	
	awk -v varHNS=$HalfNumberOfSamples 'BEGIN{count=1} {if(count>varHNS){print $1,$2,$3,$4"_skip"}else{print}} {count++;}' Samplesheet_preFinal.txt > SampleSheet_Tobeprocessed_$Comp2.txt
	
	echo "Modifying the SampleSheet_forMergeAnalysis_preFinal.txt to completely remove the last half of the sample-list ,i.e. Making SampleSheet_forMergeAnalysis_$Comp2.txt"
	
	awk -v varHNS=$HalfNumberOfSamples 'BEGIN{count=1} {if(count<=varHNS){print}} {count++;}' SampleSheet_forMergeAnalysis_preFinal.txt > SampleSheet_forMergeAnalysis_$Comp2.txt
	
	echo "Copying the SampleSheet_Tobeprocessed_$Comp2.txt to localfolder on Drive_D_$Comp2..."
	cp SampleSheet_Tobeprocessed_$Comp2.txt $localFolder_Comp2
	echo "Copying the SampleSheet_forMergeAnalysis_$Comp2.txt to localfolder on Drive_D_$Comp2..."
	cp SampleSheet_forMergeAnalysis_$Comp2.txt $localFolder_Comp2
	
	
	
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
	mkdir /MSdata/illumina/Analysis/$userInput/Pipeline$PipeVersion/$Comp1
	echo "Copying the Data from Run Folder on Drive D to /MSdata/illumina/Analysis/'Run Folder'/Pipeline$PipeVersion/$Comp1.."
	
	for i in `ls $localFolder_Comp1`
	do
	        if [ -d $localFolder_Comp1/$i ];
                then
                        cp -r $localFolder_Comp1/$i /MSdata/illumina/Analysis/$userInput/Pipeline$PipeVersion/$Comp1
                else
                        cp $localFolder_Comp1/$i /MSdata/illumina/Analysis/$userInput/Pipeline$PipeVersion/$Comp1
                fi
	
	done
			


fi


