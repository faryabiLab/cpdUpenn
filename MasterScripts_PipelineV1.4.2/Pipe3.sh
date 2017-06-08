##V1.4.2
##Run this script as ./Pipe3.sh <RunFolderName> <PipelineVersion> <Comp1> <Comp2>
#!/bin/sh
## Pipe3.sh script is introduced to eliminate NAS from the workflow. It adopts the new storage system named Isilon.
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
     		MiSeqAnalysis=MiSeqAnalysis
		MiSeqOutput=MiSeqOutput
	elif [ $MachineString = M01651 ];
	then
     		MiSeqAnalysis=MiSeqAnalysis2
		MiSeqOutput=MiSeqOutput2
     	elif [ $MachineString = M03333 ];
     	then
     		MiSeqAnalysis=MiSeqAnalysis3
		MiSeqOutput=MiSeqOutput3
	else
     		echo "Error in RunFolderName. Please mention the correct name (Example:130123_M00351_0016_000000000-A1T7P) and try again..."
     		echo "Usage $0 <RunFolderName> <PipelineVersion> <Comp1> <Comp2>"
     		exit 2
	fi

	## Defining the Run folders on each drive
	MiSeqAnalysisFolder=/$MiSeqAnalysis/$userInput
	MiSeqOutputFolder=/$MiSeqOutput/$userInput
	localFolder_Comp1=/Drive_D/illumina/Analysis/$userInput
	localFolder_Comp2=/Drive_D_illumina_$Comp2/Analysis/$userInput
	MSoutLocation=/MSdata/illumina/MSout/

	## Defining the Location of the Pipeline Scripts of the Pipeline Version mentioned on command line

	ScriptsLoc=/MSdata/Informatics/Pipeline/CurrentPipes/Illumina/Pipeline$PipeVersion

	echo "Running WaitTillComplete.sh.."

	## Waiting Till MiSeq Analysis gets completed

	./WaitTillComplete.sh $MiSeqAnalysisFolder

	## Copying Sample Fastqs and SampleSheet from run folder on MiSeqAnalysis to that on MiSeqOutput


	##cp $MiSeqAnalysisFolder/Data/Intensities/BaseCalls/*.fastq.gz $MiSeqOutputFolder/Data/Intensities/BaseCalls
	cp $MiSeqAnalysisFolder/SampleSheet.csv $MiSeqOutputFolder/Data/Intensities/BaseCalls
 

	echo "Copying Complete.."

	## Copying Data from the Run Folder on MiSeqOutput to the Run Folder on Local Drive D

	echo "Now Copying the Data from Run Folder on MiSeqOutput to the Run Folder which would be created on Local Drive D.."
	
 	echo "Making Run Folder inside /Drive_D/illumina/Analysis  .."

	mkdir $localFolder_Comp1
	cp $MiSeqOutputFolder/Data/Intensities/BaseCalls/*.fastq.gz $localFolder_Comp1
	cp $MiSeqOutputFolder/Data/Intensities/BaseCalls/SampleSheet.csv $localFolder_Comp1

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

	echo "Copying the Run Folder on MiSeqOutput to MSoutLocation"
	cp -r $MiSeqOutputFolder $MSoutLocation
	
	
	
fi


