##V1.4alphaM
##Run this script as ./Merge.sh
#!/bin/sh

#Process SampleSheet.csv

sed -n '/Data/,//p' SampleSheet.csv | sed 1,2d | awk -F "," '{print $1,$2,$6,$8,$10}' > SampleSheet_forMergeAnalysis.txt

# Create a Pipeline Script for each Sub-Directory
echo "Creating the PipelineScript for each sample ...\n"
python Merge.py SampleSheet_forMergeAnalysis.txt

# Run the Pipe for each sample
for i in `ls|grep merged`
do
if [ -d $i ];
then
	if [ -s $i/$i.PipeScript ];
	then
		echo "Running the Pipe for Sample $i ...\n"
	
		nohup $i/$i.PipeScript > $i/$i.out
	
	
		echo "Preparing raw Run Files ... \n"
	
		if [ -e $i/$i.StatSummary ];
		then
			cat $i/$i.StatSummary >> RunStats.txt
		fi
	
		if [ -s $i.canon_annotated_flagged ];
		then
			cat $i.canon_annotated_flagged >> Run_masterVar.txt
		fi
		
		
		
		## Concating the Panel-specific MeanCovByPosition file of each Sample to the final Run_(Panel)_MeanCovByPosition.txt.
		
		if [ -e $i/$i.HEME_V1.2.MeanCovByPosition ];
		then
			cat $i/$i.HEME_V1.2.MeanCovByPosition >> Run_HEME_V1.2_MeanCovByPosition.txt
		elif [ -e $i/$i.HEME_V1.1.MeanCovByPosition ];
		then
			cat $i/$i.HEME_V1.1.MeanCovByPosition >> Run_HEME_V1.1_MeanCovByPosition.txt
		elif [ -e $i/$i.FFPE_V1.1.MeanCovByPosition ];
		then
			cat $i/$i.FFPE_V1.1.MeanCovByPosition >> Run_FFPE_V1.1_MeanCovByPosition.txt
		elif [ -e $i/$i.mPPP.MeanCovByPosition ];
		then
			cat $i/$i.mPPP.MeanCovByPosition >> Run_mPPP_MeanCovByPosition.txt
		elif [ -e $i/$i.bPPP.MeanCovByPosition ];
		then
			cat $i/$i.bPPP.MeanCovByPosition >> Run_bPPP_MeanCovByPosition.txt
		elif [ -e $i/$i.CEBPa.MeanCovByPosition ];
		then
			cat $i/$i.CEBPa.MeanCovByPosition >> Run_CEBPa_MeanCovByPosition.txt
		elif [ -e $i/$i.HEME_V1.2_HEMEao_V1.3.MeanCovByPosition ];
		then
			cat $i/$i.HEME_V1.2_HEMEao_V1.3.MeanCovByPosition >> Run_HEME_V1.2_HEMEao_V1.3_MeanCovByPosition.txt
		fi
	
	
	
		if [ -s $i/$i.Depth.below250X.amplicons_sorted.summary ];
		then
			cat $i/$i.Depth.below250X.amplicons_sorted.summary >> Run_AmpliconSummary.txt
		fi


	else
		echo "No PipeScript found for Sample $i, hence skipping it...\n"


	fi
fi
done




## Creating Final Run Files

echo "Creating Final Run Stats File..\n"
sort -k 1,1 -r RunStats.txt|uniq > RunStatsFinal.txt

echo "Creating Final MasterVar file..\n"
sort -k 1,1 -r Run_masterVar.txt|uniq > Run_masterVarFinal.txt

echo "Creating Final Amplicon Summary File..\n"
sort -k 1,1 -r Run_AmpliconSummary.txt|uniq > Run_AmpliconSummaryFinal.txt





echo "Creating Final MeanCovByPosition and NormalizedMeanCovByPosition Files for each Panel..\n"

if [ -e Run_HEME_V1.2_MeanCovByPosition.txt ];
then
	sort -k 1,1 -r Run_HEME_V1.2_MeanCovByPosition.txt|uniq > Run_HEME_V1.2_MeanCovByPositionFinal.txt	
	python CombineStats.py Run_HEME_V1.2_MeanCovByPositionFinal.txt RunStatsFinal.txt Run_HEME_V1.2_NormalizedMeanCovByPosition.txt
fi

if [ -e Run_HEME_V1.1_MeanCovByPosition.txt ];
then
	sort -k 1,1 -r Run_HEME_V1.1_MeanCovByPosition.txt|uniq > Run_HEME_V1.1_MeanCovByPositionFinal.txt
	python CombineStats.py Run_HEME_V1.1_MeanCovByPositionFinal.txt RunStatsFinal.txt Run_HEME_V1.1_NormalizedMeanCovByPosition.txt
fi

if [ -e Run_FFPE_V1.1_MeanCovByPosition.txt ];
then
	sort -k 1,1 -r Run_FFPE_V1.1_MeanCovByPosition.txt|uniq > Run_FFPE_V1.1_MeanCovByPositionFinal.txt
	python CombineStats.py Run_FFPE_V1.1_MeanCovByPositionFinal.txt RunStatsFinal.txt Run_FFPE_V1.1_NormalizedMeanCovByPosition.txt
fi
if [ -e Run_mPPP_MeanCovByPosition.txt ];
then
	sort -k 1,1 -r Run_mPPP_MeanCovByPosition.txt|uniq > Run_mPPP_MeanCovByPositionFinal.txt
	python CombineStats.py Run_mPPP_MeanCovByPositionFinal.txt RunStatsFinal.txt Run_mPPP_NormalizedMeanCovByPosition.txt
fi
if [ -e Run_bPPP_MeanCovByPosition.txt ];
then
	sort -k 1,1 -r Run_bPPP_MeanCovByPosition.txt|uniq > Run_bPPP_MeanCovByPositionFinal.txt
	python CombineStats.py Run_bPPP_MeanCovByPositionFinal.txt RunStatsFinal.txt Run_bPPP_NormalizedMeanCovByPosition.txt
fi
if [ -e Run_CEBPa_MeanCovByPosition.txt ];
then
	sort -k 1,1 -r Run_CEBPa_MeanCovByPosition.txt|uniq > Run_CEBPa_MeanCovByPositionFinal.txt
	python CombineStats.py Run_CEBPa_MeanCovByPositionFinal.txt RunStatsFinal.txt Run_CEBPa_NormalizedMeanCovByPosition.txt
fi
if [ -e Run_HEME_V1.2_HEMEao_V1.3_MeanCovByPosition.txt ];
then
	sort -k 1,1 -r Run_HEME_V1.2_HEMEao_V1.3_MeanCovByPosition.txt|uniq > Run_HEME_V1.2_HEMEao_V1.3_MeanCovByPositionFinal.txt	
	python CombineStats.py Run_HEME_V1.2_HEMEao_V1.3_MeanCovByPositionFinal.txt RunStatsFinal.txt Run_HEME_V1.2_HEMEao_V1.3_NormalizedMeanCovByPosition.txt
fi

