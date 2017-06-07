#V2.1.1_hpc
## Run this script as ./ConcordanceExtract.sh <SampleName> <SampleCount> <Path_Scripts>
## Changes in PipelineV2.1.1_hpc: 1) Added the third argument, i.e. "Path for Scripts", on the command. This argument is used to point to '$Path_Scripts/FastqToTableConverter.py' script.
##                          2) Changed the argument <SampleName> to <PathOfSample> on command line to manage detection of Sample-Folder by the script. Also, added $PathOfSample variable to which the 'PathOfSample' gets assigned to. Modified the definition of 'NameOfSample' variable accordingly.
#!/bin/sh

PathOfSample=$1
NameOfSample=`echo $PathOfSample|cut -d "/" -f6`
CountOfSample=$2
R1_btrim_sum=$PathOfSample/$NameOfSample"_R1.btrim.summary.txt"
R2_btrim_sum=$PathOfSample/$NameOfSample"_R2.btrim.summary.txt"
Path_Scripts=$3

if [ -s $R1_btrim_sum ] && [ -s $R2_btrim_sum ]
then
	grep 'Pass' $R1_btrim_sum | sort -k 1,1 > $PathOfSample/$NameOfSample"_R1.btrim.summary.pass.sort.txt"
	grep 'Pass' $R2_btrim_sum | sort -k 1,1 > $PathOfSample/$NameOfSample"_R2.btrim.summary.pass.sort.txt"
	join $PathOfSample/$NameOfSample"_R1.btrim.summary.pass.sort.txt" $PathOfSample/$NameOfSample"_R2.btrim.summary.pass.sort.txt" -o 1.1,1.4,2.1,2.4 | awk '{if($2 == $4){print "@"$1}}' > $PathOfSample/$NameOfSample.concordantReadNames
	
	echo "Making btrimmed concordant fastqs"
	python $Path_Scripts/FastqToTableConverter.py $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.btrim.fastq.pe" $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.btrim.fastq.pe.tab"
	python $Path_Scripts/FastqToTableConverter.py $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.btrim.fastq.pe" $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.btrim.fastq.pe.tab"
	sort -k 1,1 $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.btrim.fastq.pe.tab" > $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.btrim.fastq.pe.tab.sort"
	sort -k 1,1 $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.btrim.fastq.pe.tab" > $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.btrim.fastq.pe.tab.sort"
	join -1 1 -2 1 $PathOfSample/$NameOfSample.concordantReadNames $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.btrim.fastq.pe.tab.sort" > $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.btrim.fastq.pe.tab.sort.cndt"
	join -1 1 -2 1 $PathOfSample/$NameOfSample.concordantReadNames $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.btrim.fastq.pe.tab.sort" > $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.btrim.fastq.pe.tab.sort.cndt"
	awk '{print $1" "$2"\n"$3"\n"$4" "$5"\n"$6}' $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.btrim.fastq.pe.tab.sort.cndt" > $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.btrim.fastq.pe.tab.sort.cndt.fastq"
	awk '{print $1" "$2"\n"$3"\n"$4" "$5"\n"$6}' $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.btrim.fastq.pe.tab.sort.cndt" > $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.btrim.fastq.pe.tab.sort.cndt.fastq"
	
	echo "Making raw concordant fastqs"
	python $Path_Scripts/FastqToTableConverter.py $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.cutAd.fastq" $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.cutAd.fastq.tab"
	python $Path_Scripts/FastqToTableConverter.py $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.cutAd.fastq" $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.cutAd.fastq.tab"
	sort -k 1,1 $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.cutAd.fastq.tab" > $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.cutAd.fastq.tab.sort"
	sort -k 1,1 $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.cutAd.fastq.tab" > $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.cutAd.fastq.tab.sort"
	join -1 1 -2 1 $PathOfSample/$NameOfSample.concordantReadNames $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.cutAd.fastq.tab.sort" > $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.cutAd.fastq.tab.sort.cndt"
	join -1 1 -2 1 $PathOfSample/$NameOfSample.concordantReadNames $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.cutAd.fastq.tab.sort" > $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.cutAd.fastq.tab.sort.cndt"
	awk '{print $1" "$2"\n"$3"\n"$4"\n"$5}' $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.cutAd.fastq.tab.sort.cndt" > $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.cutAd.fastq.tab.sort.cndt.fastq"
	awk '{print $1" "$2"\n"$3"\n"$4"\n"$5}' $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.cutAd.fastq.tab.sort.cndt" > $PathOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.cutAd.fastq.tab.sort.cndt.fastq"
else
	echo "Btrim summary files not found"
	exit 2;
fi

