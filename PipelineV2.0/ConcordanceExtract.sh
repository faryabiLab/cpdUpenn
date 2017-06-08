#V1.3.3
## Run this script as ./ConcordanceExtract.sh <SampleName> <SampleCount>
#!/bin/sh

NameOfSample=$1
CountOfSample=$2
R1_btrim_sum=$NameOfSample/$NameOfSample"_R1.btrim.summary.txt"
R2_btrim_sum=$NameOfSample/$NameOfSample"_R2.btrim.summary.txt"

if [ -s $R1_btrim_sum ] && [ -s $R2_btrim_sum ]
then
	grep 'Pass' $R1_btrim_sum | sort -k 1,1 > $NameOfSample/$NameOfSample"_R1.btrim.summary.pass.sort.txt"
	grep 'Pass' $R2_btrim_sum | sort -k 1,1 > $NameOfSample/$NameOfSample"_R2.btrim.summary.pass.sort.txt"
	join $NameOfSample/$NameOfSample"_R1.btrim.summary.pass.sort.txt" $NameOfSample/$NameOfSample"_R2.btrim.summary.pass.sort.txt" -o 1.1,1.4,2.1,2.4 | awk '{if($2 == $4){print "@"$1}}' > $NameOfSample/$NameOfSample.concordantReadNames
	
	echo "Making btrimmed concordant fastqs"
	python FastqToTableConverter.py $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.btrim.fastq.pe" $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.btrim.fastq.pe.tab"
	python FastqToTableConverter.py $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.btrim.fastq.pe" $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.btrim.fastq.pe.tab"
	sort -k 1,1 $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.btrim.fastq.pe.tab" > $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.btrim.fastq.pe.tab.sort"
	sort -k 1,1 $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.btrim.fastq.pe.tab" > $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.btrim.fastq.pe.tab.sort"
	join -1 1 -2 1 $NameOfSample/$NameOfSample.concordantReadNames $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.btrim.fastq.pe.tab.sort" > $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.btrim.fastq.pe.tab.sort.cndt"
	join -1 1 -2 1 $NameOfSample/$NameOfSample.concordantReadNames $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.btrim.fastq.pe.tab.sort" > $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.btrim.fastq.pe.tab.sort.cndt"
	awk '{print $1" "$2"\n"$3"\n"$4" "$5"\n"$6}' $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.btrim.fastq.pe.tab.sort.cndt" > $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.btrim.fastq.pe.tab.sort.cndt.fastq"
	awk '{print $1" "$2"\n"$3"\n"$4" "$5"\n"$6}' $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.btrim.fastq.pe.tab.sort.cndt" > $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.btrim.fastq.pe.tab.sort.cndt.fastq"
	
	echo "Making raw concordant fastqs"
	python FastqToTableConverter.py $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.cutAd.fastq" $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.cutAd.fastq.tab"
	python FastqToTableConverter.py $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.cutAd.fastq" $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.cutAd.fastq.tab"
	sort -k 1,1 $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.cutAd.fastq.tab" > $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.cutAd.fastq.tab.sort"
	sort -k 1,1 $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.cutAd.fastq.tab" > $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.cutAd.fastq.tab.sort"
	join -1 1 -2 1 $NameOfSample/$NameOfSample.concordantReadNames $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.cutAd.fastq.tab.sort" > $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.cutAd.fastq.tab.sort.cndt"
	join -1 1 -2 1 $NameOfSample/$NameOfSample.concordantReadNames $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.cutAd.fastq.tab.sort" > $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.cutAd.fastq.tab.sort.cndt"
	awk '{print $1" "$2"\n"$3"\n"$4"\n"$5}' $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.cutAd.fastq.tab.sort.cndt" > $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R1_001.cutAd.fastq.tab.sort.cndt.fastq"
	awk '{print $1" "$2"\n"$3"\n"$4"\n"$5}' $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.cutAd.fastq.tab.sort.cndt" > $NameOfSample/$NameOfSample"_"$CountOfSample"_L001_R2_001.cutAd.fastq.tab.sort.cndt.fastq"
else
	echo "Btrim summary files not found"
	exit 2;
fi

