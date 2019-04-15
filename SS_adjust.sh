#!/bin/bash

SeqRun=$1
usr=$2
seq_id=$(awk NR==8 /PathCPD/illumina/SeqOut/${SeqRun}/First_Base_Report.htm | awk 'BEGIN {FS="_"} ; {print $3}')
flowcell2=$(awk NR==8 /PathCPD/illumina/SeqOut/${SeqRun}/First_Base_Report.htm | awk 'BEGIN {FS="_"} ; {print $4}')
flowcell=${flowcell2%??????}

cp /home/cpdlab/SampleSheet_base.csv /data1/HiSeqRun/Lymphoma/AmpliconDS/SampleSheet_base.csv
sed -i "s/180475/${seq_id}/g" SampleSheet_base.csv

cp /PathCPD/illumina/SeqOut/${SeqRun}/SampleSheet.csv /data1/HiSeqRun/Lymphoma/AmpliconDS/SampleSheet_demux.csv
col1=$(cat SampleSheet_demux.csv | awk 'BEGIN { FS="," } ; {print $2}' | tail -n +20) # sort uniq to cut in half so only lane 1 is there?
length=$(echo "${col1}" | wc -l)
#echo "${col1}" >> SampleSheet_lowerhalf.csv
col2=$(cat SampleSheet_demux.csv | awk 'BEGIN { FS="," } ; {print $3}' | tail -n +20)
col3=$(cat SampleSheet_demux.csv | awk 'BEGIN { FS="," } ; {print $10}' | tail -n +20)
NL=$'\n';i=1;col4=""
while [ $i -le $length ]; do col4+="/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta"${NL}; i=$(( $i + 1 )); done
j=1;col5=""
while [ $j -le $length ]; do col5+="1"${NL}; j=$(( $j + 1 )); done
k=1;col6=""
while [ $k -le $length ]; do col6+="/data/AmpliconDS"${NL}; k=$(( $k + 1 )); done
col7=$(cat SampleSheet_demux.csv | awk 'BEGIN { FS="," } ; {print $1}' | tail -n +20)
m=1;col8=""
while [ $m -le $length ]; do col8+=${flowcell}","${NL}; m=$(( $m + 1 )); done

paste -d"," <(printf %s "$col1") <(printf %s "$col2") <(printf %s "$col3") <(printf %s "$col4") <(printf %s "$col5") <(printf %s "$col6") <(printf %s "$col7") <(printf %s "$col8") > SampleSheet_meat.csv
cat SampleSheet_base.csv > SampleSheet_final.csv
awk 'BEGIN { FS="," } ; $7!=2' SampleSheet_meat.csv >> SampleSheet_final.csv

mv SampleSheet_final.csv SampleSheet.csv
