#!/bin/bash

# Login as root to run the following:

SeqRun=$1
usr=$2
seq_id=$(awk NR==8 /PathCPD/illumina/SeqOut/${SeqRun}/First_Base_Report.htm | awk 'BEGIN {FS="_"} ; {print $3}')
flowcell2=$(awk NR==8 /PathCPD/illumina/SeqOut/${SeqRun}/First_Base_Report.htm | awk 'BEGIN {FS="_"} ; {print $4}')
flowcell=${flowcell2%??????}

# Demultiplexing step
nohup bcl2fastq --runfolder-dir /PathCPD/illumina/SeqOut/${SeqRun} --output-dir /data1/HiSeqRun/Lymphoma/AmpliconDS/Data/Intensities/BaseCalls/ --barcode-mismatches 1 --sample-sheet /PathCPD/illumina/SeqOut/${SeqRun}/SampleSheet_demux.csv --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0

# Create sample ID directories (A/B)
#cd /data1/HiSeqRun/Lymphoma/AmpliconDS/Analysis/Fastq

#for line in $(ls *.fastq.gz | cut -d'_' -f1 | sort -u)
#do
#  mkdir ${line}"A"
#  mkdir ${line}"B"
#done

# Move generated FastQs to their proper directory
#count=0
#declare -a array=(A B A B A B A B etc...)
#for line in $(ls *.fastq.gz | cut -d'_' -f1,2 | sort -t _ -u -k 1,1 -k 2.2n)
#do
#  sample=$(cut -d'_' -f1 <<< $line)
#  mv ${line}*.fastq.gz ${sample}${array[count]}
#  count=$((count+1))
#done

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

### MODIFY THE SAMPLESHEET TO SUIT THE DOCKER IMAGE FORMAT ###

# Window #1 | Start up Docker daemon on Isilon
#sudo su root
#cd /data1/var/lib/docker/
#rm -rf *
#cd /data1/HiSeqRun/Lymphoma
#systemctl start docker
#systemctl stop docker
#dockerd -g /data1/var/lib/docker

## Window #2 | Analysis step
cd /data1/HiSeqRun/Lymphoma
cat ampliconds-1-2-1.tar | docker load
docker run -v /data1/HiSeqRun/Lymphoma:/data -v /data1/HiSeqRun/Lymphoma/genomes:/genomes docker.illumina.com/arnoldliao/ampliconds-1-2-1 /opt/illumina/Isis/latest/Isis -r /data/AmpliconDS -a /data/AmpliconDS/Analysis -c 2

# Post-processing/move to HPC
cd /data1/HiSeqRun/Lymphoma/AmpliconDS/Analysis
mkdir genome_VCFs
mv *genome.vcf.gz genome_VCFs/
for file in *vcf.gz; do zcat $file > ${file%???}; done

rsync -a --include='*.vcf' --include='*AmpliconCoverage*' --include='*bam*' --exclude='*/' --exclude='*' . ${usr}@mercury.pmacs.upenn.edu:/project/cpdlab/Lymphoma/$SeqRun
cd ../Data/Intensities/BaseCalls/
rsync -a *.fastq.gz ${usr}@mercury.pmacs.upenn.edu:/project/cpdlab/Lymphoma/${SeqRun}/FASTQs
cd ../../..
scp SampleSheet.csv ${usr}@mercury.pmacs.upenn.edu:/project/cpdlab/Lymphoma/$SeqRun

#seq_id="SEQ-"`awk -F "," 'NR==4 {print $2}' SampleSheet.csv`
seq_clean=$(echo $seq_id | sed "s///")

echo "bsub -q cpd-dev sh docker_lymphoma_parse.sh $1 ${seq_clean}" | ssh ${usr}@consign.pmacs.upenn.edu "cat >> /project/cpdlab/Lymphoma/${seq_clean}_run_this.sh"
