##V1.1
## Run this script as ./egfrPick.sh <SampleName> <PrimerManifestFile> <PrimerFile> <PrimerRCFile>
## Changes in V1.1: 1) Removed Primer.PrimParts and Primer_RC.PrimParts files from the command-line as they are not used in this script at all.
##		    2) Removed the steps to generate HPRT1,GUCB and ACTB amplicon files as these amplicons are not covered in the current assay. Have also removed the entries relevant to these three amplicons from the Manifest.txt file.

#!/bin/sh


NameOfSample=$1
PrimManifest=$2
ActualPrimerFile=$3
ActualPrimerRCFile=$4

echo "Concating Fastqs"

zcat *_R1_* > $NameOfSample.R1.fastq
zcat *_R2_* > $NameOfSample.R2.fastq

echo "Running Btrim on R2 fastq"

/home/bdaber01/Btrim64 -p $ActualPrimerRCFile -t $NameOfSample.R2.fastq -o $NameOfSample.R2.btrim.fastq -s $NameOfSample.R2.btrim.summary.txt -e 100 -w 0 > $NameOfSample.R2.btrim.stats.txt

echo "Running Btrim on R1 fastq"
/home/bdaber01/Btrim64 -p $ActualPrimerFile -t $NameOfSample.R1.fastq -o $NameOfSample.R1.btrim.fastq -s $NameOfSample.R1.btrim.summary.txt -e 100 -w 0 > $NameOfSample.R1.btrim.stats.txt




grep -v 'no_5_found' $NameOfSample.R2.btrim.summary.txt| grep -v 'short'|awk '{print $1,$4}' | sort -k 1,1 > $NameOfSample.R2.btrim.summary.awked.sorted
grep -v 'no_5_found' $NameOfSample.R1.btrim.summary.txt| grep -v 'short'|awk '{print $1,$4}' | sort -k 1,1 > $NameOfSample.R1.btrim.summary.awked.sorted

join -1 1 -2 1 $NameOfSample.R1.btrim.summary.awked.sorted $NameOfSample.R2.btrim.summary.awked.sorted > $NameOfSample.btrim.summary.joined


grep ' 0 0' $NameOfSample.btrim.summary.joined > EGFR_exon1_intron1
grep ' 0 1' $NameOfSample.btrim.summary.joined > EGFR_exon1_exon2
grep ' 0 2' $NameOfSample.btrim.summary.joined > EGFR_exon1_exon8
grep ' 3 3' $NameOfSample.btrim.summary.joined > HPRT4
grep ' 4 4' $NameOfSample.btrim.summary.joined > SDHA
grep ' 5 5' $NameOfSample.btrim.summary.joined > RPL13A
grep ' 6 6' $NameOfSample.btrim.summary.joined > EGFR_exon9_exon12
grep ' 6 7' $NameOfSample.btrim.summary.joined > EGFR_exon9_exon10
grep ' 6 8' $NameOfSample.btrim.summary.joined > EGFR_exon9_exon9

echo "egfrPicker ended..."

echo "Making Count Summary File.."
nohup ./Count.sh > $NameOfSample.AmpCount.out
awk '{print $2,$1}' $NameOfSample.AmpCount.out > $NameOfSample.AmpCount.awked.out
echo "Making Final Count-Summary File"
python MakingOfSummary.py $NameOfSample

#echo "Now running Cut Adapt..."

#/home/bdaber01/cutadapt-1.3/bin/cutadapt -a GCGAATTTCGACGATCGTTGCATTAACTCGCGA -o $NameOfSample.R1.adapTrim.fastq $NameOfSample.R1.fastq > $NameOfSample.R1.adapTrim.out
#/home/bdaber01/cutadapt-1.3/bin/cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o $NameOfSample.R2.adapTrim.fastq $NameOfSample.R2.fastq > $NameOfSample.R2.adapTrim.out

#echo "Running BWA-MEM on adap-trimmed fastqs"
#
#/home/bdaber01/bwa-0.7.3a/bwa mem -t 6 -R '@RG\tID:1\tSM:'$NameOfSample'\tPL:Illumina\tLB:TSCA\tPU:Index_'$NameOfSample'\tCN:CPD' /Drive_D/Databases/genomes/hg19_genome.fa $NameOfSample.R1.adapTrim.fastq $NameOfSample.R2.adapTrim.fastq > $NameOfSample.sam
#
#echo "converting sam to bam"
#samtools view -bS $NameOfSample.sam > $NameOfSample.bam
#echo "Sorting bam"
#samtools sort $NameOfSample.bam $NameOfSample.sorted
#echo "Indexing Bam.."
#samtools index $NameOfSample.sorted.bam
#echo "Calling stats on sorted bam"
#samtools flagstat $NameOfSample.sorted.bam > $NameOfSample.sorted.flagstat
#
#
#echo "Running Depth Of Coverage using the SnpSeq Target Bed in the current folder"
#java -Xmx6g -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -T DepthOfCoverage -o $NameOfSample.Depth -I $NameOfSample.sorted.bam --minBaseQuality 22 -baseCounts -ct 0 -ct 1 -ct 250 -ct 1000 -L TSCA_CANCER_PANEL_SnpSeq_target.BED -R /Drive_D/Databases/genomes/hg19_genome.fa
#
#echo "$SampleName processed successfully..."