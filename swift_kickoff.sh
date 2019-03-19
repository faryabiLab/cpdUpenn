#!/bin/bash

runID=$1

cd /PathCPD/illumina/SeqOut/${runID}/Data/Intensities/BaseCalls
rm Undetermined_S0_L001_R*gz
cd /PathCPD/illumina/SeqOut/${runID}
seqraw=`awk -F "," 'NR==4 {print $2}' SampleSheet.csv`
seqrun=$(echo $seqraw | sed "s///")

cd /PathCPD/illumina/SeqOut/${runID}/Data/Intensities/BaseCalls
ssh chitturi@mercury.pmacs.upenn.edu mkdir -p /project/cpdlab/SWIFT/${seqrun}
ssh chitturi@mercury.pmacs.upenn.edu "cat /project/cpdlab/SWIFT/stats_header.txt > /project/cpdlab/SWIFT/${seqrun}/RunStatsFinal.txt"
ssh chitturi@mercury.pmacs.upenn.edu "cat /project/cpdlab/SWIFT/variants_header.txt > /project/cpdlab/SWIFT/${seqrun}/Run_masterVarFinal.txt"
echo "cd /project/cpdlab/SWIFT/${seqrun}" | ssh chitturi@consign.pmacs.upenn.edu "cat >> /project/cpdlab/SWIFT/${seqrun}/run_this.sh"

jobnum=0
for f in *_R1_001.fastq.gz
do
  fq1=$f
  fq2=${fq1%%_R1_*}_R2_001.fastq.gz
  prefix=`cut -d'_' -f1 <<< ${fq1%%_R1_001.fastq.gz}`"-SEQ-"${seqrun}

  mkdir ${prefix}
  mv $fq1 $fq2 ${prefix}
  rsync -avPR ${prefix} chitturi@mercury.pmacs.upenn.edu:/project/cpdlab/SWIFT/${seqrun}
  job_name=$(echo ${prefix} | cut -c4-10)
  echo "bsub -q cpd -n 16 -M 25000 -e /project/cpdlab/SWIFT/${seqrun}/CPD${job_name}.e -o /project/cpdlab/SWIFT/${seqrun}/CPD${job_name}.o -J "CPD${job_name}" sh /home/chitturi/swift_parallel_test.sh ${prefix} ${seqrun}" | ssh chitturi@consign.pmacs.upenn.edu "cat >> /project/cpdlab/SWIFT/${seqrun}/run_this.sh"
  jobnum=$((jobnum + 1))
done
echo "bsub -q cpd -w \"done(CPD*)\" rsync -avPR * cpdlab@170.212.141.107:/PathCPD/FromHPC/${seqrun}/" | ssh chitturi@consign.pmacs.upenn.edu "cat >> /project/cpdlab/SWIFT/${seqrun}/run_this.sh"
ssh chitturi@consign.pmacs.upenn.edu "sh /project/cpdlab/SWIFT/${seqrun}/run_this.sh"
