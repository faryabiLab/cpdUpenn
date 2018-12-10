folder1=$1 #illimuna/SeqOut folder name
folder2=$2 #hpc /cpdlab/HiSeqRun folder name

#sample_sheet=$(ssh cpdlab@170.212.141.107 find /PathCPD/illumina/SeqOut/${folder1}/SampleSheet.csv)

scp cpdlab@170.212.141.107://PathCPD/illumina/SeqOut/${folder1}/SampleSheet.csv /project/cpdlab/HiSeqRun/${folder2} 

sed '20,$!d' ./${folder2}/SampleSheet.csv | sed -r 's/^.{12}//' | sed 's/,.*//' | sort | uniq > ./${folder2}/seqout_sample_id.txt


sed -i -e 's/_/-/g' ./${folder2}/seqout_sample_id.txt

find /project/cpdlab/HiSeqRun/${folder2}/Sample_* -maxdepth 0 | sed -r 's/^.{61}//' > ./${folder2}/hpc_sample_id.txt

file1=/project/cpdlab/HiSeqRun/${folder2}/seqout_sample_id.txt
file2=/project/cpdlab/HiSeqRun/${folder2}/hpc_sample_id.txt

#wc1=$(awk -vORS=, '{ print $1 }' ${file1} | sed 's/,$/\n/')
#wc2=$(awk -vORS=, '{ print $1 }' ${file2} | sed 's/,$/\n/')

hash_1=$(md5sum ${file1} | awk '{print $1}')
hash_2=$(md5sum ${file2} | awk '{print $1}')
#summarize=$(sh /project/cpdlab/Scripts/SolidV2/solidv2_v1.2/methods/summarize.sh)

#rsync=$(sh /project/cpdlab/Scripts/SolidV2/solidv2_v1.2/methods/rsync_run.sh)

#contam_check=$(python /project/cpdlab/Scripts/SolidV2/solidv2_v1.2/methods/contamination_check.py)

#echo ${wc1}
#echo ${wc2}
if [ $hash_1 = $hash_2 ]; then
    echo The files match
    sh /project/cpdlab/Scripts/SolidV2/solidv2_v1.2/methods/summarize.sh ${folder2} /project/cpdlab/HiSeqRun/${folder2}/${folder2}
    python /project/cpdlab/Scripts/SolidV2/solidv2_v1.2/methods/contamination_check.py ${folder2}
    sh /project/cpdlab/Scripts/SolidV2/solidv2_v1.2/methods/rsync_run.sh ${folder2}
else
    echo The files are different
    comm -2 -3 <(sort file1) <(sort file2) >./${folder2}/missing_samples.txt
fi
