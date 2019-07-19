#!/bin/bash
#
#
# Archiver.sh
#
# Created by Ashkan Bigdeli 12/12/2018
#
# This script will check if fastq archives are present
# for all sequencing events. If they are not, it will
# attempt to create an archive. If an archive exists,
# the correspond sequencing directory will be removed
# from all file systems. Otherwise the file will be placed in the quarantine
# folder(s) for manual inspection. This script will only act on files
# older than 60 days. Archives are organized by year & Hi||MiSeq.

# Associated Paths
seq_out="/PathCPD/illumina/SeqOut/"
demux_out="/data1/HiSeqRun/"
miseq_base="/PathCPD/illumina/SeqOut/MiSeq_fastq_archive/"
hiseq_base="/PathCPD/illumina/SeqOut/HiSeq_fastq_archive/"
quarantine_seq="/PathCPD/illumina/SeqOut/Quarantine/quarantine.txt"
quarantine_demux="/data1/HiSeqRun/Quarantine/"

test_delete_demux="/data1/HiSeqRun/delete_test/"
test_delete_isilon="/PathCPD/illumina/SeqOut/delete_test.txt"

# Get the current year and organize logs and directories accordingly.
year=$(date -d -90days | cut -f6 -d " ")"/"
mkdir -p ${miseq_base}${year}
mkdir -p ${hiseq_base}${year}
miseq_archive=${miseq_base}${year}
miseq_log=${miseq_archive}miseq_archival_log.txt

hiseq_archive=${hiseq_base}${year}
hiseq_log=${hiseq_archive}hiseq_archival_log.txt

cd /PathCPD/illumina/SeqOut/
seq_60=$(find ${seq_out} -maxdepth 1 -mindepth 1 -type d -ctime +103) # Five seq runs should be processed. Only three should be archived.
for s in ${seq_60}; do
    # MiSeq + HiSeq Archival
    seq_dir=$(basename ${s})
    seq_type=$(echo ${seq_dir} | cut -f2 -d "_")
    if [ ${seq_type} = "M0333" ] || [ ${seq_type} = "M00351" ] || [ ${seq_type} = "M01651" ] || [ ${seq_type} = "M03738" ] || [ ${seq_type} = "M06394" ] || [ ${seq_type} = "M00552" ] || [ ${seq_type} = "M06371" ]; then
         if [[ $(find ${miseq_archive}${seq_dir}_fastq_archive.tar.gz -type f -size +500M 2>/dev/null) ]]; then # If the archive's already there...
             echo "${seq_dir},Already archived,$(date -r ${miseq_archive}${seq_dir}_fastq_archive.tar.gz)" >> ${miseq_log}
             #sudo rm -rf ${s}
             echo "${seq_dir}" >> ${test_delete_isilon}
             echo "${seq_dir},Removed,$(date)" >> ${miseq_log} # record in log, delete seq folder from Isilon
         else # If it's NOT already there...
             echo "Miseq archive not found, making one for ${seq_dir}"
             cd ${s}/Data/Intensities/BaseCalls/
             tar -czvf ${miseq_archive}${seq_dir}_fastq_archive.tar.gz $(ls -- *R*.fastq.gz SampleSheet.csv) # archive FASTQs + SS
             cd /PathCPD/illumina/SeqOut/
             if [[ $(find ${miseq_archive}${seq_dir}_fastq_archive.tar.gz -type f -size +500M 2>/dev/null) ]]; then # If it's now there + over 500 MB...
                 echo "${seq_dir},Archived,$(date -r ${miseq_archive}${seq_dir}_fastq_archive.tar.gz)" >> ${miseq_log}
                 #sudo rm -rf ${s}
                 echo "${seq_dir} has just been archived."
                 echo "${seq_dir}" >> ${test_delete_isilon} # delete seq folder from Isilon
                 echo "${seq_dir},Removed,$(date)" >> ${miseq_log} # record in log
             else # If it's there but under 500 MB in size...
                 if [[ $(find ${miseq_archive}${seq_dir}_fastq_archive.tar.gz -type f -size -500M 2>/dev/null) ]]; then
                     echo "${miseq_archive}${seq_dir}_fastq_archive.tar.gz: archive deleted for being too small, check this manually" >> ${test_delete_isilon} # delete archive + complain
                     echo "${seq_dir} archive is too small, will delete attempted archive but not the original." # leave seq folder in place on Isilon
                 fi
             fi
             echo "checking off ${seq_dir}"
             echo "${seq_dir}" >> ${quarantine_seq}
             echo "${seq_dir},Quarantined,$(date)" >> ${miseq_log}
         fi
    elif [ ${seq_type} = "D00831" ] || [ ${seq_type} = "D00832" ]; then
         if [[ $(find ${hiseq_archive}${seq_dir}_fastq_archive.tar.gz -type f -size +2000M 2>/dev/null) ]]; then
             echo "${seq_dir},Already archived,$(date -r ${hiseq_archive}${seq_dir}_fastq_archive.tar.gz)" >> ${hiseq_log}
             #sudo rm -rf ${s}
             echo "${seq_dir}" >> ${test_delete_isilon}
             echo "${seq_dir},Removed,$(date)" >> ${hiseq_log} # record in log, delete seq folder from Isilon
         else
             echo "HiSeq archive not found, will make one for ${seq_dir}"
             if sudo grep -q "Archer" ${s}/SampleSheet.csv; then
               cd ${s}/Data/Intensities/BaseCalls/ # only for Archer!
               tar -czvf ${hiseq_archive}${seq_dir}_fastq_archive.tar.gz $(ls -- *R*.fastq.gz SampleSheet.csv)
               cd /PathCPD/illumina/SeqOut/
             elif sudo grep -q "Solid" ${s}/SampleSheet.csv; then
               cd /data1/HiSeqRun/${seq_dir}/Unaligned/Project_SEQ_*
               tar -czvf ${hiseq_archive}${seq_dir}_fastq_archive.tar.gz $(ls */*R1.fastq.gz */*R2.fastq.gz */*R3.fastq.gz)
               cd /PathCPD/illumina/SeqOut/
             fi 
             if [[ $(find ${hiseq_archive}${seq_dir}_fastq_archive.tar.gz -type f -size +2000M 2>/dev/null) ]]; then # If it's now there + over 2000 MB...
                 echo "${seq_dir},Archived,$(date -r ${hiseq_archive}${seq_dir}_fastq_archive.tar.gz)" >> ${hiseq_log}
                 #sudo rm -rf ${s}
                 #sudo rm -rf /data1/HiSeqRun/${seq_dir}/Unaligned
                 echo "entering conditional 2, the ${seq_dir} newly been archived!"
                 echo "${seq_dir} gone from both /PathCPD/illumina/SeqOut/ and /data1/HiSeqRun/" >> ${test_delete_isilon}  # delete seq folder from Isilon
                 echo "${seq_dir},Removed,$(date)" >> ${hiseq_log} # record in log
             else # If it's there but under 2000 MB in size...
                 if [[ $(find ${hiseq_archive}${seq_dir}_fastq_archive.tar.gz -type f -size -2000M 2>/dev/null) ]]; then
                     echo "${hiseq_archive}${seq_dir}_fastq_archive.tar.gz" >> ${test_delete_isilon} # delete archive + complain
                     echo "${seq_dir} archive is too small, is there anything there?" >> ${hiseq_log} # leave seq folder in place on Isilon
                 fi
             fi
             echo "checking off ${seq_dir}"
             echo "${seq_dir}" >> ${quarantine_seq}
             echo "${seq_dir},Quarantined,$(date)" >> ${hiseq_log}
         fi 
    fi
done
