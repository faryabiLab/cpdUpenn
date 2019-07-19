#!/bin/bash
#
# Analysis archiver
#
# This script will check if analysis archives are present
# for all sequencing events. If they are not, it will
# attempt to create an archive. If an archive exists,
# the corresponding analysis directory will be removed
# from Isilon. Otherwise the file will be placed in the quarantine
# folder(s) for manual inspection. This script will only act on files
# older than 60 days. Archives are organized by year.

# Associated Paths (note: Archer not included for now)
from_hpc="/PathCPD/FromHPC/"
miseq_base="/PathCPD/FromHPC/Analysis_Archives/MiSeq_analysis_archive/"
hiseq_base="/PathCPD/FromHPC/Analysis_Archives/HiSeq_analysis_archive/"
quarantine="/PathCPD/FromHPC/quarantine.txt"

test_delete="/PathCPD/FromHPC/test_delete.txt"

# Get the current year and organize logs and directories accordingly.
year=$(date -d -90days | cut -f6 -d " ")"/"
mkdir -p ${miseq_base}${year}
mkdir -p ${hiseq_base}${year}
miseq_archive=${miseq_base}${year}
miseq_log=${miseq_archive}miseq_archival_log.txt

hiseq_archive=${hiseq_base}${year}
hiseq_log=${hiseq_archive}hiseq_archival_log.txt

cd /PathCPD/FromHPC/
#seq_60=$(find ${from_hpc} -maxdepth 1 -mindepth 1 -type d -ctime +120)
seq_60=$(find ${from_hpc}1* -maxdepth 0 -mindepth 0 -type d -ctime +135)
for s in ${seq_60}; do
  # This round is just for MiSeq runs! HiSeq runs are one level further down.
  seq_dir=$(basename ${s})
  if [[ $(find ${miseq_archive}${seq_dir}_analysis_archive.tar.gz -type f -size +500M 2>/dev/null) ]]; then # If the archive's already there...
    echo "${seq_dir},Already archived,$(date -r ${miseq_archive}${seq_dir}_analysis_archive.tar.gz)" >> ${miseq_log}
    sudo chown -R svc_cpdlab:svc_cpdlab ${seq_dir}
    #sudo rm -rf ${s}
    ssh chitturi@mercury.pmacs.upenn.edu "mv /project/cpdlab/*/${seq_dir} /project/cpdlab/quarantine/"
    echo "${seq_dir}" >> ${test_delete}
    echo "${seq_dir},Removed,$(date)" >> ${miseq_log}
  else
    echo "Miseq analysis archive not found, going to go ahead and make one for ${seq_dir}" # archive everything in analysis folder
    sudo chown -R svc_cpdlab:svc_cpdlab ${seq_dir}
    tar -czvf ${miseq_archive}${seq_dir}_analysis_archive.tar.gz $(find ${s}/*)
    if [[ $(find ${miseq_archive}${seq_dir}_analysis_archive.tar.gz -type f -size +500M 2>/dev/null) ]]; then # If it's now there + over 50 MB...
      echo "${seq_dir},Archived,$(date -r ${miseq_archive}${seq_dir}_analysis_archive.tar.gz)" >> ${miseq_log}
      ssh chitturi@mercury.pmacs.upenn.edu "mv /project/cpdlab/*/${seq_dir} /project/cpdlab/quarantine/"
      #sudo rm -rf ${s}
      echo "entering conditional 2, the ${seq_dir} newly been archived!"
      echo "${seq_dir}" >> ${test_delete} # delete folder from Isilon
      echo "${seq_dir},Removed,$(date)" >> ${miseq_log} # record in log
    else # If it's there but under 500 MB in size...
      if [[ $(find ${miseq_archive}${seq_dir}_analysis_archive.tar.gz -type f -size -500M 2>/dev/null) ]]; then
        echo "${miseq_archive}${seq_dir}_analysis_archive.tar.gz: archive deleted for being too small, check this manually" >> ${test_delete} # delete archive + complain
        echo "${seq_dir} archive is too small, will delete attempted archive but not the original." >> ${miseq_log} # leave seq folder in place on Isilon
      fi
    fi
    echo "quarantining ${seq_dir}"
    echo "${seq_dir}" >> ${quarantine}
    echo "${seq_dir},Quarantined,$(date)" >> ${miseq_log}
  fi
done

cd /PathCPD/FromHPC/HiSeqRun
hiseq_60=$(find /PathCPD/FromHPC/HiSeqRun/Project_SEQ_* -maxdepth 0 -mindepth 0 -type d -mtime +80)
#hiseq_60=$(find /PathCPD/FromHPC/HiSeqRun/ -maxdepth 1 -mindepth 1 -type d -mtime +40)
for s in ${hiseq_60}; do
  hiseq_dir=$(basename ${s})
  if [[ $(find ${hiseq_archive}${hiseq_dir}_analysis_archive.tar.gz -type f -size +500M 2>/dev/null) ]]; then
    echo "${hiseq_dir},Already archived,$(date -r ${hiseq_archive}${hiseq_dir}_analysis_archive.tar.gz)" >> ${hiseq_log}
    sudo chown -R svc_cpdlab:cpdlab ${hiseq_dir}
    #sudo rm -rf ${s}
    ssh chitturi@mercury.pmacs.upenn.edu "mv /project/cpdlab/HiSeqRun/${hiseq_dir} /project/cpdlab/quarantine/"
    echo "${hiseq_dir}" >> ${test_delete}
    echo "${hiseq_dir},Removed,$(date)" >> ${hiseq_log}
  else
    echo "Hiseq analysis archive not found, going to go ahead and make one for ${hiseq_dir}"
    sudo chown -R svc_cpdlab:cpdlab ${hiseq_dir}
    tar -czvf ${hiseq_archive}${hiseq_dir}_analysis_archive.tar.gz $(find ${s}/*) # archive everything in given analysis folder
    ssh chitturi@mercury.pmacs.upenn.edu "mv /project/cpdlab/HiSeqRun/${hiseq_dir} /project/cpdlab/quarantine/"
    if [[ $(find ${hiseq_archive}${hiseq_dir}_analysis_archive.tar.gz -type f -size +500M 2>/dev/null) ]]; then
      echo "${hiseq_dir},Archived,$(date -r ${hiseq_archive}${hiseq_dir}_analysis_archive.tar.gz)" >> ${hiseq_log}
      #sudo rm -rf ${s}
      echo "entering conditional 2, the ${hiseq_dir} newly been archived!"
      echo "${hiseq_dir}" >> ${test_delete} # delete folder from Isilon
      echo "${hiseq_dir},Removed,$(date)" >> ${hiseq_log} # record in log
    else
      if [[ $(find ${hiseq_archive}${hiseq_dir}_analysis_archive.tar.gz -type f -size -500M 2>/dev/null) ]]; then
        echo "${hiseq_archive}${hiseq_dir}_analysis_archive.tar.gz: archive deleted for being too small, check this manually" >> ${test_delete} # delete archive + complain
        echo "${hiseq_dir} archive is too small, will delete attempted archive but not the original." >> ${hiseq_log} # leave seq folder in place on Isilon
      fi
    fi
    echo "quarantining ${hiseq_dir}"
    echo "${hiseq_dir}" >> ${quarantine}
    echo "${hiseq_dir},Quarantined,$(date)" >> ${hiseq_log}
  fi
done
