#!/bin/bash
#
# look for a way to not have to run this as root
# Archiver.sh
#
# Created by Akshay Chitturi 7/9/2019
#
# This script will check if analysis archives are present
# for all sequencing events. If they are not, it will
# attempt to create an archive. If an archive exists,
# the corresponding analysis directory will be removed
# from the VM. Otherwise the file will be placed in the quarantine
# folder(s) for manual inspection. This script will only act on files
# older than 60 days. Archives are organized by year.

# Associated Paths
analysis="/var/www/analysis/"
archer_base="/PathCPD/FromHPC/Analysis_Archives/Archer_analysis_archive/"
quarantine="/var/www/analysis/quarantine.txt"

test_delete="/var/www/analysis/test_delete.txt"

# Get the current year and organize logs and directories accordingly.
year=$(date -d -90days | cut -f6 -d " ")"/"
mkdir -p ${archer_base}${year}
archer_archive=${archer_base}${year}
archer_log=${archer_base}${year}archer_archival_log.txt

cd /var/www/analysis/
arch_60=$(find 4* -maxdepth 0 -mindepth 0 -type d -ctime +60) # expand to numeric check
for s in ${arch_60}; do
  seq_dir=$(basename ${s})
  if [[ $(find ${archer_archive}${seq_dir}_analysis_archive.tar.gz -type f -size +500M 2>/dev/null) ]]; then
    echo "${seq_dir},Already archived,$(date -r ${archer_archive}${seq_dir}_analysis_archive.tar.gz)" >> ${archer_log}
    #sudo rm -rf ${s}
    echo "${seq_dir}" >> ${test_delete}
    echo "${seq_dir},Removed,$(date)" >> ${archer_log}
  else
    echo "Archer analysis archive not found, going to go ahead and make one for ${seq_dir}"
    find ${s}/* > filelist_temp.txt
    tar -cvz -T filelist_temp.txt -f ${archer_archive}${seq_dir}_analysis_archive.tar.gz
    if [[ $(find ${archer_archive}${seq_dir}_analysis_archive.tar.gz -type f -size +500M 2>/dev/null) ]]; then
      echo "${seq_dir},Archived,$(date -r ${archer_archive}${seq_dir}_analysis_archive.tar.gz)" >> ${archer_log}
      #sudo rm -rf ${s}
      echo "entering conditional 2, the ${seq_dir} newly been archived!"
      echo "${seq_dir}" >> ${test_delete} # delete folder from Isilon
      echo "${seq_dir},Removed,$(date)" >> ${archer_log} # record in log
    else
      if [[ $(find ${archer_archive}${seq_dir}_analysis_archive.tar.gz -type f -size -500M 2>/dev/null) ]]; then
        echo "${archer_archive}${seq_dir}_analysis_archive.tar.gz: archive deleted for being too small, check this manually" >> ${test_delete}
        echo "${seq_dir} archive is too small, will delete attempted archive but not the original." >> ${archer_log} # leave seq folder in place on Isilon
      fi
    fi
    echo "checking off ${seq_dir}"
    echo "${seq_dir}" >> ${quarantine}
    echo "${seq_dir},Quarantined,$(date)" >> ${archer_log}
  fi
done
