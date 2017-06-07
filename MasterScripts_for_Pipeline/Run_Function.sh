#!/bin/sh
##Run this script as : ./Run_Function.sh <RunFolderName> <PipelineVersion> <Minimym_Memory> <Maximum_Memory> <Num_Of_Threads> <Status> <config_file>
rfld=$1
pipeV=$2
m_mem=$3
M_mem=$4
N_td=$5
rfld_status=$6
config_file=$7
if [ -d /project/cpdlab/RunFolders/$rfld ];
then
	if [ -s /project/cpdlab/RunFolders/$rfld/AnalysisLog.txt ];
	then
		#sed -i "s/$rfld,$pipeV,$m_mem,$M_mem,$N_td,$rfld_status/$rfld,$pipeV,$m_mem,$M_mem,$N_td,sleeping_20min/g" $config_file
		#sleep 20m
		sed -i "s/$rfld,$pipeV,$m_mem,$M_mem,$N_td,$rfld_status/$rfld,$pipeV,$m_mem,$M_mem,$N_td,performing_copying_and_submitting_jobs/g" $config_file
		echo "Running Create_Run_Merge.sh ..."
		bsub -e /project/cpdlab/RunFolders/$rfld/$rfld"_"$pipeV.e -o /project/cpdlab/RunFolders/$rfld/$rfld"_"$pipeV.o sh Create_Run_Merge.sh /project/cpdlab/RunFolders/$rfld /project/cpdlab/Scripts/Pipeline$pipeV $m_mem $M_mem $N_td
		sed -i "s/$rfld,$pipeV,$m_mem,$M_mem,$N_td,performing_copying_and_submitting_jobs/$rfld,$pipeV,$m_mem,$M_mem,$N_td,job_submitted/g" $config_file
		echo "Copying $config_file to /PathCPD/FromHPC ..."
		rsync --stats --progress --bwlimit=10000 -pgt $config_file cpdlab@170.212.141.107:/PathCPD/FromHPC/
	else
		echo "Copying AnalysisLog.txt from /PathCPD/illumina/SeqOut/$rfld to /project/cpdlab/RunFolders/$rfld/ ..."
		rsync --stats --progress --bwlimit=10000 -pgt cpdlab@170.212.141.107:/PathCPD/illumina/SeqOut/$rfld/AnalysisLog.txt /project/cpdlab/RunFolders/$rfld/
	fi
else
	mkdir /project/cpdlab/RunFolders/$rfld
	echo "Copying RunFolder to /PathCPD/FromHPC ..."
	rsync --stats --progress --bwlimit=10000 -rpgt /project/cpdlab/RunFolders/$rfld cpdlab@170.212.141.107:/PathCPD/FromHPC/
fi
