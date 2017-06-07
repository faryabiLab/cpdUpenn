#!/bin/sh
##Run this script as : ./HPC_Master.sh <Run_Config.txt>
config_file=$1

for i in `awk '{print}' $config_file|grep -v '^#'|grep 'pending'`
do
        RunFolder=`echo $i|cut -d "," -f 1`
        pipe_version=`echo $i|cut -d "," -f 2`
        Minimum_Memory=`echo $i|cut -d "," -f 3`
        Maximum_Memory=`echo $i|cut -d "," -f 4`
        Num_Of_Threads=`echo $i|cut -d "," -f 5`
	RunFolderStatus=`echo $i|cut -d "," -f 6`
        bsub -o /project/cpdlab/RunFolders/$RunFolder/Run_Function_$RunFolder.o -e /project/cpdlab/RunFolders/$RunFolder/Run_Function_$RunFolder.e sh Run_Function.sh $RunFolder $pipe_version $Minimum_Memory $Maximum_Memory $Num_Of_Threads $RunFolderStatus $config_file
	sleep 1m
done
