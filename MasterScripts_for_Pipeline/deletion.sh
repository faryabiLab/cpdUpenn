#!/bin/sh
##Run this script as : ./HPC_Master.sh <PATH>
work_path=$1

for i in `ls $work_path|grep '1*_M0*-*'`
do
if [ -d $work_path/$i ];
then
rm -r -f $work_path/$i
fi
done
