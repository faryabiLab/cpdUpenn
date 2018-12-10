#!/bin/sh

dir=$1
run_dir=/PathCPD/illumina/SeqOut/${dir}
out_dir=/data1/HiSeqRun/${dir}
hpc_dir='/project/cpdlab/HiSeqRun'
hpc_s2_script='/project/cpdlab/Scripts/SolidV2/Halo_v1.3.sh'
variants_header='/data1/HiSeqRun/headers_filemaker_upload/variants_header.txt'
amps_header='/data1/HiSeqRun/headers_filemaker_upload/amps_header.txt'
stats_header='/data1/HiSeqRun/headers_filemaker_upload/stats_header.txt'
depths_header='/data1/HiSeqRun/headers_filemaker_upload/depths.txt'
failed_header='/data1/HiSeqRun/headers_filemaker_upload/failed_samples.txt'


#fastq archive location
archive_dir='/PathCPD/illumina/SeqOut/HiSeq_fastq_archive'



mkdir ${out_dir}
cp ${run_dir}/SampleSheet.csv ${run_dir}/Data/Intensities/BaseCalls/
mv ${run_dir}/Data/Intensities/BaseCalls/SampleSheet.csv ${run_dir}/Data/Intensities/BaseCalls/SampleSheet.tmp.csv
sed -i 1,18d ${run_dir}/Data/Intensities/BaseCalls/SampleSheet.tmp.csv
awk -F ',' 'BEGIN{OFS=",";}{NR>1(gsub("_","-",$3))}1' ${run_dir}/Data/Intensities/BaseCalls/SampleSheet.tmp.csv > ${run_dir}/Data/Intensities/BaseCalls/SampleSheet.tmp2.csv
cut -d ',' -f1,2,3,4,6,9,10,11,12,13 ${run_dir}/Data/Intensities/BaseCalls/SampleSheet.tmp2.csv > ${run_dir}/Data/Intensities/BaseCalls/SampleSheet.csv 
sed -i 's/Sample-ID/Sample_ID/g' ${run_dir}/Data/Intensities/BaseCalls/SampleSheet.csv
rm ${run_dir}/Data/Intensities/BaseCalls/SampleSheet.tmp*

configureBclToFastq.pl --input-dir ${run_dir}/Data/Intensities/BaseCalls --output-dir ${out_dir}/Unaligned --sample-sheet ${run_dir}/Data/Intensities/BaseCalls/SampleSheet.csv --no-eamss --use-bases-mask Y150n,I8,Y10,Y150n --mismatches 1

cd ${out_dir}/Unaligned
nohup make -j 4 > demux.out


project=Project_$(awk -F ',' 'NR==2 {print $10}' ${run_dir}/Data/Intensities/BaseCalls/SampleSheet.csv | grep -v 'Sample_Project' | sort -u | tr '\r' ' ')
project="$(echo -e "${project}" | tr -d '[[:space:]]')"
project_dir=${out_dir}/Unaligned/${project}
cd ${project_dir}
#where read 1 is read 1, read2 is the UMI, and read3 is read2
for j in `awk -F "," '{print $3}' ${run_dir}/Data/Intensities/BaseCalls/SampleSheet.csv | grep -v 'Sample_ID' | sort |uniq`
do
zcat Sample_$j/$j*R1*gz > Sample_$j/$j.R1.fastq
echo "R1 done for Sample_$j" >> demux.out
zcat Sample_$j/$j*R2*gz > Sample_$j/$j.R2.fastq
echo "R2 done for Sample_$j" >> demux.out
zcat Sample_$j/$j*R3*gz > Sample_$j/$j.R3.fastq
echo "R3 done for Sample_$j" >> demux.out
gzip Sample_$j/*fastq
echo "bsub -q cpd -n 10 -M 45000 -R "span[hosts=1]" -e ${hpc_dir}/${project}/Sample_$j/$j.e -o ${hpc_dir}/${project}/Sample_$j/$j.o sh ${hpc_s2_script} ${hpc_dir}/${project}/Sample_$j/" >> ${project_dir}/bsubs_${project}.sh
#echo "bsub -n 20 -M 80000 -R "span[hosts=1]" -e ${hpc_dir}/${project}/Sample_$j/$j.e -o ${hpc_dir}/${project}/Sample_$j/$j.o sh ${hpc_s2_flfix_script} ${hpc_dir}/${project}/Sample_$j/" >> ${project_dir}/bsubs_${project}.flfix.sh

done
chmod 775 ${project_dir}/bsubs_${project}.sh
rsync -a --include '*.R*.fastq.gz' --include='SampleSheet.csv' --include='bsubs*' --include '*/' --exclude='*' ${project_dir}  bigdelia@mercury.pmacs.upenn.edu:/${hpc_dir}
ssh bigdelia@consign.pmacs.upenn.edu "bash -ic" ${hpc_dir}/${project}/bsubs_${project}.sh

#generate files for concat post hpc run
hpc_output=/PathCPD/FromHPC/HiSeqRun/${project}

mkdir ${hpc_output}
cp ${variants_header} ${hpc_output}/${project}_Run_masterVarFinal.txt
cp ${amps_header} ${hpc_output}/${project}_Amplifications.txt
cp ${stats_header} ${hpc_output}/${project}_RunStatsFinal.txt
cp ${depths_header} ${hpc_output}/${project}_MeanCovByPostion.txt
cp ${failed_header} ${hpc_output}/${project}_FailedSamples.txt


rsync -a --include '*.R*.fastq.gz' --include='SampleSheet.csv' --include='bsubs*' --include '*/' --exclude='*' ${project_dir} ${archive_dir}
tar --remove-files -czvf ${project_dir}_fastq_archive.tar.gz -C ${archive_dir} ${project_dir}
