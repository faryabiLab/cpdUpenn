#!/bin/sh

run_dir=/PathCPD/illumina/SeqOut/$1
out_dir=/data1/HiSeqRun/ARCHER/$1
samplesheet=$2

mkdir ${out_dir}


cp ${run_dir}/${samplesheet} ${run_dir}/Data/Intensities/BaseCalls/
mv ${run_dir}/Data/Intensities/BaseCalls/${samplesheet} ${run_dir}/Data/Intensities/BaseCalls/${samplesheet}.tmp
sed -i 1,18d ${run_dir}/Data/Intensities/BaseCalls/${samplesheet}.tmp
awk -F ',' 'BEGIN{OFS=",";}{NR>1(gsub("_","-",$3))}1' ${run_dir}/Data/Intensities/BaseCalls/${samplesheet}.tmp > ${run_dir}/Data/Intensities/BaseCalls/${samplesheet}.tmp2
awk -F "," '{OFS=","; print $1,$2,$3,$4,$6"-"$8,$9,$10,$11,$12,$13}' ${run_dir}/Data/Intensities/BaseCalls/${samplesheet}.tmp2 > ${run_dir}/Data/Intensities/BaseCalls/${samplesheet}
sed -i 's/Sample-ID/Sample_ID/g' ${run_dir}/Data/Intensities/BaseCalls/${samplesheet}
rm ${run_dir}/Data/Intensities/BaseCalls/${samplesheet}.tmp*

configureBclToFastq.pl --input-dir ${run_dir}/Data/Intensities/BaseCalls --output-dir ${out_dir}/Unaligned --sample-sheet ${run_dir}/Data/Intensities/BaseCalls/${samplesheet} --no-eamss --use-bases-mask Y150n,I8,I8,Y150n --mismatches 1

cd ${out_dir}/Unaligned
nohup make > demux.out


project=Project_$(awk -F ',' 'NR==2 {print $10}' ${run_dir}/Data/Intensities/BaseCalls/${samplesheet} | grep -v 'Sample_Project' | sort -u | tr '\r' ' ')
project="$(echo -e "${project}" | tr -d '[[:space:]]')"
project_dir=${out_dir}/Unaligned/${project}
cd ${project_dir}
for j in `awk -F "," '{print $3}' ${run_dir}/Data/Intensities/BaseCalls/${samplesheet} | grep -v 'Sample_ID' | sort |uniq`
do
zcat Sample_$j/$j*R1*gz > Sample_$j/$j'_R1.fastq'
echo "R1 done for Sample_$j" >> demux.out
zcat Sample_$j/$j*R2*gz > Sample_$j/$j'_R2.fastq'
echo "R2 done for Sample_$j" >> demux.out
cat Sample_$j/$j*R*.fastq | echo $j $((`wc -l`/4)) >> ${project}_read_counts.txt
gzip Sample_$j/*fastq
done

mkdir /PathCPD/FromHPC/HiSeqRun/ARCHER_fastqs/${project}

scp ${out_dir}/Unaligned/${project}/*/*_R1.fastq.gz /PathCPD/FromHPC/HiSeqRun/ARCHER_fastqs/${project}/
scp ${out_dir}/Unaligned/${project}/*/*_R2.fastq.gz /PathCPD/FromHPC/HiSeqRun/ARCHER_fastqs/${project}/
scp ${out_dir}/Unaligned/${project}/${project}_read_counts.txt /PathCPD/FromHPC/HiSeqRun/ARCHER_fastqs/${project}/

rsync -a --include '*.R*.fastq.gz' --include='SampleSheet.csv' --include '*/' --exclude='*' ${project_dir} ${archive_dir}
tar --remove-files -czvf ${project_dir}_fastq_archive.tar.gz -C ${archive_dir} ${project_dir}

