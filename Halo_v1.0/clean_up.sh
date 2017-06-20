home='/project/cpdlab/HiSeqRun'
run_dir=$1

fastq_dir=${run_dir}_fastq_archive
analysis_dir=${run_dir}_analysis_archive
mkdir ${fastq_dir}
mkdir ${analysis_dir}
#sync fastqs
rsync -a --include='*.R*.fastq.gz' --include='SampleSheet.csv' --exclude='*logs/' --exclude='*fastqc/' --include='*/' --exclude='*'  ${run_dir}/ ${fastq_dir}

cd ${fastq_dir} 
find . -maxdepth 1 -mindepth 1 -type d -exec tar zcvf {}.tar.gz {} --remove-files \;

cd ${home}
rsync -a --include='*bsub*' --include='*Run*' --include='*MeanCovByPostion.txt' --include='*Amplifications.txt' --include='*_final.bam*' --include='*amps.tsv' --include='*Stats' --include='*_final.vcf' --include='*target.depth.MeanCoverageByPosition' --include='*TCGAMAG' --include='*_parsed.tsv' --exclude='*logs/' --exclude='*fastqc/' --include='*/' --exclude='*' ${run_dir}/ ${analysis_dir}

cd ${analysis_dir}
find . -maxdepth 1 -mindepth 1 -type d -exec tar zcvf {}.tar.gz {} --remove-files \;

cd ${home}
tar -zcvf ${fastq_dir}.tar.gz ${fastq_dir}
tar -zcvf ${analysis_dir}.tar.gz ${analysis_dir}
