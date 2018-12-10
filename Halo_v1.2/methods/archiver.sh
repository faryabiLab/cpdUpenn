seq_runs=($(find /PathCPD/illumina/SeqOut/ -maxdepth 1 -type d -mtime +30))


for dir in ${seq_runs[*]}
do
    tar cvfz ${dir}_fastq_archive.tar.gz $(find ${dir} -type f -name "bsub*"-name "*.R1.fastq.gz" -or -name "*.R2.fastq.gz" -or -name "*.R3.fastq.gz" -or -name "SampleSheet.csv")
    tar cvfz ${dir}_analysis_archive.tar.gz $(find ${dir} -type f -name "bsub*" -or -name "*final.bam*" -or -name "*TCGAMAF" -or -name "*Run*")
    rsync ${dir}_fastq_archive.tar.gz cpdlab@170.212.141.107:/PathCPD/illumina/SeqOut/HiSeq_fastq_archive/
    rsync ${dir}_fastq_archive.tar.gz cpdlab@170.212.141.107:/PathCPD/illumina/SeqOut/HiSeq_fastq_archive/
    rsync ${dir}_analysis_archive.tar.gz cpdlab@170.212.141.107:/PathCPD/FromHPC/HiSeq_analysis_archive/
    rsync ${dir}_analysis_archive.tar.gz cpdlab@170.212.141.107:/PathCPD/FromHPC/HiSeq_analysis_archive/
    rm -rf ${dir}_fastq_archive.tar.gz
    rm -rf ${dir}_analysis_archive.tar.gz 
    #rm -rf ${dir}
done
