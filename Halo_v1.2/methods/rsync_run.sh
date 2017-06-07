rsync -a --include='*_final.bam*' --include='SEQ*' --exclude='*logs/' --exclude='*/tmp' --exclude='*fastqc/' --include='*/' --exclude='*' $1  cpdlab@170.212.141.107:/PathCPD/FromHPC/HiSeqRun/
