#pass the miseq run folder name and the HPC username
#syncs relevant vcf's and bams to user's home folder

# Input arguments: SeqOut input directory, CONSIGN username

rsync -a --include='*.vcf' --include='*AmpliconCoverage*' --include='*bam*' --exclude='*/' --exclude='*'  /PathCPD/illumina/SeqOut/$1/Data/Intensities/BaseCalls/Alignment/ $2@consign.pmacs.upenn.edu:/project/cpdlab/Lymphoma/$1
rsync -a --include='*.fastq.gz' --exclude='*/' --exclude='*' /PathCPD/illumina/SeqOut/$1/Data/Intensities/BaseCalls/ $2@consign.pmacs.upenn.edu:/project/cpdlab/Lymphoma/$1/FASTQs
scp /PathCPD/illumina/SeqOut/$1/SampleSheet.csv $2@consign.pmacs.upenn.edu:/project/cpdlab/Lymphoma/$1/

cd /PathCPD/illumina/SeqOut/$1
seq_id="SEQ-"`awk -F "," 'NR==4 {print $2}' SampleSheet.csv`
seq_clean=$(echo $seq_id | sed "s///")

echo "bsub -q cpd-dev sh new_lymphoma_parse.sh $1 ${seq_clean}" | ssh $2@consign.pmacs.upenn.edu "cat >> /project/cpdlab/Lymphoma/${seq_clean}_run_this.sh"

