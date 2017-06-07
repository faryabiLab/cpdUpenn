##V1.1
## Changes in V1.1 : 1) 'EGFR_exon_8' has been replaced with 'EGFR_exon1_exon8' in the respective file-names. This change is to be in concordance with the amplicon naming scheme adopted in Manifest.txt and egfrPick.sh.
##                   2) 'EGFR_exon1_exon2' files are also created by searching the sequence present in EGFRwt_ToSearch.txt.
##                   3) The name of the sequence file containing EGFRvIII sequence has been changed from 'seqToSearch.txt' to 'EGFRvIII_ToSearch.txt'.
##                   4) R2 pivot final files are also generated.
##                   5) R1 and R2 fasta files are generated so that they could be utilized to do any further in-depth analysis if required.
##                   6) Extra comments are added to the script.
##                   7) 'OFS="\t"' was removed from the awk command operated on seq.R1(orR2).txt files as mentioning of this parameter was really not required for the generation of desired output. Even when it was present in the previous version it didnt have any impact on the output generated.
##                   8) This script also generates pivot.noMatch.txt files that contain the sequences that has no exact match with EGFRvIII sequence (for EGFRvIII no-match files) and EGFRwt sequence (for EGFRwt no-match files).
#!/bin/sh

for i in `awk '{print $1}' SampleSheet_Tobeprocessed.txt`
do
	## Converting btrimmed fastq files into tab-delemited text format
	
	python FastqToTableConverter.py $i/$i.R1.btrim.fastq $i/$i.R1.btrim.fastq.tab.txt
	sed 's/^@//g' $i/$i.R1.btrim.fastq.tab.txt|sort -k 1,1 > $i/$i.R1.btrim.fastq.tab.ToJoin.txt
	
	python FastqToTableConverter.py $i/$i.R2.btrim.fastq $i/$i.R2.btrim.fastq.tab.txt
	sed 's/^@//g' $i/$i.R2.btrim.fastq.tab.txt|sort -k 1,1 > $i/$i.R2.btrim.fastq.tab.ToJoin.txt
	
	## R1 pivot file generation for EGFRvIII sequences
	
	join -1 1 -2 1 $i/EGFR_exon1_exon8 $i/$i.R1.btrim.fastq.tab.ToJoin.txt > $i/EGFR_exon1_exon8_seq.R1.txt
	awk '{print $5}' $i/EGFR_exon1_exon8_seq.R1.txt|uniq -c > $i/EGFR_exon1_exon8_seq.R1.pivot.txt
	grep -f EGFRvIII_ToSearch.txt $i/EGFR_exon1_exon8_seq.R1.pivot.txt > $i/EGFR_exon1_exon8_seq.R1.pivot.final.txt
	grep -v -f EGFRvIII_ToSearch.txt $i/EGFR_exon1_exon8_seq.R1.pivot.txt > $i/EGFR_exon1_exon8_seq.R1.pivot.noMatch.txt
	
	awk '{print $2}' $i/EGFR_exon1_exon8_seq.R1.pivot.txt > $i/EGFR_exon1_exon8_seq.R1.pivot.awked.txt
	python seqToFasta.py $i/EGFR_exon1_exon8_seq.R1.pivot.awked.txt $i/EGFR_exon1_exon8_seq.R1.pivot.fasta
	
	## R2 pivot file generation for EGFRvIII sequences
	
	join -1 1 -2 1 $i/EGFR_exon1_exon8 $i/$i.R2.btrim.fastq.tab.ToJoin.txt > $i/EGFR_exon1_exon8_seq.R2.txt
	awk '{print $5}' $i/EGFR_exon1_exon8_seq.R2.txt > $i/EGFR_exon1_exon8_seq.R2.awked.txt
	python seqToFasta.py $i/EGFR_exon1_exon8_seq.R2.awked.txt $i/EGFR_exon1_exon8_seq.R2.awked.fasta
	perl ReverseComplement.pl $i/EGFR_exon1_exon8_seq.R2.awked.fasta $i/EGFR_exon1_exon8_seq.R2.awked.RC.fasta
	grep -v '>' $i/EGFR_exon1_exon8_seq.R2.awked.RC.fasta > $i/EGFR_exon1_exon8_seq.R2.awked.RC.txt
	uniq -c $i/EGFR_exon1_exon8_seq.R2.awked.RC.txt > $i/EGFR_exon1_exon8_seq.R2.awked.RC.pivot.txt
	grep -f EGFRvIII_ToSearch.txt $i/EGFR_exon1_exon8_seq.R2.awked.RC.pivot.txt > $i/EGFR_exon1_exon8_seq.R2.awked.RC.pivot.final.txt
	grep -v -f EGFRvIII_ToSearch.txt $i/EGFR_exon1_exon8_seq.R2.awked.RC.pivot.txt > $i/EGFR_exon1_exon8_seq.R2.awked.RC.pivot.noMatch.txt
	
	## R1 pivot file generation for EGFRwt sequences
	
	join -1 1 -2 1 $i/EGFR_exon1_exon2 $i/$i.R1.btrim.fastq.tab.ToJoin.txt > $i/EGFR_exon1_exon2_seq.R1.txt
	awk '{print $5}' $i/EGFR_exon1_exon2_seq.R1.txt|uniq -c > $i/EGFR_exon1_exon2_seq.R1.pivot.txt
	grep -f EGFRwt_ToSearch.txt $i/EGFR_exon1_exon2_seq.R1.pivot.txt > $i/EGFR_exon1_exon2_seq.R1.pivot.final.txt
	grep -v -f EGFRwt_ToSearch.txt $i/EGFR_exon1_exon2_seq.R1.pivot.txt > $i/EGFR_exon1_exon2_seq.R1.pivot.noMatch.txt
	
	awk '{print $2}' $i/EGFR_exon1_exon2_seq.R1.pivot.txt > $i/EGFR_exon1_exon2_seq.R1.pivot.awked.txt
	python seqToFasta.py $i/EGFR_exon1_exon2_seq.R1.pivot.awked.txt $i/EGFR_exon1_exon2_seq.R1.pivot.fasta
	
	## R2 pivot file generation for EGFRwt sequences
	
	join -1 1 -2 1 $i/EGFR_exon1_exon2 $i/$i.R2.btrim.fastq.tab.ToJoin.txt > $i/EGFR_exon1_exon2_seq.R2.txt
	awk '{print $5}' $i/EGFR_exon1_exon2_seq.R2.txt > $i/EGFR_exon1_exon2_seq.R2.awked.txt
	python seqToFasta.py $i/EGFR_exon1_exon2_seq.R2.awked.txt $i/EGFR_exon1_exon2_seq.R2.awked.fasta
	perl ReverseComplement.pl $i/EGFR_exon1_exon2_seq.R2.awked.fasta $i/EGFR_exon1_exon2_seq.R2.awked.RC.fasta
	grep -v '>' $i/EGFR_exon1_exon2_seq.R2.awked.RC.fasta > $i/EGFR_exon1_exon2_seq.R2.awked.RC.txt
	uniq -c $i/EGFR_exon1_exon2_seq.R2.awked.RC.txt > $i/EGFR_exon1_exon2_seq.R2.awked.RC.pivot.txt
	grep -f EGFRwt_ToSearch.txt $i/EGFR_exon1_exon2_seq.R2.awked.RC.pivot.txt > $i/EGFR_exon1_exon2_seq.R2.awked.RC.pivot.final.txt
	grep -v -f EGFRwt_ToSearch.txt $i/EGFR_exon1_exon2_seq.R2.awked.RC.pivot.txt > $i/EGFR_exon1_exon2_seq.R2.awked.RC.pivot.noMatch.txt
done

