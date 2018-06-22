#!/bin/bash

# Note: This script is typically run from parent "SEQ_XXXXXX_run_this.sh" in /project/cpdlab/Lymphoma/
# snpeff VCF --> Oncotator --> anno/anno_in (MAF) --> modified onco_parse --> parsed vcf for FileMaker

# Dependencies
snpeff='/project/cpdlab/Tools/snpEff/4.1l/snpEff.jar'
snpeff_conf='/project/cpdlab/Tools/snpEff/4.1l/snpEff.config'
java7='/project/cpdlab/Tools/jdk1.7.0_80/bin/java'
onco_parse='/project/cpdlab/Scripts/SolidV2/solidv2_v1.2/methods/onco_parse.lymphoma.py'
Oncotator='/project/cpdlab/Tools/oncotator/oncotator-1.8.0.0/oncotator/Oncotator.py'
db_onco='/project/cpdlab/Databases/oncotator_db/oncotator_v1_ds_April052016'
exon_bed='/project/cpdlab/Lymphoma/lymphoma.merge.bed' 
stats_header='/project/cpdlab/Lymphoma/stats_header.txt'

# Inputs
vcfs=$1
dir=$2

# Check for correct input (if running script directly)
if [ $# -eq 0 ] || [ -z "$1" ]; then
  echo "Two arguments must be supplied:"
  echo " "
  echo "ARG1 = full run name, found in /PathCPD/illumina/SeqOut/"
  echo "ARG2 = run ID, e.g. SEQ_XXXXXX"
  echo " "
  echo "Note: Run this in interactive session or via a non-interactive job."
  echo " "
  exit 1
fi

#remove genome vcfs from processing
rm ${vcfs}/*genome.vcf

# Annotate VCFs with Oncotator, parse into FileMaker-recognizable format.
for file in ${vcfs}/*.vcf
do
  ${java7} -jar ${snpeff} -c ${snpeff_conf} -canon hg19 -i vcf -noStats $file > $file.snpeff.lymphoma.vcf
  python ${Oncotator} -v -i VCF --db-dir ${db_onco} -o TCGAMAF $file.snpeff.lymphoma.vcf $file.anno.TCGAMAF hg19
  python ${onco_parse} $file.anno.TCGAMAF $file.parsed.tsv ${dir}
 echo $file
done

# Generate Run_masterVarFinal table for FileMaker
variants=${dir}_Run_masterVarFinal.txt
echo "Adding variants for "$file
awk '!/#/' ${vcfs}/*parsed.tsv >> ${variants}
head -n1 $file.parsed.tsv | cat - ${variants} > temp && mv temp ${variants}
sed -i 's/TSCA/LYMPH/g' ${variants}

# Generate RunStatsFinal table for FileMaker
stats=${dir}_RunStatsFinal.txt
cp ${stats_header} ${vcfs}/${stats} # write loop results to $stats from here forward
cd ${vcfs}
count=0
declare -a array=(A B A B A B A B A B A B A B A B A B A B A B A B A B A B A B A B A B A B A B A B) # A/B marker (Manifest column in SampleSheet)
cp SampleSheet.csv SampleSheet.cp.csv
sed -i 's/TSCA/LYMPH/g' SampleSheet.csv AmpliconCoverage_M1.tsv AmpliconCoverage_M2.tsv

for line in $(ls *.bam | sed 's/\..*//' | sort -t _ -u -k 1,1 -k 2.2n) # sorts by sample, then by SampleSheet position (A always comes first)
do
  cd FASTQs/
  total_reads=$(($(zcat $line"_L001_R1_001".fastq.gz | echo $((`wc -l`/4))) + $(zcat $line"_L001_R2_001".fastq.gz | echo $((`wc -l`/4))))) # manual read count
  cd ..
  total_aligned_reads="NA"
  per_lost="NA"
  mapped_reads=$(samtools view -c -F 260 $line.bam)
  per_mapped=$(echo "scale=4; ($mapped_reads/$total_reads)*100" | bc)
  true_name=`cut -d'_' -f1 <<< $line | sed 's/-/_/g'` # Sample name value in SampleSheet
  # Finds full sample ID by matching sample name and manifest
  name=$(awk -v var1=$true_name 'BEGIN {FS=","} $2==var1' SampleSheet.csv | awk -v var2=${array[count]} 'BEGIN {FS=","} $9==var2 {print $1}')
  echo "Collecting statistics for "$name" ("${array[count]}")"
    
  if [[ ${array[count]} == A ]];
  then
    sum=0 
    reads1=$(cut AmpliconCoverage_M1.tsv -f `head -1 AmpliconCoverage_M1.tsv | tr "\t" "\n" | grep -n "^$name"'$' | cut -f 1 -d :`)
    for num in ${reads1};
    do
      sum=$((sum+num))
    done
  fi
  
  if [[ ${array[count]} == B ]];
  then
    sum=0
    reads2=$(cut AmpliconCoverage_M2.tsv -f `head -1 AmpliconCoverage_M2.tsv | tr "\t" "\n" | grep -n "^$name"'$' | cut -f 1 -d :`)
    for num in ${reads2};
    do
      sum=$((sum+num))
    done
  fi
  sum_reads=$sum

  on_target_reads="NA"
  per_on_target="NA"
  qual_filter_reads=$sum_reads
  per_after_qual=$(echo "scale=4; ($qual_filter_reads/$mapped_reads)*100" | bc)
  per_usable=$(echo "scale=4; ($qual_filter_reads/$total_reads)*100" | bc)
 
  samtools depth $line.bam > $line.coverage.txt
  mean_cov=$(awk '{sum+=$3} END { print sum/NR}' $line.coverage.txt)
  cov_array=$(awk '{print $3}' $line.coverage.txt)
  
  counter1=0
  counter250=0
  counter1K=0
  countertotal=0
  for entry250 in ${cov_array};
  do
    countertotal=$((countertotal+1))
    if [ "$entry250" -lt 250 ]; then
      counter250=$((counter250+1))
    fi
  done

  for entry1K in ${cov_array};
  do
    if [ "$entry1K" -lt 1000 ]; then
      counter1K=$((counter1K+1))
    fi
  done

  for entry in ${cov_array};
  do
    if [ "$entry" -lt 1 ]; then
      counter1=$((counter1+1))
    fi
  done

  per_0=1
  per_1=$(echo "scale=4; (1-($counter1/$countertotal))*100" | bc)
  per_250=$(echo "scale=4; (1-($counter250/$countertotal))*100" | bc)
  per_1K=$(echo "scale=4; (1-($counter1K/$countertotal))*100" | bc)

  bedtools intersect -a $line.bam -b ${exon_bed} > $line.int.bam
  bedtools genomecov -trackline -bg -ibam $line.int.bam > $line.int.bedgraph
  
  ex_counter250=`cat $line.int.bedgraph | awk '$4<250' | wc -l`
  ex_counter150=`cat $line.int.bedgraph | awk '$4<150' | wc -l`
  clip_count="NA"
  
  echo -e "$name"-"$dir\t$total_reads\t$total_aligned_reads\t$per_lost\t$mapped_reads\t$per_mapped\t$on_target_reads\t$per_on_target\t$qual_filter_reads\t$per_after_qual\t$per_usable\t$mean_cov\t$per_0\t$per_1\t$per_250\t$per_1K\t$ex_counter250\t$clip_count\t$ex_counter150" >> ${stats}
  
  # Long one-liner for dummy line variant insertion here, perhaps? Something like...
  echo -e "$name"-"$dir\tchr1\t1\t-\t-\t-\tNA\t0\t0\tNA\t0\tXXXX\tNA\tXXXX\tNA\t\tXXXX\tNA\t\tNA\tNA\tNA\tNA\tNA\t\tc.-\t\tXXXX\tXXXX\tXXXX\tXXXX\t\t\tchr1\t1\t\t-\t-\t0\t\tNA\tNA\tNA\t\t\t\t\t\t\t\tXXX\t0\t0\t0\t0\tXXXX\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t0\t0\tXXXX\t\t\tXXXX\t\tNA\t\t\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t\t\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t\tNA\tNA\tNA\tNA\tNA\tNA" >> ../${variants}
  mv $line.bam $name"-"$dir.bam
  mv $line.bam.bai $name"-"$dir.bam.bai

  cat $line.int.bedgraph | awk '$4<250' > $line.exonic_reg_bel250.bed
  count=$((count+1))
done
cd ..

mv ${variants} ${vcfs}/${variants}
cd ${vcfs}
rm *coverage.bed
rm *int.bam
rm -rf FASTQs/
cd ..

rsync -avP ${vcfs} cpdlab@170.212.141.107:/PathCPD/FromHPC/

rm oncotator.log
rm Oncotator_profile.bin
rm profile_stats.txt
rm ${dir}_run_this.sh
