#!/bin/bash

# example Swift Accel-Amplicon high-sensitivity variant calling analysis workflow

# Original assay: S. Sandhu, S. Chaluvadi and J. Irish 20171002
# Bioinfo script: Akshay Chitturi 5/29/2018 | Last modified 3/7/2019

# USAGE: provide two directories as arguments:
#        1. provide sample directory containing FASTQ pair(s)
#        2. provide SEQ name (appended to end of sample name for FileMaker)

#export DISPLAY=:1.0

# args specified on command line when calling script:
coremaster='/project/cpdlab/SWIFT/cp078_masterfile_171205.txt'
sidmaster='/project/cpdlab/SWIFT/sampleID_masterfile_170113.txt'
fastq_dir="/project/cpdlab/SWIFT/$2/$1"
stats_header='/project/cpdlab/Lymphoma/stats_header.txt'
seq=$2

# start organization
rundirnameroot="$(echo ${coremaster} | tr '_' '\t' | awk '{print $1}')"
runstart="$(date +%y%m%d_%H.%M.%S)"
#rundir="${rundirnameroot}"_"${runstart}"
rundir=$fastq_dir

cd ${fastq_dir}
#mv ../${coremaster} .
#mv ../${sidmaster} .

# start organization
mkdir -p tmp fastq fastqc bed bam vcf metrics gatk lofreq logs

# common paths and script-specific aliases # Tried with reference fastas from UCSC and Broad, undercalling problem persists.
ref='/project/cpdlab/Databases/BWA_ucsc/Homo_sapiens_assembly19broad.fasta'
#ref='/project/cpdlab/Databases/BWA_ucsc/ucsc.hg19.fasta'
#anno='/project/cpdlab/Tools/snpEff/'
bedtools='/home/chitturi/bedtools/bin/bedtools'
samtools='/home/chitturi/samtools-1.9/samtools'

# "alias" variables
java8='/project/cpdlab/Tools/java8/jdk1.8.0_77/bin/java'
picard="${java8} -jar /project/cpdlab/Tools/picard-tools-2.18.2/picard.jar" # apt package lacks PcrMetrics!
gatk="${java8} -jar /project/cpdlab/Tools/GATK_3.6/GenomeAnalysisTK.jar" # using GATK version 3.6

###### code to generate BED files from master files ######

# non-merged target BED from coremaster
awk '{print $1,$2,$3,$4}' OFS="\t" "${coremaster}" | sort -k1,1n -k2,2n > core_nonmerged_targets_temp.bed

# non-merged target BED from sidmaster
awk '{print $1,$2,$3,$4}' OFS="\t" "${sidmaster}" | sort -k1,1n -k2,2n > sID_targets.bed

# merged target BED for core
sort -k1,1 -k2,2n core_nonmerged_targets_temp.bed > core_nonmerged_targets.bed
${bedtools} merge -nms -i core_nonmerged_targets.bed | sed 's/;.*//' | sort -k1,1n -k2,2n > core_merged_targets.bed

### overlap-removed non-merged target BED ###
# get all overlapping regions for target bedfile
# NOTE: awk step isolates overlapping regions and removes full regions
#       (intersectBed outputs all full regions as they self-intersect)
bedtools intersect -a core_nonmerged_targets.bed -b core_nonmerged_targets.bed | awk 'a[$1FS$2FS$3]++' OFS="\t" > core_overlapped_regions.bed

# New edge case
if [ ! -s core_overlapped_regions.bed ]
then
    echo "No overlapping amplicons found, progress with loop."
    cp core_nonmerged_targets.bed core_overlapped_regions.bed
fi

bedtools subtract -a core_nonmerged_targets.bed -b core_overlapped_regions.bed > core_nonmerged_noolaps_targets.bed

mv core_overlapped_regions.bed bed

# handle both old and current master file formats (as of 170217)
colcnts=$(awk 'NR==1{print NF}' ${coremaster})
if [ $colcnts -eq 12 ]
then
    # primer BED (for primer trimming)
    awk '{print $1,$5,$6,$7;print $1,$8,$9,$10}' OFS="\t" "${coremaster}" > core_primers.bed
else
    # primer BED (for primer trimming)
    awk '{print $1,$6,$7,$8;print $1,$9,$10,$11}' OFS="\t" "${coremaster}" > core_primers.bed
fi

# Convert target BED files to 5-column format for analysis workflow
for f in *targets.bed
do
    awk '{print $1,$2,$3,"+",$4}' OFS="\t" "$f" > "${f%%.bed}_5col.bed"
done

mv *targets.bed bed

# variables to make changes easier
bedfile='/project/cpdlab/SWIFT/scripts/core_merged_targets_5col.bed'
nomergebed=core_nonmerged_targets_5col.bed
olapfreebed=core_nonmerged_noolaps_targets_5col.bed
sidbedfile=sID_targets_5col.bed

# create the core + sid masterfile and bed files for ot-checking
cat "${coremaster}" "${sidmaster}" > totalmaster.tmp

awk '{print $1,$2,$3,$4}' OFS="\t" totalmaster.tmp | sort -k1,1n -k2,2n > total_nonmerged_targets_temp.bed

sort -k1,1 -k2,2n total_nonmerged_targets_temp.bed > total_nonmerged_targets.bed
${bedtools} merge -nms -i total_nonmerged_targets.bed | sed 's/;.*//' | sort -k1,1n -k2,2n > total_merged_targets.bed

awk '{print $1,$5,$6,$7;print $1,$8,$9,$10}' OFS="\t" totalmaster.tmp > total_primers.bed

totalmergedbed=total_merged_targets.bed
totalprimers=total_primers.bed

###############################################################################
###### Begin workflow for each pair of FASTQ files in working directory #######
###############################################################################

fastqc='/project/cpdlab/Tools/FastQC/fastqc'

for f in *_R1_001.fastq.gz
do

fq1=$f
fq2=${fq1%%_R1_*}_R2_001.fastq.gz
# changed to test proper .bam naming for FileMaker
prefix=`cut -d'_' -f1 <<< ${fq1%%_R1_001.fastq.gz}`"-SEQ-"$2
pref_var=${prefix}

### run FASTQC ###
echo "running fastqc"
${fastqc} -t 8 $fq1 $fq2
mv *fastqc.zip fastqc

### trim adapters ###
# Illumina adapter trimming setup (Trimmomatic)
trimdir='/project/cpdlab/Tools/Trimmomatic-0.36'

echo "Trimming Illumina adapters"
# NOTE: custom adapter file for Accel-amplicon Illumina adapter trimming
java -Xmx24g -Xms16g -jar ${trimdir}/trimmomatic-0.36.jar PE -threads 6 -trimlog ${prefix}_trimmatic_trimlog.log $fq1 $fq2 ${prefix}_R1_atrimd.fq.gz ${prefix}_unpaired_R1.fq.gz ${prefix}_R2_atrimd.fq.gz ${prefix}_unpaired_R2.fq.gz ILLUMINACLIP:${trimdir}/adapters/TruSeq3-PE-JI.fa:2:30:10 MINLEN:30 2> ${prefix}_01_atrim.log

fqt1=${prefix}_R1_atrimd.fq.gz
fqt2=${prefix}_R2_atrimd.fq.gz

bwa='/project/cpdlab/Tools/BWA/bwa-0.7.17' # Latest version of bwa mem

### align reads ###
echo "Aligning with bwa b37"
${bwa}/bwa mem $ref $fqt1 $fqt2 -U 17 -M -t 24 > ${prefix}_nontrimd.sam 2> ${prefix}_02_bwa.log

### name-sort alignment file
echo "Name-sorting SAM file"
${samtools} sort -@ 12 -n -O SAM "${prefix}_nontrimd.sam" > "${prefix}_nontrimd_namesrtd.sam" 2> "${prefix}_02_2_namesort.log"

primerclip='/project/cpdlab/Tools/primerclip'

### trim primers ###
echo "Trimming primers"
${primerclip}/primerclip ~/totalmaster.tmp ${prefix}_nontrimd_namesrtd.sam ${prefix}_ptrimd.sam 2> ${prefix}_03_ptrim.log

#${samtools} sort -@ 24 -O sam "${prefix}_ptrimd.sam" > "${prefix}_ptrimd_srtd.sam" 2> ${prefix}_ptrim_coordsrt.log

echo "sorting and adding read groups"
${picard} AddOrReplaceReadGroups I=${prefix}_ptrimd.sam O=${prefix}_RG_b37.bam SO=coordinate RGID=snpID LB=swift SM=${prefix} PL=illumina PU=miseq VALIDATION_STRINGENCY=LENIENT 2> ${prefix}_04_addRGs.log

# save non-primertrimd BAM file for debugging and inspection
${picard} SortSam I=${prefix}_nontrimd.sam O=${prefix}_nontrimd.bam CREATE_INDEX=true SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT 2> ${prefix}_05_makenonptrimdbam.log

echo "indexing bam file"
${samtools} index ${prefix}_RG_b37.bam 2> ${prefix}_06_index.log

### calculate coverage metrics ###
echo "calculating coverage metrics"
/project/cpdlab/anaconda2/bin/bedtools coverage -abam ${prefix}_RG_b37.bam -b $bedfile -d > ${prefix}.covd # bedtools v2.26.0
awk '{sum+=$7}END{m=(sum/NR); b=m*0.2; print m, b}' ${prefix}.covd > ${prefix}_covd.tmp 2> ${prefix}_06_cov1.log
awk 'BEGIN{n=0}NR==FNR{m=$1;b=$2;next}{if($7>=b)n++}END{print m,b,(n/FNR*100.0)}' OFS="\t" ${prefix}_covd.tmp ${prefix}.covd > ${prefix}_covMetrics.txt
# Adds a third number, but it does generate...

# make intervals file for CollectTargetedPcrMetrics
${samtools} view -H ${prefix}_RG_b37.bam > ${prefix}_header.txt
# concatenate core and sampleID target BED files for on-target
cat ${bedfile} ${sidbedfile} | sort -k1,1n -k2,2n > core_sID_concatd_targets_5col.bed

ontargetbed='core_sID_concatd_targets_5col.bed'
cat ${prefix}_header.txt ${ontargetbed} > ${prefix}_fullintervals
cat ${prefix}_header.txt ${ontargetbed} > ${prefix}_noprimerintervals

# find on-target metrics using picard-tools
echo "Running CollectTargetedPcrMetrics"
${picard} CollectTargetedPcrMetrics I=${prefix}_RG_b37.bam O=${prefix}_targetPCRmetrics.txt AI=${prefix}_fullintervals TI=${prefix}_noprimerintervals R=$ref PER_TARGET_COVERAGE=${prefix}_perTargetCov.txt VALIDATION_STRINGENCY=LENIENT 2> ${prefix}_09_pcrmetrics.log

#rm ${prefix}_*intervals ${prefix}_header.txt
mv ${ontargetbed} bed

# sampleID spike-in coverage metrics
echo "calculating coverage metrics for sampleID spike-in"
/project/cpdlab/anaconda2/bin/bedtools coverage -abam ${prefix}_RG_b37.bam -b $sidbedfile -d > ${prefix}_sID.covd
awk '{sum+=$7}END{m=(sum/NR); b=m*0.2; print m, b}' ${prefix}_sID.covd > ${prefix}_sID_covd.tmp 2> ${prefix}_sID_10_cov1.log

awk 'BEGIN{n=0}NR==FNR{m=$1;b=$2;next}{if($7>=b)n++}END{print m,b,(n/FNR*100.0)}' OFS="\t" ${prefix}_sID_covd.tmp ${prefix}_sID.covd > ${prefix}_sID_covMetrics.txt

###############################################################################
########################  Variant Calling  ####################################
###############################################################################

echo "Starting variant calling with GATK"
# NEW 171002 add back in alignment post-processing before variant calling
db_snp='/project/cpdlab/Databases/dbsnp/dbsnp141.vcf'

${gatk} -T RealignerTargetCreator -R $ref -I ${prefix}_RG_b37.bam -o ${prefix}_indelRealign.intervals 2> ${prefix}_gatk_realgntargcreat.log

${gatk} -T IndelRealigner -R $ref -I ${prefix}_RG_b37.bam --targetIntervals ${prefix}_indelRealign.intervals -o ${prefix}_realigned.bam 2> ${prefix}_gatk_indelrealgn.log

${gatk} -T BaseRecalibrator -R $ref -I ${prefix}_realigned.bam --knownSites ${db_snp} -nct 12 -o ${prefix}_recal_data_table.txt 2> ${prefix}_gatk_baserecal.log

${gatk} -T PrintReads -R $ref -I ${prefix}_realigned.bam -BQSR ${prefix}_recal_data_table.txt -o ${prefix}_realigned_bqsrCal.bam 2> ${prefix}_gatk_printrds.log

#gatk HC
${gatk} -T HaplotypeCaller -R $ref -I ${prefix}_realigned_bqsrCal.bam -stand_call_conf 20 --dontUseSoftClippedBases -stand_emit_conf 20 -mbq 20 -L $bedfile -o ${prefix}_gatkHC.vcf 2> ${prefix}_gatkHC.log

#low-frequency variant calls
echo "Starting variant calling with LoFreq"
lofreq='/project/cpdlab/anaconda2/bin/lofreq_2'
${lofreq} call -q 20 -Q 20 -m 30 -C 50 --call-indels -f $ref ${prefix}_realigned_bqsrCal.bam -l $bedfile -o ${prefix}_notmerged_lf.vcf 2> ${prefix}_lofreq.log

# Skipping merging step for now, variant review wants to see duplicate calls as well.
cp ${prefix}_notmerged_lf.vcf ${prefix}_lf_b37.vcf

###############################################################################
####################  sampleID variant calling  ###############################
###############################################################################

echo "Starting sampleID spike-in variant calling with GATK"
#gatk HC
${gatk} -T HaplotypeCaller -R $ref -I ${prefix}_RG_b37.bam -stand_call_conf 20 -mbq 20 -L $sidbedfile -o ${prefix}_gatkHC_sID.vcf 2> ${prefix}_gatkHC_sID.log

#low freq variant calls
echo "Starting sampleID spike-in variant calling with LoFreq"
${lofreq} call -q 20 -Q 20 -m 30 -C 50 --call-indels -f $ref ${prefix}_RG_b37.bam -l $sidbedfile -o ${prefix}_lf_b37_sID.vcf 2> ${prefix}_lofreq_sID.log

#annotate LoFreq VCFs with snpEff, Oncotator
snpeff='/project/cpdlab/Tools/snpEff/4.1l/snpEff.jar'
snpeff_conf='/project/cpdlab/Tools/snpEff/4.1l/snpEff.config'
java7='/project/cpdlab/Tools/jdk1.7.0_80/bin/java'
Oncotator='/project/cpdlab/Tools/oncotator/oncotator-1.8.0.0/oncotator/Oncotator.py'
db_onco='/project/cpdlab/Databases/oncotator_db/oncotator_v1_ds_April052016'
stats_header='/project/cpdlab/Lymphoma/stats_header.txt'
onco_parse='/project/cpdlab/Scripts/SolidV2/solidv2_v1.2/methods/onco_parse.oldswift.py'

#------------------------------------------------
#--- Annotate VCF using 4% limit of detection ---
#------------------------------------------------

# Annotate VCF with SnpEff, Oncotator
${java7} -jar ${snpeff} -c ${snpeff_conf} -canon hg19 -i vcf -noStats ${prefix}_lf_b37.vcf > ${prefix}_combined_snpeff.vcf
/project/cpdlab/anaconda2/bin/python ${Oncotator} -v -i VCF --db-dir ${db_onco} -o TCGAMAF ${prefix}_combined_snpeff.vcf ${prefix}.anno.TCGAMAF hg19
/project/cpdlab/anaconda2/bin/python ${onco_parse} ${prefix}.anno.TCGAMAF ${prefix}.parsed.tsv ${pref_var}

# End of loop, master variant table creation occurs outside.
done

# Master variant table creation
variants=$1_Run_masterVarFinal.txt
awk '!/#/' *parsed.tsv >> ${variants}
#head -n1 ${prefix}.parsed.tsv | cat - ${variants} > temp && mv temp ${variants}

# Summarize on-target and coverage metrics for all samples into a single report
for f in *targetPCRmetrics.txt
do
    #awk -v n=${f%%_L001*} 'NR==1{print n,$11,$14*100.0,$21*100.0}' OFS="\t" | (tail -n 3 $f) > ${f%%.txt}_summary.txt
    awk -v n=${f%%_L001*} 'NR==8{print n,$11,$14*100.0,$22*100.0}' OFS="\t" $f > ${f%%.txt}_summary.txt
    f2=${f%%_target*}_covMetrics.txt
    paste ${f%%.txt}_summary.txt $f2 > ${f2%%_cov*}_combined_cov_metrics.txt
done

#----------------------------------------------------------------------------------
#----- Parse all results to match format of Penn Hospital's FileMaker system. -----
#----------------------------------------------------------------------------------

# Generate RunStatistics
stats=$1_RunStatsFinal.txt
#cp ${stats_header} ${stats}

for line in $(ls *R1_001.fastq.gz)
do
  sample=`cut -d'_' -f1 <<< ${line%%_R1_001.fastq.gz}`
  #sample=${line%????????????????} 
  total_reads=$(${samtools} view -c ${sample}*nontrimd.sam)
  echo "Total reads: "${total_reads}
  reads_for_alnment=$(${samtools} flagstat ${sample}*nontrimd.sam | awk 'NR==6 {print $1}')
  per_lost=$(echo "scale=4; 100-($reads_for_alnment/$total_reads)*100" | bc)
  mapped_reads=$(${samtools} view -c -F 260 ${sample}*_b37.bam)
  per_mapped=$(echo "scale=4; ($mapped_reads/$total_reads)*100" | bc)
  on_target_reads=$(${samtools} view -c -F 260 -L $bedfile ${sample}*nontrimd.sam)
  per_ontarget=$(echo "scale=4; ($on_target_reads/$total_reads)*100" | bc)
  reads_post_filter=$on_target_reads
  per_after_qual="NA"
  per_useable=$(echo "scale=4; ($reads_post_filter/$total_reads)*100" | bc)
  mean_cov=$(${samtools} depth -m 1000000 -b $bedfile ${sample}*realigned.bam | awk '{sum+=$3} END {print sum/NR}')

  # coverage-array loops

  cov_array=$(${samtools} depth -m 1000000 -b $bedfile ${sample}*realigned.bam | awk '{print $3}')
  per_0=$(echo "scale=4; (1-$(echo "$cov_array" | awk '$0<1' | wc -l)/$(echo "$cov_array" | wc -l))*100" | bc)
  per_1=$(echo "scale=4; (1-$(echo "$cov_array" | awk '$0<2' | wc -l)/$(echo "$cov_array" | wc -l))*100" | bc)
  per_250=$(echo "scale=4; (1-$(echo "$cov_array" | awk '$0<250' | wc -l)/$(echo "$cov_array" | wc -l))*100" | bc)
  per_1K=$(echo "scale=4; (1-$(echo "$cov_array" | awk '$0<1000' | wc -l)/$(echo "$cov_array" | wc -l))*100" | bc)
  
  # No. of exonic regions below cutoff

  ${gatk} -T DepthOfCoverage -R $ref -I ${sample}*realigned.bam --countType COUNT_FRAGMENTS -U ALLOW_SEQ_DICT_INCOMPATIBILITY -o ${sample}_output -L /home/chitturi/swift_exonic_regions_final.bed
  ex_counter250=`cat ${sample}_output.sample_interval_summary | awk '$3<250' | wc -l`
  ex_counter150=`cat ${sample}_output.sample_interval_summary | awk '$3<150' | wc -l`

  clip_count="NA"
  #fname=${sample%?????????} de-limit sample instead!
  fname=$(cut -d'_' -f1 <<< $sample)"-SEQ-"$2

  echo -e "$fname\t$total_reads\t$reads_for_alnment\t$per_lost\t$mapped_reads\t$per_mapped\t$on_target_reads\t$per_ontarget\t$reads_post_filter\t$per_after_qual\t$per_useable\t$mean_cov\t$per_0\t$per_1\t$per_250\t$per_1K\t$ex_counter250\t$clip_count\t$ex_counter150" >> ${stats}
done

echo "Sample Reads_Aligned %Reads_Aligned %Aligned_Bases_OnTarget" "Mean_Coverage 20%Mean_Coverage %Coverage_Uniformity" | tr ' ' '\t' > final_metrics_report.txt
cat *_combined_cov_metrics.txt >> final_metrics_report.txt

mv final_metrics_report.txt metrics
rm *sam
mv *.log logs
mv *.bed bed
mv *.ba* bam
rm *.covd
mv *covMetrics.txt metrics
mv *targetPCRmetrics.txt metrics
mv *.fastq.gz fastq
mv *.fq.gz fastq
mv *combined_cov_metrics.txt metrics
mv *PCRmetrics_summary.txt metrics
mv *perTargetCov.txt metrics
mv *.vcf vcf
mv *.idx vcf
cat ${stats} >> ../RunStatsFinal.txt
cat ${variants} >> ../Run_masterVarFinal.txt
cd ..
#rsync -avPR ${fastq_dir} cpdlab@170.212.141.107:/PathCPD/FromHPC/$2
echo " "
echo "Swift 60G analysis workflow finished. Results on Isilon ready to upload."
echo " "
