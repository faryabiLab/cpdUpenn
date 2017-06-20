#!/bin/sh

#  exome_seq_v1.sh
#
#
#  Created by babak on 6/2/16.
#

################## Data Bases ########
## Databases
db_fa='/project/cpdlab/ashkan/muTect/ucsc.hg19.fasta'
db_snp='/project/cpdlab/ashkan/muTect/dbsnp_137.hg19.vcf'
db_cosmic='/project/cpdlab/ashkan/muTect/cosmic_v67.hg19.vcf'
db_onco='/project/cpdlab/Databases/oncotator_db/oncotator_v1_ds_April052016'

## Languages
# Java
java6='java'
java7='/project/cpdlab/Tools/jdk1.7.0_80/bin/java'
java8='/project/cpdlab/Tools/java8/jdk1.8.0_77/bin/java'

###################################
# Trim Galore - 0.4.1
trim_galore='/project/cpdlab/Tools/trim_galore/trim_galore'

# Cut Adapt - 1.3
cut_adapt='/project/cpdlab/Tools/cutadapt/1.3/bin/cutadapt'

#bwa
bwa='/project/cpdlab/Tools/BWA/bwa-0.7.15'

# Samtools - 1.2
samtools='/project/cpdlab/Tools/samtools-1.3.1/samtools'

# Bedtools - 2.25
bedtools='/project/cpdlab/Tools/bedtools/2.25.0/bin/'

# VCF Tools
vcftools='/project/cpdlab/Tools/vcftools/0.1.14/src/perl/'
bcftools='/project/cpdlab/Tools/bcftools-1.3.1'
bgzip='/project/cpdlab/Tools/samtools-1.3.1/htslib-1.3.1'
tabix='/project/cpdlab/Tools/samtools-1.3.1/htslib-1.3.1'

# Picard -  2.4.1
picard='/project/cpdlab/Tools/picard-tools-2.4.1'

# GATK
GATK='/project/cpdlab/Tools/GATK/3.4-46/GenomeAnalysisTK.jar'
GATK2='/project/cpdlab/Tools/GATK_3.5/GenomeAnalysisTK.jar'
GATKkey='/project/cpdlab/Tools/GATK/3.4-46/CPD_upenn_gatk.key'

# AgilentDedup - 1
MBCdedup='/project/cpdlab/Tools/agilent_dedup/AgilentMBCDedup.jar'

#variant_callers
freebayes='/project/cpdlab/Tools/freebayes/bin/freebayes'
varscan2='/project/cpdlab/Tools/varscan-master/VarScan.v2.4.1.jar'
vardict='/project/cpdlab/Tools/VarDictJava/VarDict/'
scalpel='/project/cpdlab/Tools/scalpel-0.5.3'

#annotations
alamut='/project/cpdlab/Tools/alamut/1.4.2/./alamut-batch'
oncotator='/project/cpdlab/Tools/oncotator/oncotator-1.8.0.0/oncotator'
#snpEff/SnpSift
snpeff='/project/cpdlab/Tools/snpEff/4.1l/snpEff.jar'
snpeff_conf='/project/cpdlab/Tools/snpEff/4.1l/snpEff.config'

#custom parsers
coverage_calc='/project/cpdlab/HiSeqRun/coverage_summary.py'
onco_parse='/project/cpdlab/HiSeqRun/onco_parse.py'
sample_stats='/project/cpdlab/HiSeqRun/sample_stat.py'
amp_calc='/project/cpdlab/HiSeqRun/amp_calc.py'

###################################
out_dir=$1

cd ${out_dir}

SampleName=$(awk -F',' 'NR==2 {print $3}' SampleSheet.csv)
read_index=$(awk -F',' 'NR==2 {print $5}' SampleSheet.csv)


read1=${out_dir}/${SampleName}.R1.fastq.gz
read2=${out_dir}/${SampleName}.R3.fastq.gz
index2=${out_dir}/${SampleName}.R2.fastq.gz
tmp='/project/cpdlab/HiSeqRun/tmp'

cpu=24
mem=72
adapter1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
adapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
amplicon_bed='/project/cpdlab/cpdUpenn/solid2_dev/bed_files/04818-1457701567_Amplicons.bed'
target_bed='/project/cpdlab/cpdUpenn/solid2_dev/bed_files/04818-1457701567_Regions_clean.bed'
coverage_bed='/project/cpdlab/cpdUpenn/solid2_dev/bed_files/04818-1457701567_Covered.bed'
hotspots_WG_bed='/project/cpdlab/cpdUpenn/solid2_dev/bed_files/solid2_hotspots_whole_gene.bed'
hotspots_target_bed='/project/cpdlab/cpdUpenn/solid2_dev/bed_files/solidv2_hotspot_regions.bed'
normal_amps='/project/cpdlab/cpdUpenn/solid2_dev/bed_files/s2_normals_amp.tsv'



FW=$(basename $read1 ".fastq.gz")
BW=$(basename $read2 ".fastq.gz")
B=$(basename $read1 ".R1.fastq.gz")


${trim_galore} -q 20 --phred33 --fastqc -a ${adapter1} -a2 ${adapter2} --stringency 3 -e 0.1 --length 20 --paired ${read1} ${read2} -o ${out_dir} --path_to_cutadapt ${cut_adapt}

#### Align ####
trim_fq1=${FW}_val_1.fq.gz
trim_fq2=${BW}_val_2.fq.gz
out_sam=${B}.align.sam
sam_stat=${B}.sam.stat

${bwa}/bwa mem -t ${cpu} -P ${db_fa} ${trim_fq1} ${trim_fq2}  > ${out_sam}

${samtools} flagstat ${out_sam} > ${sam_stat}

##### Deduplicate (should be replaced by Picard removeDuplicates) #######

dedup_out=${B}.nodup.bam
dedup_fix_out=${B}.nodup.fix.bam
dedup_fix_sort=${B}.nodup.fix.srt.bam
dedup_out_stat=${B}.nodup.stat
dedup_fix_sort_vld=${B}.nodup.fix.srt.vr

${java8} -Xmx${mem}g -jar ${MBCdedup} -X ${tmp} -b ${amplicon_bed} -o ${dedup_out} ${out_sam} ${index2}

## fix the bam file ##
${java8} -Xmx${mem}g -jar ${picard}/picard.jar AddOrReplaceReadGroups TMP_DIR=${tmp} I=${dedup_out} O=${dedup_fix_out} RGID="halo" RGLB=${B} RGPL=Illumina RGPU=${read_index} RGSM=${B}

${samtools} sort -m 32000000000 ${dedup_fix_out} -o ${dedup_fix_sort}
${samtools} index ${dedup_fix_sort}

#${samtools} flagstat ${dedup_fix_sort} > ${dedup_out_stat}

### for the picard dedup ####

#out_bam=${B}.align.bam
#${samtools} view -Su ${out_sam} | ${samtools} sort -m 30000000000 - -o ${out_bam}

#dedup_out=${B}.nodup.bam
#dedup_fix_out=${B}.nodup.fix.bam
#dedup_out_sort=${B}.nodup.srt.bam
#dup_file_qc=${B}.dup.qc
#dedup_out_stat=${B}.nodup.stat
#${java8} -Xmx${mem}g -jar ${picard}/picard.jar MarkDuplicates TMP_DIR=${tmp} M=${dup_file_qc} I=${out_bam} O=${dedup_out} REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT

## fix the bam file ##
#${java8} -Xmx${mem}g -jar ${picard}/picard.jar AddOrReplaceReadGroups TMP_DIR=${tmp} I=${dedup_out} O=${dedup_fix_out} RGID="exome-seq" RGLB=${B} RGPL=Illumina RGPU=${read_index} RGSM=${B}


#${samtools} sort -m 32000000000 ${dedup_fix_out} -o ${dedup_fix_sort}
#${samtools} index ${dedup_fix_sort}

#${samtools} flagstat ${dedup_fix_sort} > ${dedup_out_stat}

############################################################################

${java8} -Xmx${mem}g -jar ${picard}/picard.jar ValidateSamFile I=${dedup_fix_sort} O=${dedup_fix_sort_vld}

############################################################################


############################################################################
##

dedup_fix_clean=${B}.nodup.fix.clean.bam
dedup_fix_clean_sort=${B}.nodup.fix.clean.srt.bam

${java8} -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${GATK2} -T PrintReads -I ${dedup_fix_sort} -o ${dedup_fix_clean} -R ${db_fa} -nct ${cpu} -L ${amplicon_bed}

${samtools} sort -m 32000000000 ${dedup_fix_clean} -o ${dedup_fix_clean_sort}
${samtools} index ${dedup_fix_clean_sort}

################################## Realigner #################################
realign_intv=${B}.realign.intervals
realigned_bam=${B}.realign.bam
realigned_bam_sort=${B}.realign.srt.bam

# Realignment target creator
${java8} -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${GATK2} -T RealignerTargetCreator -R ${db_fa} -I ${dedup_fix_clean_sort} -o ${realign_intv} -L ${target_bed} -known /project/cpdlab/ashkan/muTect/Mills_and_1000G_gold_standard.indels.ucsc.hg19.sites.vcf -known /project/cpdlab/Databases/GATK/bundle/2.8/hg19/1000G_phase1.indels.hg19.sites.vcf

# realignment
${java8} -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${GATK2} -T IndelRealigner -R ${db_fa} -L ${target_bed} -I ${dedup_fix_clean_sort} -known /project/cpdlab/ashkan/muTect/Mills_and_1000G_gold_standard.indels.ucsc.hg19.sites.vcf -known /project/cpdlab/Databases/GATK/bundle/2.8/hg19/1000G_phase1.indels.hg19.sites.vcf -targetIntervals ${realign_intv} -o ${realigned_bam}

${samtools} sort -m 32000000000 ${realigned_bam} -o ${realigned_bam_sort}
${samtools} index ${realigned_bam_sort}

################################## Base Recallibration #################################

recal_table=${B}.grp
recal_bam=${B}.realign.recal.bam
recal_sort_bam=${B}.realign.recal.srt.bam
after_recal_table=${B}.after_recal.grp
recal_plot=${B}_recal.pdf

${java8} -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${GATK2} -T BaseRecalibrator -R ${db_fa} -nct ${cpu} -knownSites /project/cpdlab/ashkan/muTect/Mills_and_1000G_gold_standard.indels.ucsc.hg19.sites.vcf -knownSites /project/cpdlab/Databases/GATK/bundle/2.8/hg19/1000G_phase1.indels.hg19.sites.vcf -knownSites /project/cpdlab/Databases/GATK/bundle/2.8/hg19/dbsnp_138.hg19.vcf -L ${target_bed} -I ${realigned_bam_sort} -o ${recal_table}

${java8} -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${GATK2} -T PrintReads -R ${db_fa} -L ${target_bed} -nct ${cpu} -I ${realigned_bam_sort} -BQSR ${recal_table} -o ${recal_bam}

${samtools} sort -m 32000000000 ${recal_bam} -o ${recal_sort_bam}
${samtools} index ${recal_sort_bam}

# plot recallibration effect

${java8} -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${GATK2} -T BaseRecalibrator -R ${db_fa} -nct ${cpu} -knownSites /project/cpdlab/ashkan/muTect/Mills_and_1000G_gold_standard.indels.ucsc.hg19.sites.vcf -knownSites /project/cpdlab/Databases/GATK/bundle/2.8/hg19/1000G_phase1.indels.hg19.sites.vcf -knownSites /project/cpdlab/Databases/GATK/bundle/2.8/hg19/dbsnp_138.hg19.vcf -L ${target_bed} -I ${realigned_bam_sort} -BQSR ${recal_table} -o ${after_recal_table}

${java8} -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${GATK2} -T AnalyzeCovariates -R ${db_fa} -before ${recal_table} -after ${after_recal_table} -plots ${recal_plot}

################################# rename file #######################################
final_bam=${B}_final.bam
final_bam_stat=${B}_final.stat

cp ${recal_sort_bam} ${final_bam}
${samtools} index ${final_bam}

${samtools} flagstat ${final_bam} > ${final_bam_stat}

################################# stats  #######################################
final_bam=${B}_final.bam
target_depth_out=${B}.target.depth
coverage_depth_out=${B}.coverage.depth
hotspots_target_out=${B}.hotspots.depth
${java8} -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${GATK2} -T DepthOfCoverage -R ${db_fa} -K ${GATKkey} -et NO_ET -I ${final_bam} -o ${target_depth_out} --minBaseQuality 22 -baseCounts -ct 0 -ct 1 -ct 150 -ct 250 -ct 1000 -L ${target_bed}
${java8} -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${GATK2} -T DepthOfCoverage -R ${db_fa} -K ${GATKkey} -et NO_ET -I ${final_bam} -o ${coverage_depth_out} --minBaseQuality 22 -baseCounts -ct 0 -ct 1 -ct 150 -ct 250 -ct 1000 -L ${coverage_bed}
${java8} -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${GATK2} -T DepthOfCoverage -R ${db_fa} -K ${GATKkey} -et NO_ET -I ${final_bam} -o ${hotspots_target_out} --minBaseQuality 22 -baseCounts -ct 0 -ct 1 -ct 150 -ct 250 -ct 1000 -L ${hotspots_target_bed}

python ${coverage_calc} ${hotspots_target_out}
python ${coverage_calc} ${target_depth_out}
python ${coverage_calc} ${coverage_depth_out}

#################################### Variant Calls ####################################
mutect_vcf=${B}.mutect.vcf
mutect_vcf_gz=${B}.mutect.vcf.gz

${java8} -Xmx${mem}g -Djava.io.tmpdir=${tmp} -jar ${GATK2} -nct ${cpu} -T MuTect2 -R ${db_fa} -I:tumor ${final_bam} -L ${target_bed} -o ${mutect_vcf}

sed -i -e 's/##fileformat=VCFv4.2/##fileformat=VCF4.1/g' ${mutect_vcf}
sed -i '/PASS\|##\|#/!d' ${mutect_vcf}
${bgzip}/bgzip -c ${mutect_vcf} > ${mutect_vcf_gz}
${tabix}/tabix -p vcf ${mutect_vcf_gz}

##freebayes call DOCS: https://github.com/ekg/freebayes
freebayes_vcf=${B}.freebayes.vcf
freebayes_vcf_gz=${B}.freebayes.vcf.gz

${freebayes} -f ${db_fa} ${final_bam} --min-mapping-quality 10  --min-alternate-count 2 --min-alternate-fraction 0.01 > ${freebayes_vcf}
${bgzip}/bgzip -c ${freebayes_vcf} > ${freebayes_vcf_gz}
${tabix}/tabix -p vcf ${freebayes_vcf_gz}


#varscan2 call DOCS: http://varscan.sourceforge.net/using-varscan.html
#requires pile up
final_pile=${B}.piled
varscan2_snp_vcf=${B}.varscan2_snp.vcf
varscan2_indel_vcf=${B}.varscan2_indel.vcf
varscan2_combined_vcf=${B}.varscan2_combined.vcf
varscan2_combined_sort_vcf=${B}.varscan2_combined_sort.vcf
varscan2_combined_sort_vcf_gz=${B}.varscan2_combined_sort.vcf.gz

${samtools} mpileup -A -C50 -f ${db_fa} ${final_bam} -l ${target_bed} > ${final_pile}
${java7} -jar ${varscan2} mpileup2snp ${final_pile} --output-vcf 1 --p-value 0.01 --min-coverage 10 --minreads2 0 --min-var-freq 0.01 > ${varscan2_snp_vcf}
${java7} -jar ${varscan2} mpileup2indel ${final_pile} --output-vcf 1 --p-value 0.01 --min-coverage 10 --minreads2 0 --min-var-freq 0.01 > ${varscan2_indel_vcf}

#post call you need to combine/sort/replace sample column header in vcf
${vcftools}vcf-concat ${varscan2_snp_vcf} ${varscan2_indel_vcf} > ${varscan2_combined_vcf}
cat ${varscan2_combined_vcf}|${vcftools}vcf-sort -c ${varscan2_combined_vcf} > ${varscan2_combined_sort_vcf}
sed -i -e 's/Sample1/'${SampleName}'/g' ${varscan2_combined_sort_vcf}
${bgzip}/bgzip -c ${varscan2_combined_sort_vcf} > ${varscan2_combined_sort_vcf_gz}
${tabix}/tabix -p vcf ${varscan2_combined_sort_vcf_gz}

#vardict call DOCS: https://github.com/AstraZeneca-NGS/VarDictJava
vardict_vcf=${B}.vardict.vcf
vardict_sort_vcf=${B}.vardict_sort.vcf
vardict_sort_vcf_gz=${B}.vardict_sort.vcf.gz

#turn off hex filtering, turn off read2 requirments, quality score 15, #of reads for strand bias
${vardict}vardict -F 0 -r 2 -B 2 -q 15 -G ${db_fa} -f 0.01 -N ${SampleName} -b ${final_bam} -c 1 -S 2 -E 3 -g 4 ${target_bed} | ${vardict}teststrandbias.R | ${vardict}var2vcf_valid.pl -N ${SampleName} -E -f 0.01 > ${vardict_vcf}

${vcftools}./vcf-sort -c ${vardict_vcf} > ${vardict_sort_vcf}
${bgzip}/bgzip -c ${vardict_sort_vcf} > ${vardict_sort_vcf_gz}
${tabix}/tabix -p vcf ${vardict_sort_vcf_gz}

#intersect 3 of 4 callers
consensus_vcf=${B}.consensus.vcf
consensus_vcf_gz=${B}.consensus.vcf.gz

${vcftools}./vcf-isec -n +3 ${freebayes_vcf_gz} ${mutect_vcf_gz} ${varscan2_combined_sort_vcf_gz} ${vardict_sort_vcf_gz} > ${consensus_vcf}
${bgzip}/bgzip -c ${consensus_vcf} > ${consensus_vcf_gz}
${tabix}/tabix -p vcf ${consensus_vcf_gz}

#split multi-allelic sites for annotation
consensus_split_vcf=${B}.consensus_split.vcf
consensus_split_vcf_gz=${B}.consensus_split.vcf.gz

${bcftools}/bcftools norm -N -m -any -O v -o ${consensus_split_vcf} ${consensus_vcf_gz}
${bgzip}/bgzip -c ${consensus_split_vcf} > ${consensus_split_vcf_gz}
${tabix}/tabix -p vcf ${consensus_split_vcf_gz}


#call indels
scalpel_vcf=${B}.scalpel.variants.vcf
scalpel_vcf_gz=${B}.scalpel.variants.vcf.gz

${scalpel}/./scalpel-discovery --single --bam ${final_bam} --bed ${target_bed} --ref ${db_fa} --format vcf --dir ${out_dir} --numprocs 24
mv variants.indel.vcf ${scalpel_vcf}
sed -i -e 's/sample/'${SampleName}'/g' ${scalpel_vcf}
${bgzip}/bgzip -c ${scalpel_vcf} > ${scalpel_vcf_gz}
${tabix}/tabix -p vcf ${scalpel_vcf_gz}

scalpel_filt_vcf=${B}.scalpel.variants.filt.vcf
scalpel_filt_vcf_gz=${B}.scalpel.variants.filt.vcf.gz
#filter indels to match other variant caller params AND Pass
${bcftools}/bcftools filter -i 'FILTER="PASS" && (AD[0]+AD[1] > 10) && (AD[1]/(AD[0]+AD[1]) > 0.01) && (AD[1] >4)' -O v ${scalpel_vcf_gz} -o ${scalpel_filt_vcf}
${bgzip}/bgzip -c ${scalpel_filt_vcf} > ${scalpel_filt_vcf_gz}
${tabix}/tabix -p vcf ${scalpel_filt_vcf_gz}

#take the compliment to generate a vcf with only uniq indels, merge these calls
indel_only=${B}.indel_compliment.vcf
final_vcf=${B}.consensus_and_indel.vcf
final_sort_vcf=${B}_final.vcf
final_sort_vcf_gz=${B}_final.vcf.gz
#merge indels wtih previous calls
${bcftools}/bcftools isec -C -p ${out_dir} -O v ${scalpel_filt_vcf_gz} ${consensus_split_vcf_gz}
mv 0000.vcf ${indel_only}
${vcftools}vcf-concat ${consensus_split_vcf_gz} ${indel_only} > ${final_vcf}
${vcftools}./vcf-sort -c ${final_vcf} > ${final_sort_vcf}
${bgzip}/bgzip -c ${final_sort_vcf} > ${final_sort_vcf_gz}
${tabix}/tabix -p vcf ${final_sort_vcf_gz}

snpeff_anno=${B}.snpeff.vcf

${java7} -jar ${snpeff} -c ${snpeff_conf} -canon hg19 -i vcf -noStats ${final_sort_vcf} > ${snpeff_anno}

onco_anno=${B}_final.anno.TCGAMAF
variants=${B}.onco_parsed.tsv

python ${oncotator}/Oncotator.py -v -i VCF --db-dir ${db_onco}  -o TCGAMAF ${snpeff_anno} ${onco_anno} hg19
python ${onco_parse} ${onco_anno} ${variants}

################################# File Generation and parsing for File Maker Upload #########################################

#generates sample stats ... out puts file with suffix .Stats
stats_out=${B}.Stats
python ${sample_stats} ${B} ${B}.sam.stat ${B}.nodup.stat ${B}_final.stat ${B}.target.depth.sample_summary ${B}.target.depth.sample_interval_summary ${stats_out}

amps=${B}.amps.tsv
python ${amp_calc} ${B}.hotspots.depth.sample_interval_summary ${B}.target.depth.sample_summary ${normal_amps} ${amps} ${B}


#previously generated files contain "#" in header to be awked for combination
variants=${B}.onco_parsed.tsv
stats=${B}.Stats
depths=${B}.target.depth.MeanCovByPosition

variants_upload=${B}.variants.upload
stats_upload=${B}.stats.upload
depths_upload=${B}.depths.upload


amps_upload=${B}.amps.upload

awk '!/#/' ${amps} > ${amps_upload}
#python /project/cpdlab/HiSeqRun/sample_hotspot_norm.py ${B}.hotspots.depth.sample_interval_summary ${B}.target.depth.sample_summary ${B}.normals.sample_avgeraged.interval.avg

awk '!/#/' ${variants} > ${variants_upload}
awk '!/#/' ${stats} > ${stats_upload}
awk '!/#/' ${depths} > ${depths_upload}
