## Utils Used
# Trim Galore - 0.4.1
trim_galore = '/project/cpdlab/Tools/trim_galore/trim_galore'
# Cut Adapt - 1.3
cut_adapt = '/project/cpdlab/Tools/cutadapt/1.3/bin/cutadapt'
# Novalign - 3.02.07
novoalign = '/project/cpdlab/Tools/novocraft/3.02.07/novoalign'
novosort = '/project/cpdlab/Tools/novocraft/3.02.07/novosort'
# AgilentDedup - 1
MBCdedup = '/project/cpdlab/Tools/agilent_dedup/AgilentMBCDedup.jar'
# Samtools - 1.2
samtools = '/project/cpdlab/Tools/samtools/1.2/samtools'
# Picard -  1.46
picard ='/project/cpdlab/Tools/picard/1.46/'
# Bedtools - 2.25
bedtools = '/project/cpdlab/Tools/bedtools/2.25.0/bin/'
# GATK
GATK = '/project/cpdlab/Tools/GATK/3.4-46/GenomeAnalysisTK.jar'
GATK2 = '/home/bigdelia/GATK_3.5/GenomeAnalysisTK.jar'
GATKkey='/project/cpdlab/Tools/GATK/3.4-46/CPD_upenn_gatk.key'
# VCF Tools
vcftools = '/project/cpdlab/Tools/vcftools/0.1.14/src/cpp/vcftools'

#snpEff/SnpSift
snpeff = '/project/cpdlab/Tools/snpEff/4.1l/snpEff.jar'
snpeff_conf =  '/project/cpdlab/Tools/snpEff/4.1l/snpEff.config'
snpsift = '/project/cpdlab/Tools/snpEff/4.1l/SnpSift.jar'
#annovar annotation
annovar_table = '/project/cpdlab/Tools/annovar/17_June_2015/table_annovar.pl'
annovar_humandb ='/project/cpdlab/Tools/humandb/'
#alamut
alamut = '/project/cpdlab/Tools/alamut/1.4.2/./alamut-batch'

## CPD Scripts
FilterAlignScore = '/project/cpdlab/Scripts/PipelineV2.1.1_hpc/FilterAlignScore.py'

## Databases
# hg19
db_nix = '/project/cpdlab/Databases/genomes/hg19_genome.nix'
#db_fa = '/project/cpdlab/Databases/genomes/hg19_genome.fa'
db_fa = '/project/cpdlab/ashkan/muTect/ucsc.hg19.fasta'
db_snp = '/project/cpdlab/Databases/dbsnp/hg19_dbsnp137_correct.vcf'
db_cosmic ='/home/bigdelia/b37_cosmic_v54_120711.vcf'
db_alamut='/project/cpdlab/Databases/Alalmut/hg19.alalmut'

#Mills and 1000G indels
db_indel ='/project/cpdlab/Databases/indels/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf'

## Languages
# Java
java6 = 'java'
java7 = '/project/cpdlab/Tools/jdk1.7.0_80/bin/java'
java8 = '/project/cpdlab/Tools/java8/jdk1.8.0_77/bin/java'
