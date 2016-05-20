# CPD_ETR.py
#
# Ashkan Bigdeli 3/21/2016
#
# A library of methods for extracting, transforming, and reporting next-gen
# sequencing data.

import subprocess, logging, traceback, glob, sys, os
from util import Paths

#global variable of temp directory
tmp = '/project/cpdlab/tmp'
# trim - trims adapters from fastq files using trim_galore
#
# @param1/2 = adapter sequences to be trimmed
# @param3/4 = read sequences *fastq.gz
#
# @return out_fqs = single string with ' ' seperating filenames
def trim(adapt1, adapt2, read1, read2, out_dir):
    try:
        LOG_FILE = out_dir + 'trim_ERROR.log'
        subprocess.call(Paths.trim_galore + ' -q 20 --phred33 --fastqc -a ' + adapt1 + ' -a2 ' + adapt2 +
        ' --stringency 3 -e 0.1 --length 20 --paired ' + read1 + ' ' + read2 +  ' -o ' + out_dir + ' --path_to_cutadapt ' + Paths.cut_adapt, shell=True)

        out_fqs =''
        count = 0
        for filename in glob.glob(out_dir + '*_val_*'):
            out_fqs = out_fqs + filename +' '
            count += 1

        # Return will be inaccurate if trim_galore has already run in destintation folder.
        if ( count > 2):
           raise 'Trim_Galore has already run in this folder!'

    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return out_fqs

# align - aligns fastq's to a reference genome using novocraft novoalign
#         highest alignment score acceptable set to 95 *actually this broke the aligner
#
# @param1/2 = read sequences *.fq.gz
# @param4 = size of fragments for alignment
#
# @return out_sam = string of sam filename
def align(read1, read2, frag_size, method):
    try:
        sample_name = read1.split('.')
        sample_name = sample_name[0]
        LOG_FILE = sample_name + '.align_ERROR.log'
        out_sam = sample_name + '.' + method + 'align.sam'
        subprocess.call(Paths.novoalign + ' -d ' + Paths.db_nix + ' -f ' + read1 + ' ' + read2 + ' -i PE ' + frag_size +
        ' -c 32 -o FullNW -o SAM  > '  + out_sam, shell=True)
    except :
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return out_sam
    
# dedup - removes non-unique reads from sam files by molecular barcode using AgilentMBCdedup
#
# @param1 = aligned sam file
# @param2 = index file associated with read ex: 'sample_name:I2.fastq'
# @param3 = bed file of amplicons
#
# return = string of sorted BAM filename containing only unique reads
def dedup(aligned_sam, index_file, amplicon_bed):
    try:
        LOG_FILE = aligned_sam + '.dedup_ERROR.log'
        dedup_out = aligned_sam.replace('sam', '.bam')
        subprocess.call(Paths.java8 + ' -Xmx72g -jar ' + Paths.MBCdedup + ' -X ' + tmp + ' -b ' +
        amplicon_bed + ' -o ' + dedup_out + ' ' + aligned_sam + ' ' + index_file + ' > dedup_out', shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return dedup_out

# sam2bam - converts sam files to bam files using samtools
#
# @param1 = sam to be convereted
#
# @return = string of bam filename
def sam2bam(sam):
    try:
        LOG_FILE = sam + '.sam2bam_ERROR.log'
        bam_out = sam.replace('sam', 'bam')
        subprocess.call(Paths.samtools + ' view -bS ' + sam + ' > ' + bam_out, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return bam_out

# sort - sorts sam files
#
# @param1 = bam file that requires sorting
def sort(bam):
    try:
        LOG_FILE = bam + '.sort_ERROR.log'
        subprocess.call(Paths.novosort + ' -t ' +tmp + ' -c 32 ' + bam, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit

# fix - uses picard tools to fix read groups using picard-tools
#
# @param1 = bam file to be fixed
# @param2 = bed file of amplicons
# @param3 = index sequence of the read group
# @param4 = name of the sample
def fix(bam, amplicon_bed, read_index, sample_name, lib_name):
    try:
        LOG_FILE = bam + '.fix_read_ERROR.log'
        fix_out = bam.replace ('bam', 'fix.bam')
        subprocess.call(Paths.java6 + ' -Xmx72g -jar ' +  Paths.picard + 'AddOrReplaceReadGroups.jar TMP_DIR=' + tmp + ' I=' + bam + ' O=' +
        fix_out + ' RGID=1 RGLB=' + lib_name  + ' RGPL=Illumina RGPU=' + read_index + ' RGSM=' + sample_name, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return fix_out

        
# index - creats .bai  companion file for bams using samtools
#
# @param1 = bam file that requires index
def index(bam):
    try:
        LOG_FILE = bam + '.index_ERROR.log'
        subprocess.call(Paths.samtools + ' index ' +  bam, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit

# intersect - screens for overlaps in bam files using targets in bed file using bedtools
#
# @param1 = bam file that requires intersection
# @param2 = bed file of amplified regions
#
# return = string of bam filename that has been interesected
def intersect(bam, amplicon_bed):
    try:
        LOG_FILE = bam + '.intersect_ERROR.log'
        intersect_out = bam.replace('bam', 'intersect.bam')
        subprocess.call(Paths.bedtools + 'intersectBed -abam ' + bam + ' -b ' + amplicon_bed + ' > ' + intersect_out, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return intersect_out



# intersect - screens for overlaps in bam files using targets in bed file using bedtools
#
# @param1 = bam file that requires intersection
# @param2 = bed file of amplified regions
#
# return = string of bam filename that has been interesected
def coverage(bam, amplicon_bed):
    try:
        LOG_FILE = bam + '.coverage_ERROR.log'
        coverage_out = bam.replace('bam', 'coverage.calc')
        subprocess.call(Paths.bedtools + 'coverageBed -b ' + bam + ' -a ' + amplicon_bed + ' > ' + coverage_out, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return coverage_out
    
    
    
# filter40 - filters reads with a mapping quality below 40 using samtools
#
# @param1 = bam file to filtered
#
# @return = string of bam filename that has been filtered
def filter40(sam):
    try:
        LOG_FILE = sam + '.filter40_ERROR.log'
        filter40_out = sam.replace('sam', 'q40.sam')
        subprocess.call(Paths.samtools + ' view -h -q 40 ' + sam + ' > ' + filter40_out, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return filter40_out

# filter95 - calls inhouse FilterAlignScore.py to remove alignment qualities below 95
#
# @param1 = sam file to be filtered
#
# @return = string of sam filename that has been flitered
def filter95(sam):
    try:
        LOG_FILE = sam + '.filter95_ERROR.log'
        filter95_out = sam.replace('q40.sam', 'q40.as95.sam')
        #This method Paths needs to be chagned for hpc!
        subprocess.call('python FilterAlignScore.py ' + sam + ' ' + filter95_out , shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return filter95_out

# clip - clips the end of reads using GATK
#
# @param1 = bam file to have reads clipped
# @param2 = name of the sample being clipped
# @param3 = reference database being used
#
# @return = string of bam filename that has been clipped
def clip(bam):
    try:
        LOG_FILE = bam + '.clips_ends_ERROR.log'
        clip_out = bam.replace('bam', 'clip.bam')
        subprocess.call(Paths.java7 + ' -Xmx72g -Djava.io.tmpdir=' + tmp + ' -jar '+ Paths.GATK + ' -T ClipReads -K ' +
        Paths.GATKkey + ' -et NO_ET -I ' + bam + ' -o ' + clip_out + ' -R ' + Paths.db_fa + ' -CR SOFTCLIP_BASES -QT 22', shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return clip_out

# depth - calculates the depth of coverage at 250x
#         generates a .depth file
#
# @param1 = bam file for analysis
# @param2 = sample name
# @param3 = bed file of target regions
def depth(bam, out_dir, sample_name, amplicon_bed):
    try:
        LOG_FILE = bam + '.depth_ERROR.log'
        depth_out = out_dir + sample_name 
        subprocess.call(Paths.java7 + ' -Xmx72g -Djava.io.tmpdir=' + tmp + ' -jar ' + Paths.GATK +
        ' -T DepthOfCoverage -K ' + Paths.GATKkey + ' -et NO_ET -I ' + bam + ' -o ' + depth_out +
        '.depth --minBaseQuality 22 -baseCounts -ct 0 -ct 1 -ct 100 -ct 250 -ct 1000 -L ' + amplicon_bed + ' -R ' + Paths.db_fa, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit

# mpile - generates vcf file
#
# @param1 = bam file to be analyzed
# @param2 = correspond bed file
#
# @return - vcf file
def mpile (bam, amplicon_bed):
    try:
        LOG_FILE = bam + '.mpile_ERROR.log'
        vcf_out = bam.replace('bam', 'piled')
        subprocess.call(Paths.samtools + ' mpileup -f ' + Paths.db_fa + ' ' + bam + ' -l' + amplicon_bed +
        ' >' + vcf_out, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit

 
# haplotyper - calls variants with HaploTypeCaller
#
# @param1 = bam file to have variants called
# @param2 = bed file of amplicons associated with panel
#
# @return = generated vcf file
def haplotyper(bam, amplicon_bed):
    try:
        LOG_FILE = bam + '.discover_ERROR.log'
        vcf_out = bam.replace('bam', 'haplo.vcf')
        subprocess.call(Paths.java7 + ' -Xmx72g -Djava.io.tmpdir=' + tmp + ' -jar ' + Paths.GATK + ' -R ' + Paths.db_fa +
        ' -K ' + Paths.GATKkey + ' -T' + ' HaplotypeCaller -I '+ bam +  ' --dbsnp ' + Paths.db_snp + 
        ' -stand_call_conf 30.0 -stand_emit_conf 10.0 -o ' + vcf_out +
        ' -ERC BP_RESOLUTION --variant_index_type LINEAR --variant_index_parameter 128000 -L ' + amplicon_bed, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return vcf_out

# genotyper - genotypes variants with GATK genotyper
#
# @param1 = vcf file to be genotyped
#
# @return = genotyped vcf file
def genotyper(vcf_in):
    try:
        LOG_FILE = vcf_in + '.genotyper_ERROR.log'
        vcf_out = vcf_in.replace('vcf', 'geno.vcf')
        subprocess.call(Paths.java7 + ' -Xmx72g -Djava.io.tmpdir=' + tmp + ' -jar ' + Paths.GATK + ' -R ' + Paths.db_fa +
        ' -K ' + Paths.GATKkey + ' -T GenotypeGVCFs --max_alternate_alleles 2 -stand_call_conf 30 -stand_emit_conf 10 --variant ' +
        vcf_in + ' -o ' + vcf_out, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return vcf_out

    
# discover - calls variants in GATK discovery mode with Unified Genotyper
#
# @param1 = bam file to have variants called
# @param2 = bed file of amplicons associated with panel
# @param3 = minIndelCount of panel
# @param4 = minIndelFrac of panel
#
# @return = generated vcf file
def uni_discover(bam, amplicon_bed, min_indel_cnt, min_indel_frac ):
    try:
        LOG_FILE = bam + '.discover_ERROR.log'
        vcf_out = bam.replace('bam', 'discover.vcf')
        subprocess.call(Paths.java7 + ' -Xmx72g -Djava.io.tmpdir=' + tmp + ' -jar ' + Paths.GATK + ' -R ' + Paths.db_fa +
        ' -K ' + Paths.GATKkey + ' -et NO_ET -T' + ' UnifiedGenotyper -I '+ bam + ' --dbsnp ' + Paths.db_snp +
        ' -o ' + vcf_out + ' -nct 32 -glm BOTH -mbq 22 -dt NONE -minIndelCnt ' + min_indel_cnt +
        ' -minIndelFrac ' + min_indel_frac + ' -nda -maxAltAlleles 5 -stand_emit_conf 10.0 -L ' + amplicon_bed, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return vcf_out

# alleles
#
# @param1 = bam file to have variants called
# @param2 = bed file of amplicons associated with panel
#
# @return = generated vcf file
def uni_alleles (bam, amplicon_bed, raw_var_vcf, out_dir ):
    try:
        LOG_FILE = bam + '.genotype_alleles_ERROR.log'
        vcf_out = bam.replace('bam', 'gatk_alleles.vcf')
        subprocess.call(Paths.java7 + ' -Xmx72g -Djava.io.tmpdir=' + tmp + ' -jar ' + Paths.GATK + ' -R ' + Paths.db_fa + ' -K ' + Paths.GATKkey +
        ' -et NO_ET -T UnifiedGenotyper -I ' + bam + ' --dbsnp ' + Paths.db_snp + ' -o ' + vcf_out + ' -nct 32 -glm SNP -mbq 22 -dt NONE -alleles:VCF ' +
        raw_var_vcf + ' -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_emit_conf 10.0 -L '+ amplicon_bed, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return vcf_out

# gatk_combine
#
# @param1 = first vcf for combination
# @param2 = second vcf for combination
#
# @return = combined vcf file
def vcf_combine(gatk_discover, gatk_alleles):
    try:
        LOG_FILE = gatk_discover + '.gatk_combine_ERROR.log'
        vcf_out = gatk_discover.replace('vcf', 'gatk_combined.vcf')
        subprocess.call(Paths.java7 +' -Xmx72g -Djava.io.tmpdir=' + tmp + ' -jar ' + Paths.GATK + ' -R ' + Paths.db_fa +
        ' -T CombineVariants -K ' + Paths.GATKkey +' -et NO_ET --variant:variant1 ' + gatk_discover + ' --variant:variant2 ' +
        gatk_alleles + ' -genotypeMergeOptions PRIORITIZE -priority variant2,variant1 -o ' + vcf_out, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return vcf_out

# vcf2table
#
# @param1 = vcf to convert
# @param2 = amplicon bed to intersect
#
# @return = vcf in table format with method specified columns
def vcf2table (vcf_in, amplicon_bed ):
    try:
        LOG_FILE = vcf_in + '.vcf2table_ERROR.log'
        table_out = vcf_in.replace('vcf', 'variant_table')
        subprocess.call(Paths.java7 + ' -Xmx72g -Djava.io.tmpdir=' + tmp + ' -jar ' + Paths.GATK + ' -R ' + Paths.db_fa + ' -K ' + Paths.GATKkey +
        ' -T VariantsToTable -V ' + vcf_in + ' -AMD -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F AC -F AF -GF GT -GF AD -GF DP -GF PL  -o ' +
        table_out + ' -L '+ amplicon_bed, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return table_out

# intervals - generates intervals for indels
#
# @param1 = bam file 
# @param2 = amplicon bed file
# @param3 = sample name - this input is for file naming
#
# @return = a file with intervals of indels
def intervals(bam, amplicon_bed,sample_name, out_dir):
    try:
        LOG_FILE = bam + '.intervals_ERROR.log'
        intervals_out = out_dir + sample_name + '.intervals'
        subprocess.call(Paths.java7 + ' -Xmx72g -Djava.io.tmpdir=' + tmp + ' -jar ' + Paths.GATK + ' -R ' + Paths.db_fa + ' -K ' + Paths.GATKkey +
        ' -T RealignerTargetCreator -I ' + bam + ' -o ' + intervals_out + ' -L '+ amplicon_bed, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return intervals_out

# realigner - realings bam around indels
#
# @param1 = bam file
# @param2 = amplicon file
# @param3 = intervals file generated by GATK RealignerTargetCreator
#
# @return = bam file locally realigned around indels   
def realigner(bam, amplicon_bed, intervals):
    try:
        LOG_FILE = bam + '.realigner_ERROR.log'
        bam_out = bam .replace('bam', 'realigned.bam')
        subprocess.call(Paths.java7 + ' -Xmx72g -Djava.io.tmpdir=' + tmp + ' -jar ' + Paths.GATK + ' -R ' + Paths.db_fa + ' -K ' + Paths.GATKkey +
        ' -T IndelRealigner -I ' + bam + ' -targetIntervals ' + intervals + ' -o ' + bam_out + ' -L '+ amplicon_bed, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return bam_out

# recal - generates table to account to compensate for basecalling errors
#
# @param1 = bam file
# @param2 = ampliconfile
#
# @return = table report of basecalling errors cross ref vs known site vcf's
def recal(bam, amplicon_bed ):
    try:
        LOG_FILE = bam + '.recal_ERROR.log'
        grp_out = bam .replace('bam', 'recal_report.grp')
        subprocess.call(Paths.java7 + ' -Xmx72g -Djava.io.tmpdir=' + tmp + ' -jar ' + Paths.GATK + ' -R ' + Paths.db_fa + ' -K ' + Paths.GATKkey +
        ' -T BaseRecalibrator -knownSites ' + Paths.db_snp + ' -knownSites ' + Paths.db_indel + ' -L '+ amplicon_bed + 
        ' -I ' + bam + ' -o ' + grp_out, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return grp_out

# print_encoded gatk's print function specific to misencoded bases
#
# @param1 = bam file
# @param2 = amplicon bed
# @param1 = gatk recalibration report
# 
# @return = re-encoded bam file
def print_misencoded(bam, amplicon_bed):
    try:
        LOG_FILE = bam + '.print_misencoded_ERROR.log'
        bam_out = bam.replace('bam', 'fix_misencode.bam')
        subprocess.call(Paths.java7 + ' -Xmx72g -Djava.io.tmpdir=' + tmp + ' -jar ' + Paths.GATK + ' -R ' + Paths.db_fa + ' -K ' + Paths.GATKkey +
        ' -nct 32 -T PrintReads -I ' + bam + ' -o ' + bam_out + ' -L '+ amplicon_bed, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return bam_out

# print_recal gatk's print function specific to post base recalibration
#
# @param1 = bam file
# @param2 = amplicon bed
# @param1 = gatk recalibration report
#
# @return = recalibrated bam
def print_recal(bam, amplicon_bed, recal_report):
    try:
        LOG_FILE = bam + '.print_recall_ERROR.log'
        bam_out = bam.replace('bam', 'recal.bam')
        subprocess.call(Paths.java7 + ' -Xmx72g -Djava.io.tmpdir=' + tmp + ' -jar ' + Paths.GATK + ' -R ' + Paths.db_fa + ' -K ' + Paths.GATKkey +
        ' -nct 32 -T PrintReads -BQSR ' + recal_report + ' -I ' + bam + ' -o ' + bam_out + ' -L '+ amplicon_bed, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return bam_out


# filter_vcf gatk's filter function for vcf's. Filter for depth & mapping qual
#
# @param1 = vcf file
#
# @return = filterered vcf
def filter_vcf(vcf_in):
    try:
        LOG_FILE = vcf_in + '.vcf_filter_ERROR.log'
        vcf_out = vcf_in.replace('vcf', 'filtered.vcf')
        subprocess.call(Paths.java7 + ' -Xmx72g -Djava.io.tmpdir=' + tmp + ' -jar ' + Paths.GATK + ' -R ' + Paths.db_fa + ' -K ' + Paths.GATKkey +
        ' -T VariantFiltration -o ' + vcf_out + ' --variant ' + vcf_in + ' --filterExpression "MQ < 40" --filterName "QDandMQ"', shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return vcf_out

# recal_variant recalibrates viriants
#
#  METHOD INCOMPLETE
#
def recal_variant(vcf_in, sample_name, out_dir):
    try:
        LOG_FILE = 'print_recall_error.log'
        tranches = out_dir + sample_name +'.tranches'
        recal = out_dir + sample_name + '.recal '
        apply_recal_input = '-tranchesFile ' + tranches + ' -recalFile ' +recal
        subprocess.call(Paths.java7 + ' -Xmx72g -Djava.io.tmpdir=' + tmp + ' -jar ' + Paths.GATK2 + ' -R ' + Paths.db_fa + ' -K ' + Paths.GATKkey +
        ' -nct 32 -T VariantRecalibrator -I ' + vcf_in + ' --resource:hapmap,VCF,known=false,training=true,prior=15.0 ' + Paths.db_hapmap, +
        ' --resource:dbsnp,VCF,known=false,training=true,prior=12.0' + Paths.dbsnp + ' -an QD -an ReadPosRankSum -an MQRankSum' +
        ' -recalFile ' + recal + ' -tranchesFile ' + tranches + ' -rscriptFile ' + sample_name + 'recal.Plots.R --maxGaussians 4 -mode BOTH',  shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return apply_recal_input

#  mutect2 variant call 
#
#  METHOD INCOMPLETE - CALL WORKS, REF & COSMIC VCF DO NOT MATCH
#   
def mutect2(bam, amplicon_bed):
    try:
        LOG_FILE = bam + '.mutect2_ERROR.log'
        vcf_out = bam.replace('bam', 'mutect.vcf')
        subprocess.call(Paths.java8 + ' -Xmx72g -Djava.io.tmpdir=' + tmp + ' -jar ' + Paths.GATK2 + ' -R ' + Paths.db_fa +
        ' -K ' + Paths.GATKkey + ' -nct 32 -T MuTect2 -I:tumor ' + bam + ' --dbsnp ' + Paths.db_snp + ' --cosmic ' + Paths.db_cosmic +
        ' --artifact_detection_mode -stand_call_conf 30.0 -stand_emit_conf 10.0 -mbq 22 -L ' + amplicon_bed  + ' -o ' + vcf_out, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return vcf_out
    
#  allele_depth 
#
# @param1 = bam file
# @param2 = vcf file
# 
# @return -vcf file that shows allele freq based upon population  
def allele_depth(bam, vcf_in):
    try:
        LOG_FILE = bam + '.allele_depth_ERROR.log'
        vcf_out = vcf_in.replace('vcf', 'allele_depth.vcf')
        subprocess.call(Paths.java7 + ' -Xmx72g -Djava.io.tmpdir=' + tmp + ' -jar ' + Paths.GATK + ' -R ' + Paths.db_fa +
        ' -K ' + Paths.GATKkey + ' -T VariantAnnotator -I ' + bam + ' -A DepthPerAlleleBySample -V '+ vcf_in + 
        ' -o ' + vcf_out, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return vcf_out

#  interval_analysis
#
# @param1 = bam file
# @param2 = vcf file
# 
# @return -vcf file that shows allele freq based upon population  
def uncovered_intervals(bam, depth):
    try:
        LOG_FILE = bam + '.uncovered_intervals_ERROR.log'
        intervals_out = bam.replace('bam', (str(depth) + '_uncovered_intervals'))
        subprocess.call(Paths.java7 + ' -Xmx72g -Djava.io.tmpdir=' + tmp + ' -jar ' + Paths.GATK + ' -R ' + Paths.db_fa +
        ' -K ' + Paths.GATKkey + ' -T FindCoveredIntervals -I ' + bam + ' -u -cov ' + str(depth) + ' -o ' + intervals_out, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return intervals_out
# snpeff
#
# @param1 = vcf to be annotated
#
# @return = annotated vcf file
def snpeff(vcf_in):
    try:
        LOG_FILE = vcf_in + '.snpeff_ERROR.log'
        vcf_out = vcf_in.replace('vcf', 'snpeff.vcf')
        subprocess.call(Paths.java7 + ' -jar ' + Paths.snpeff + ' -c ' + Paths.snpeff_conf + ' hg19 -i vcf -o vcf -noStats ' + vcf_in + 
        ' > ' + vcf_out, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return vcf_out

# snpeff
#
# @param1 = vcf to be annotated
#
# @return = annotated vcf file   
def snpsift(vcf_in):
    try:
        LOG_FILE = vcf_in + '.snpsift_ERROR.log'
        vcf_out = vcf_in.replace('vcf', 'snpsift.vcf')
        subprocess.call(Paths.java7 + ' -jar ' + Paths.snpsift + ' varType ' + vcf_in + ' > ' + vcf_out, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return vcf_out

# snpeff
#
# @param1 = vcf to be annotated
#
# @return = annotated vcf file   
def snpsift_filter(vcf_in):
    try:
        LOG_FILE = vcf_in + '.snpsift_ERROR.log'
        vcf_out = vcf_in.replace('vcf', 'snpsift_filter.vcf')
        subprocess.call( 'cat ' + vcf_in + ' | ' + Paths.java7 + ' -jar ' + Paths.snpsift + ' filter "( intron_variant )" > ' + vcf_out, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return vcf_out


# snpeff
#
# @param1 = vcf to be annotated
#
# @return = annotated vcf file    
def snpsift_extract(vcf_in):
    try:
        LOG_FILE = vcf_in + '.snpsift_ERROR.log'
        vcf_out = vcf_in.replace('vcf', 'extract.vcf')
        subprocess.call(Paths.java7 + ' -jar ' + Paths.snpsift + ' extractFields ' + vcf_in +
        ' CHROM POS ID REF ALT QUAL FILTER AF DP ANN HET HOM VARTYPE > ' + vcf_out, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return vcf_out

# annovar_table
#
# @param1 = vcf file to be converted to tab delimited
#
# @return = tab delimited version of vcf
def annovar_table(vcf_in):
    try:
        LOG_FILE= vcf_in + '.annovar_annotate_ERROR.log'
        table_out = vcf_in.replace('vcf', 'annovar_table')
        subprocess.call(Paths.annovar_table + ' '  + vcf_in + ' ' +Paths.annovar_humandb + ' -buildver hg19 -out ' + table_out + ' -remove -protocol' +
        ' refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all' +
        ' -operation g,r,r,f,f,f,f,f,f,f -nastring . -vcfinput' ,shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return table_out 

       
#alamut - annotates vcf with alamut annotations
#
# @param1 = vcf to be annotated
#
# @return = .vcf 
def alamut(vcf_in):
    try:
        LOG_FILE = vcf_in + '.alamut_ERROR.log'
        vcf_out = vcf_in.replace( 'vcf', 'alamut.vcf')
        subprocess.call(Paths.alamut + ' --in ' + vcf_in + ' --ann ' + vcf_out  + ' --unann ' + vcf_out + '.unann > alamut.out', shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit       

#flagstats - generates statistics for bam files using samtools
#           does not return, but generates a .flagstats file
#
# @param1 = file to be analyzed
def flagstats(sam_bam):
    try:
        LOG_FILE = sam_bam + '.flagstats_ERROR.log'
        stat_out = sam_bam + '.flagstat'
        subprocess.call(Paths.samtools + ' flagstat ' + sam_bam + ' > ' + stat_out, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit

def freebayes(bam):
    try:
        LOG_FILE = bam + '.freebayes_Error.log'
        vcf_out = bam.replace('bam', 'freebayes.vcf')
        subprocess.call(Paths.freebayes + ' -f ' + Paths.db_fa + ' ' + bam + ' > ' + vcf_out, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit

def varscan2_SNP(mpile):
    try:
        LOG_FILE = mpile + '.varscan2_SNP_Error.log'
        vcf_out = mpile.replace('piled',  'varscan2.SNP.vcf')
        subprocess.call(Paths.varscan2 + ' mpileup2snp ' + mpile +  + ' > ' + vcf_out, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
        
def varscan2_INDEL(mpile):
    try:
        LOG_FILE = mpile + '.varscan2_INDEL.log'
        vcf_out = mpile.replace('piled', 'varscan2.INDEL.vcf')
        subprocess.call(Paths.varscan2 + ' mpileup2indel ' + mpile + ' > ' + vcf_out, shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
# vcf_af
#
#@param1 = VCF for allele frequency calc
#
#@return = .freq file with allele frequencies
def vcf_af(vcf_in):
    try:
        LOG_FILE = vcf_in + '.vcf_af_ERROR.log'
        vcf_out = vcf_in.replace('vcf', 'af.vcf')
        subprocess.call(Paths.vcftools + ' --vcf ' + vcf_in + ' --freq', shell=True)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return vcf_out

#this parser is terrible, please dont ever let this be used in a production pipeline.
# it will break @ 3 alt alleles.  
def parse_VAF(vcf_table):
    vcf_out = vcf_table + '.VAF'
    with open(vcf_table) as f1:
        with open(vcf_out,'a') as f2:
            lines = f1.readlines()
            header = lines[:1]
            header = header[0].strip('\n')
            lines = lines[1:]
            f2.write(header + '\tVAF\n')
            for i, line in enumerate(lines):
                alt_write = '0.00'
                line = line.split('\t')
                line[-1] = line[-1].strip('\n')
                if len(line) > 9:
                    vaf_line = line[9].split(',')
                    if len(vaf_line) == 2:
                        ref_allele = float(vaf_line[0])
                        alt_allele = float(vaf_line[1])
                        if alt_allele != 0:
                            alt_freq = (alt_allele/(ref_allele + alt_allele)) * 100
                            alt_write = str(alt_freq)
                    #known bug, hanlding of 3 alts needs to be refined
                    elif len(vaf_line) == 3:
                        ref_allele = float(vaf_line[0])
                        alt_allele = float(vaf_line[1])
                        alt_allele2 = float(vaf_line[2])
                        if alt_allele != 0:
                            alt_freq = (alt_allele/(ref_allele + alt_allele + alt_allele2)) * 100
                            alt_write = str(alt_freq)
                f2.write('\t'.join(map(str,line)) + '\t' + alt_write + '\n')
        

# read_sheet - will have to read excel and return array of sample name & index
def read_sheet():
    print 'i am a placeholder'


# del_file - deletes a file when passed the name
def del_file(file_name):
    try:
        subprocess.call('rm ' + file_name, shell=True)
    except: IOError ('Cannot Remove file' + file_name)

# del_directory - deletes a directory when passed the name
def del_directory(dir_name):
    try:
        subprocess.call('rm -rf ' + dir_name, shell=True)
    except: IOError ('Cannot Remove directory!' + dir_name)

def rename_file (old_name, new_name):
    try:
        subprocess.call('mv ' + old_name + ' ' + new_name, shell=True)
    except: IOError ('Cannot Rename file ! ' + old_name)
    return new_name
# empty_file - checks to see if a file is empty
#
# @param1 = filename/location
#
# return = True if file is empty, False otherwise
def check_empty (filename):
        if (os.stat(filename).st_size== 0):
            raise IOError(filename + ' is empyt!')
