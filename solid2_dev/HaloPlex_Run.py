# HaloPlex_Run.py

# Ashkan Bigdeli 3/30/2016

# The pipeline from which all Solid Version 2 samples will run...one day.

from methods import CPD_ETR
from panels import Solid2
import logging, traceback, sys, glob, os

def sample_align (run):
    LOG_FILE = run.out_dir + run.sample_name + '_align_error.log'
    try:
        LOG_FILE = run.out_dir + run.sample_name + '_run_error.log'
        trimmed_fqs = CPD_ETR.trim(run.adapter1, run.adapter2, run.read1, run.read2, run.out_dir)
        trimmed_files = trimmed_fqs.split(' ')
        trimmed_fqs = CPD_ETR.trim(run.adapter1, run.adapter2, run.read1, run.read2, run.out_dir)
        trim_fq1 = trimmed_files[0]
        trim_fq2 = trimmed_files[1]
        aligned_sam = CPD_ETR.align(trim_fq1, trim_fq2, run.frag_size)
    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit
    return aligned_sam
    
def sample_run(aligned_sam, run, method):
    LOG_FILE = run.out_dir + run.sample_name + '_run_error.log'
    try:        
        if (method == 'dedup'):
            bam = CPD_ETR.dedup(aligned_sam, run.index2, run.amplicon_bed)
            bam = CPD_ETR.fix(bam, run.amplicon_bed, run.index2, run.sample_name, run.lib_name)
            bam = CPD_ETR.index(bam)           
            CPD_ETR.flagstats(bam)
        else:
            bam = CPD_ETR.sam2bam(aligned_sam)
            bam = CPD_ETR.sort(bam)
            bam = CPD_ETR.fix(bam, run.amplicon_bed, run.index2, run.sample_name, run.lib_name)
            bam = CPD_ETR.index(bam)        
            CPD_ETR.flagstats(bam)

        #intersect
        intersect_bam = CPD_ETR.intersect(bam, run.target_bed)
        intersect_bam = CPD_ETR.fix(intersect_bam, run.target_bed, run.index2, run.sample_name, run.lib_name)
        CPD_ETR.sort(intersect_bam)
        CPD_ETR.index(intersect_bam)
        
        #begin indel realignment
        recoded_bam = CPD_ETR.print_misencoded(intersect_bam, run.target_bed)
        CPD_ETR.sort(recoded_bam)        
        CPD_ETR.index(recoded_bam)

        intervals = CPD_ETR.intervals(recoded_bam, run.target_bed, run.sample_name, run.out_dir)
        realigned_bam = CPD_ETR.realigner(recoded_bam, run.target_bed, intervals)
        CPD_ETR.sort(realigned_bam)
        CPD_ETR.index(realigned_bam)

        recal_report = CPD_ETR.recal(realigned_bam, run.target_bed)
        recal_bam = CPD_ETR.print_recal(realigned_bam, run.target_bed, recal_report)
        CPD_ETR.sort(realigned_bam)
        CPD_ETR.index(recal_bam)
         
        #soft clip bam
        clip_bam = CPD_ETR.clip(recal_bam) 
        
        #rename generate stats of final bam
        final_bam = CPD_ETR.rename_file(clip_bam, (run.out_dir + run.sample_name + '.' + method + '.final.bam'))        
        CPD_ETR.sort(final_bam)                              
        CPD_ETR.index(final_bam)
        CPD_ETR.flagstats(final_bam)
        CPD_ETR.depth(final_bam, run.out_dir, run.sample_name, run.target_bed) 
        CPD_ETR.coverage(final_bam, run.amplicon_bed)        
        CPD_ETR.coverage(final_bam, run.covered_bed)        

        #call variants mutect2 and annotate
        mutect_vcf = CPD_ETR.mutect2(final_bam, run.target_bed)        
        mfiltered_vcf =  CPD_ETR.filter_vcf(mutect_vcf)
        snp_vcf = CPD_ETR.snpeff(mfiltered_vcf)
        snp_sift_vcf = CPD_ETR.snpsift(snp_vcf)
        annovar_vcf =CPD_ETR.annovar_table(mfiltered_vcf)
        alamut_vcf = CPD_ETR.alamut(mfiltered_vcf)       
        
        #call variants with freebayes and annotate
        freebayes_vcf = CPD_ETR.freebayes(final_bam)
        
        #call variants with varscan
        mpile = CPD_ETR.mpile(final_bam, run.amplicon_bed)        
        varscan2_snp = CPD_ETR.varscan2_SNP(mpile)
        varscan2_INDEL = CPD_ETR.varscan2_INDEL(mpile)

        
        files = glob.glob(run.out_dir)
        for f in files:
            if ("final","flagstat","depth", "profile", "log", "txt") not in f:
                CPD_ETR.del_file(f)
                

    except:
        logging.basicConfig(filename=LOG_FILE)
        logging.critical(traceback.format_exc())
        sys.exit

def main():
    sample_name = sys.argv[1]
    read1 = sys.argv[2]
    read2 = sys.argv[3]
    read_index = sys.argv[4]
    index2 = sys.argv[5]
    out_dir = sys.argv[6]
    run = Solid2.Solid2(sample_name, read1, read2, read_index, index2, out_dir)
    aligned_sam = sample_align(run)
    sample_run(aligned_sam, run, 'dedup')
    sample_run(aligned_sam, run, 'allseq')

main()
