# HaloPlex_Run.py

# Ashkan Bigdeli 3/30/2016

# The pipeline from which all Solid Version 2 samples will run...one day.

from methods import CPD_ETR
from panels import Solid2
import logging, traceback, sys, glob


def sample_run(sample_name, read1, read2, read_index, index2, out_dir):
    try:
        #create solid2 object
        run = Solid2.Solid2(sample_name, read1, read2, read_index, index2, out_dir)
        LOG_FILE = run.out_dir + run.sample_name + '_run_error.log'
        
        #trim -> align -> deduplicate reads
        #trimmed_fqs = CPD_ETR.trim(run.adapter1, run.adapter2, run.read1, run.read2, run.out_dir)
        #trimmed_files = trimmed_fqs.split(' ')
        #trimmed_fqs = CPD_ETR.trim(run.adapter1, run.adapter2, run.read1, run.read2, run.out_dir)
        #trim_fq1 = trimmed_files[0]
        #trim_fq2 = trimmed_files[1]               
        #aligned_sam = CPD_ETR.align(trim_fq1, trim_fq2, run.frag_size)
        #bam = CPD_ETR.dedup(aligned_sam, run.index2, run.amplicon_bed)
        #bam = CPD_ETR.fix(bam, run.amplicon_bed, run.index2, run.sample_name, run.lib_name)
        #bam = '/project/cpdlab/cpdUpenn/solid2_dev/HiSeqSamples/CPDV150986-35_ucsc19/CPDV150986-35.align.consensus.fix.reordered.bam'
        #CPD_ETR.index(bam)
        #CPD_ETR.flagstats(bam)

        #begin indel realignment
        # recoded_bam = CPD_ETR.print_misencoded(bam, run.amplicon_bed)
        #recoded_bam = '/project/cpdlab/cpdUpenn/solid2_dev/HiSeqSamples/CPDV150986-35_ucsc19/CPDV150986-35.align.consensus.fix.reordered.fix_misencode.NCTTEST.bam'
        #CPD_ETR.index(recoded_bam)
        #intervals = CPD_ETR.intervals(recoded_bam, run.amplicon_bed, run.sample_name, run.out_dir)
        #realigned_bam = '/project/cpdlab/cpdUpenn/solid2_dev/HiSeqSamples/CPDV150986-35_ucsc19/CPDV150986-35.align.consensus.fix.reordered.fix_misencode.realigned.bam' #CPD_ETR.realigner(recoded_bam, run.amplicon_bed, intervals)
        #CPD_ETR.index(realigned_bam)
        #recal_report = CPD_ETR.recal(realigned_bam, run.amplicon_bed)
        #recal_bam = CPD_ETR.print_recal(realigned_bam, run.amplicon_bed, recal_report)
        #CPD_ETR.index(recal_bam)
         
        #soft clip bam
        #clip_bam = CPD_ETR.clip(recal_bam) 
        #CPD_ETR.index(clip_bam)
        #CPD_ETR.flagstats(clip_bam)

        #intersect bam with bed
        #intersect_bam = CPD_ETR.intersect(clip_bam, run.amplicon_bed)
        #intersect_bam = CPD_ETR.fix(intersect_bam, run.amplicon_bed, run.index2, run.sample_name, run.lib_name)

        #rename 
        #final_bam = CPD_ETR.rename_file(intersect_bam, (run.out_dir + run.sample_name + '.final.bam'))
        #CPD_ETR.sort(final_bam)                              
        #CPD_ETR.index(final_bam)
        #CPD_ETR.flagstats(final_bam)
        #CPD_ETR.depth(final_bam, run.out_dir, run.sample_name, run.amplicon_bed)        
        #final_bam = '/project/cpdlab/cpdUpenn/solid2_dev/HiSeqSamples/CPDV150986-35_ucsc19/CPDV150986-35.final.bam'
        #mutect_vcf = CPD_ETR.mutect2(final_bam, run.amplicon_bed)
        # call variants - mapping quality set to 40
        #snp_indels_vcf = CPD_ETR.haplotyper(final_bam, run.amplicon_bed)
        filtered_vcf =  '/project/cpdlab/cpdUpenn/solid2_dev/HiSeqSamples/CPDV150986-35_ucsc19/CPDV150986-35.final.mutect.filtered.vcf'#CPD_ETR.filter_vcf(mutect_vcf)
        #recal_vcf = CPD_ETR.recal_variant(filtered_vcf)
        #genotyped_vcf = CPD_ETR.genotyper(filtered_vcf)

        #perform annotations
        snp_vcf = CPD_ETR.snpeff(filtered_vcf)
        snp_sift_vcf = CPD_ETR.snpsift(filtered_vcf)
        snp_sift_extract = CPD_ETR.snpsift_extract(filtered_vcf)
        #annovar_vcf =CPD_ETR.annovar_table(filtered_vcf)
        #alamut_vcf = CPD_ETR.alamut(filtered_vcf)

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
    sample_run(sample_name, read1, read2, read_index, index2, out_dir)

main()
