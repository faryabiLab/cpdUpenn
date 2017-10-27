#!/usr/bin/env python
#Ashkan Bigdeli 9/19/2017
#
# archer_api_parse.py
#
# A program to call Archer VM's API and produce CPD FileMaker upload ready files.


__author__ = 'Ashkan Bigdeli'
__version__ = "1.1"
__maintainer__ = "Ashkan Bigdeli"
__email__ = "ashkan.bigdeli@uphs.upenn.edu"
__status__ = "Production"
__modname__ = "archer_api_parse.py"

from argparse import ArgumentParser
import requests
import os

#main api url
url = 'http://170.212.137.168/'
#using admin ashkan bigdeli
auth = ('ashkan.bigdeli@uphs.upenn.edu', 'Sixfingerer20!')


# phone_archer
#
# @param1 = run number as defined by archer VM
#
# Makes API call and pass relevant URL's to appropriate helper method
def phone_archer(job_id, stat_file, variants_file):

    #intial call to get run data
    run_data = requests.get(url + 'rest_api/jobs/' + job_id, auth=auth).json()
    
    #extract relevant urls for each samples
    for i in range (0, len(run_data['samples'])):
        sample_data = requests.get(run_data['samples'][i]['detail_url'], auth=auth).json()
        stats_url = sample_data['results']['read_stats_url']
        sample_name = sample_data['name'].replace('_R1','')
        isoforms_url = sample_data['results']['isoforms']['list_url']
        
        
        fusion_count = parse_fusions(sample_name, isoforms_url, 'fusions', variants_file)
        isoform_count = parse_fusions(sample_name, isoforms_url, 'oncogenic_isoforms', variants_file)
        
        variant_count = fusion_count + isoform_count
        
        # add blank values to variant file for no variants
        if variant_count == 0:
            print_normals(sample_name, variants_file)
        parse_read_stats(sample_name, stats_url, stat_file, variant_count)

 
  
              
# parse_read_stats
#
# @param1 = name of the sample
# @param2 = read stats url, contains all relevant metrics
# @param3 = stats file to store output
#
# Makes API call and parses relevant data
# NOTE: reads are 'normalized' aka downsampled to 1500000 reads
def parse_read_stats(sample_name, stats_url, stats_file, variant_count):
    
    read_stats = requests.get(stats_url, auth=auth).json()
        
    # total reads is post adapter trimming, 0 index of dictionaries list is DNA, 1 is unique
    total_reads_dna = read_stats['read_stat_types'][0]['total_reads']
    mapped_num_dna = read_stats['read_stat_types'][0]['mapped_num']
    mapped_percent_dna = float(read_stats['read_stat_types'][0]['mapped_percent'])
    alignment_loss = 100.00 - mapped_percent_dna
   
    total_unique = read_stats['read_stat_types'][1]['total_reads']
    mapped_percent_unique = read_stats['read_stat_types'][1]['mapped_percent']

    # grab RNA metrics from 'total_stats' dictionary as list element 0 and gsp2 as
    # list element  
    rna_reads = read_stats['total_stats'][0]['rna_reads']
    percent_rna = read_stats['total_stats'][0]['rna_reads_percent']
    gsp2_rna_reads_percent = read_stats['total_stats'][4]['rna_reads_percent']  
    
    sites_per_gsp2_control = read_stats['qc_stats']['avg_unique_rna_start_sites_per_gsp2_control']

    #the below portion of code is an unfortunate necessacity to upload to filemaker
    with open (stats_file, 'a') as stats:
            stats.write(sample_name + '\t' + str(total_reads_dna) + '\t' + str(mapped_num_dna) + '\t' + "{0:.2f}".format(alignment_loss) + '\t' + '\t' + "{0:.2f}".format(mapped_percent_dna) 
                        + '\t' + str(total_unique) + '\t' + "{0:.2f}".format(mapped_percent_unique) + '\t' + str(rna_reads) + '\t' + "{0:.2f}".format(percent_rna)
                        + '\t' + "{0:.2f}".format(gsp2_rna_reads_percent) + '\t' + str(sites_per_gsp2_control) + '\t' + str(variant_count) + '\t\t\t\n') 
    

# parse_fusions
#
# @param1 = sample name
# @param2 = isoform url, contains all archer fusions
# @param3 = file to store variant output
# @param4 = fusion
#
# Makes API call and parses relevant data for fusions depending on type
def parse_fusions(sample_name, isoforms_url, fusion_type, variants_file):
        
    fusions = requests.get(isoforms_url, auth=auth).json()
    
    detection_count = 0

    for x in range (0, len(fusions[fusion_type])):
        if fusions[fusion_type][x]['strong_evidence_aberration'] is True:
           
           detection_count += 1 
           
           unique_starts = fusions['fusions'][x]['unique_start_sites']
           num_reads = fusions[fusion_type][x]['num_reads']
           gsp2_read = fusions[fusion_type][x]['total_gsp2_reads']
           percent_gsp2 = fusions[fusion_type][x]['percent_gsp2_reads']
          
           effect = fusions[fusion_type][x]['novel_type']
           if effect is None:
               effect = 'Fusion'
           #grab the dominant gsp2 for the call, consider its location reference
           gsp2 = fusions[fusion_type][x]['gsp2s'][0]['name']
           ref = gsp2.split('_')[0]
           chrom = gsp2.split('_')[1]
           start_pos = gsp2.split('_')[2]
           
           #extract and create single string for fusion
           gene_arms = fusions[fusion_type][x]['gene_arms']
           gene_arms = ':'.join(gene_arms)
           
           
           break_points = fusions[fusion_type][x]['breakpoint_signature'].split('|')
           
           # Intitiate instance variables for appending
           strand_breaks = ''
           gene_breaks = ''
           genome_num_breaks = ''
           chrom_breaks = ''
           transcript_breaks = ''
           
           #dictionary to house all breakpoints
           # both break points are recorded in case of future development requirments
           bp_dict = dict()
           for i in range(0, len(break_points)):
               signature = break_points[i].split(',')
               chrom1 = signature[0].split(':')[0]
               break1 = signature[0].split(':')[1]
               chrom2 = signature[1].split(':')[0]
               break2 = signature[1].split(':')[1]
               bp_dict[chrom1] = [break1]
               bp_dict[chrom2] = [break2] 
               
           #find relevant break points amongst all anntotations
           for y in range(0, len(fusions[fusion_type][x]['annotations'])):
               transcript_id = fusions[fusion_type][x]['annotations'][y]['transcript_id']
               if fusions[fusion_type][x]['annotations'][y]['start_position'] or  fusions[fusion_type][x]['annotations'][y]['end_position'] in bp_dict.values():
                   gene_breaks = gene_breaks + ':' + fusions[fusion_type][x]['annotations'][y]['gene']
                   genome_num_breaks = genome_num_breaks + ':' + fusions[fusion_type][x]['annotations'][y]['transcripts'][transcript_id]['text']
                   transcript_breaks = transcript_breaks + ':' + transcript_id
                   chrom_breaks = chrom_breaks + ':' + fusions[fusion_type][x]['annotations'][y]['chromosome']
            
           
           # write them to legacy output format   
           with open (variants_file, 'a') as master_var:
                master_var.write(sample_name + '\t' + chrom + '\t' + start_pos + '\t' + ref + '\t' + gene_arms + '\t' + effect + '\t' + '' + '\t' + '' + '\t' +
                             str(num_reads) + '\t' + 'NA' + '\t' + 'NA' + '\t' + ref + '\t' + '' + '\t' + transcript_id + '\t' + '' + '\t' + '' +
                             '\t' + effect + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + '' + '\t' + genome_num_breaks[1:] + '\t' + transcript_breaks[1:] + '\t' + fusion_type + '\t' + ref + '\t' + effect + '\t' +
                             '\t' + gene_breaks[1:] + '\t' + chrom_breaks[1:] + '\t' +  chrom + '\t' + start_pos + '\t' + 'NA' + '\t' + ref + '\t' + gene_arms + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' +
                             'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'Fusion' + '\t' + str(num_reads) + '\t' + str(gsp2_read) + '\t' + str(unique_starts) + '\t' + str(percent_gsp2) + '\t'
                             + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + start_pos + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                            '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\n')

    return detection_count
# print_normals
#
# @param1 = sample name
# @param3 = file to store variant output
#
# If a variant shows no strong fusions or oncogenic isoforms
# blank information should be input for the variant review to acknowledge
# the sample is indeed normal and output in legacy format.
def print_normals(sample_name, variants_file):
           with open (variants_file, 'a') as master_var:
                master_var.write(sample_name + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'No Known Fusion' + '\t' + 'Normal' + '\t' + 'NA' + '\t' + 'NA' + '\t' +
                             'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'No Known Fusion' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'No Known Fusion' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' +
                             'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'Normal' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t'
                             + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                            '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\n')

                                
# make headers
#
# @param1 = File for Run Statistics
# @Param2 = File for Variants
#
# Creates headers necessary for filemaker input
# so helper methods can concat results for legacy required output.
def print_headers(stats_file, variant_file):

    with open (stats_file, 'w+') as stats:
        stats.write('SampleName\tTotalStartingReads\tTotalReadsInputForAlignment\tPercantageLost\tReadsMapped\tPercentageReadMapping\tTotalOnTarget\tPercentageOnTarget\t' +
                         'TotalOnTargetFilterReads\tPercentageOnTargetFilter\tPercentageUsable\tMeanCoverage\tPercentBases_above_0\tPercentBases_above_1\tPercentBases_above_250\t' +
                         'PercentBases_above_1000\tNo._of_amplicons(exonic_part)_below_250x\tClipCount\tNo._of_amplicons(exonic_part)_below_150x\n')
   
    with open (variant_file, 'w+') as master_var:
        master_var.write('Samplename\tChrom\tPos\tRef\tAlt\tVariant Type(SnpEff)\tHomozygous\tQuality\tCoverage\tWarnings\tGene_ID\tGene_name\tBio_type\tTrancript_ID\tExon_ID\tExon_Rank\t' +
                             'Effect\told_AA/new_AA\tOld_codon/New_codon\tCodon_Num(CDS)\tCodon_Degeneracy\tCDS_size\tCodons_around\tAAs_around\tCustom_interval_ID\tc_Change(SnpEff)\t '+ 
                             'p_Change(SnpEff)\tLocation(annovar)\tGene\tConsequence\tTranscript\tcDNA change\tProtein Change\tCHROM(vcf)\tPOS\tdbsnp137_ID\tRef\tAlt\tQUAL\tFILTER\t' + 
                             'INFO\tGenotypeFormat\tGenotype\tAC(1000g)\tAN(1000g)\tAF(1000g)\tAMR_AF\tASN_AF\tEUR_AF\tAFR_AF\tVariantType\tFDP\tFRD\tFAD\tFAF\tgene\tgeneId\ttranscript\t' + 
                             'strand\ttransLen\tprotein\tUniprot\tvarType\tcodingEffect\tvarLocation\tassembly\tgDNAstart\tgDNAend\tgNomen\tcDNAstart\tcDNAend\tcNomen\tpNomen\talt_pNomen\t' +
                             'exon\tintron\tomimId\tdistNearestSS\tnearestSSType\twtSSFScore\twtMaxEntScore\twtNNSScore\twtGSScore\twtHSFScore\tvarSSFScore\tvarMaxEntScore\tvarNNSScore\t' + 
                             'varGSScore\tvarHSFScore\tnearestSSChange\tlocalSpliceEffect\tproteinDomain1\tproteinDomain2\tproteinDomain3\tproteinDomain4\trsId\trsValidated\trsSuspect\t' + 
                             'rsValidations\trsValidationNumber\trsAncestralAllele\trsHeterozygosity\trsClinicalSignificance\trsMAF\trsMAFAllele\trsMAFCount\t1000g_AF\t1000g_AFR_AF\t' + 
                             '1000g_SAS_AF\t1000g_EAS_AF\t1000g_EUR_AF\t1000g_AMR_AF\texacAllFreq\texacAFRFreq\texacAMRFreq\texacEASFreq\texacSASFreq\texacNFEFreq\texacFINFreq\t' + 
                             'exacOTHFreq\texacAFRHmz\texacAMRHmz\texacEASHmz\texacSASHmz\texacNFEHmz\texacFINHmz\texacOTHHmz\texacFilter\texacReadDepth\tespRefEACount\tespRefAACount\t' + 
                             'espRefAllCount\tespAltEACount\tespAltAACount\tespAltAllCount\tespEAMAF\tespAAMAF\tespAllMAF\tespEAAAF\tespAAAAF\tespAllAAF\tespAvgReadDepth\tclinVarIds\t' + 
                             'clinVarOrigins\tclinVarMethods\tclinVarClinSignifs\tclinVarReviewStatus\tclinVarPhenotypes\thgmdId\thgmdPhenotype\thgmdPubMedId\thgmdSubCategory\tcosmicIds\t' + 
                             'cosmicTissues\tcosmicFreqs\tcosmicSampleCounts\tinsNucs\tdelNucs\tsubstType\twtNuc\tvarNuc\tnucChange\tphastCons\tphyloP\twtAA_1\twtAA_3\twtCodon\twtCodonFreq\t' + 
                             'varAA_1\tvarAA_3\tvarCodon\tvarCodonFreq\tposAA\tnOrthos\tconservedOrthos\tconservedDistSpecies\tBLOSUM45\tBLOSUM62\tBLOSUM80\twtAAcomposition\tvarAAcomposition\t' +
                             'wtAApolarity\tvarAApolarity\twtAAvolume\tvarAAvolume\tgranthamDist\tAGVGDclass\tAGVGDgv\tAGVGDgd\tSIFTprediction\tSIFTweight\tSIFTmedian\tMAPPprediction\t' + 
                             'MAPPpValue\tMAPPpValueMedian\tAdjacentSNPflag\n')
                 

    
def main():
    
    #Take in arguements via Job Hook or User Input
    parser = ArgumentParser(description='Summarize Archer Transolocation Results\n' + 
                                'for the Center of Personalized Diagnostics FileMaker Database\n')
    parser.add_argument("-j", "--job_id", required=True, help="Job ID")
    args = parser.parse_args()
    
    #hardcoded output root dir
    archer_dir = '/PathCPD/FromHPC/Archer/'
    job_id = args.job_id

    stats_file = os.path.join(archer_dir, job_id, job_id + '_RunStatsFinal.txt')
    variants_file = os.path.join(archer_dir, job_id, job_id + '_Run_masterVarFinal.txt')
    
    #generate header files
    print_headers(stats_file, variants_file)
    
    #do all the things
    phone_archer(job_id, stats_file, variants_file)
    
main()

if __name__ == "__main__":
    main()
