import sys, os


# This method takes an ArcherDX VM v 5.06 summary file as input and generates
# CPD Filemaker consumable files. 
def archer_variants(run_dir, run_name):
    #create file list of all files to concatenate
    summary_fusions = run_dir + '/summaries/Summary_fusion.csv'
    master_var = run_dir + '/'+ run_name + '_Run_masterVarFinal.txt'
    with open(summary_fusions, 'r') as fusions:
        with open(master_var, 'w+') as master_var:
            # skip first line of summary, print ridiculous headers
            master_var.write("Samplename\tChrom\tPos\tRef\tAlt\tVariant Type(SnpEff)\tHomozygous\tQuality\tCoverage\tWarnings\tGene_ID\tGene_name\tBio_type\tTrancript_ID\tExon_ID\tExon_Rank\tEffect\told_AA/new_AA\tOld_codon/New_codon\tCodon_Num(CDS)\tCodon_Degeneracy\tCDS_size\tCodons_around\tAAs_around\tCustom_interval_ID\tc_Change(SnpEff)\tp_Change(SnpEff)\tLocation(annovar)\tGene\tConsequence\tTranscript\tcDNA change\tProtien Change\tCHROM(vcf)\tPOS\tdbsnp137_ID\tRef\tAlt\tQUAL\tFILTER\tINFO\tGenotypeFormat\tGenotype\tAC(1000g)\tAN(1000g)\tAF(1000g)\tAMR_AF\tASN_AF\tEUR_AF\tAFR_AF\tVariantType\tFDP\tFRD\tFAD\tFAF\tgene\tgeneId\ttranscript\tstrand\ttransLen\tprotein\tUniprot\tvarType\tcodingEffect\tvarLocation\tassembly\tgDNAstart\tgDNAend\tgNomen\tcDNAstart\tcDNAend\tcNomen\tpNomen\talt_pNomen\texon\tintron\tomimId\tdistNearestSS\tnearestSSType\twtSSFScore\twtMaxEntScore\twtNNSScore\twtGSScore\twtHSFScore\tvarSSFScore\tvarMaxEntScore\tvarNNSScore\tvarGSScore\tvarHSFScore\tnearestSSChange\tlocalSpliceEffect\tproteinDomain1\tproteinDomain2\tproteinDomain3\tproteinDomain4\trsId\trsValidated\trsSuspect\trsValidations\trsValidationNumber\trsAncestralAllele\trsHeterozygosity\trsClinicalSignificance\trsMAF\trsMAFAllele\trsMAFCount\t1000g_AF\t1000g_AFR_AF\t1000g_SAS_AF\t1000g_EAS_AF\t1000g_EUR_AF\t1000g_AMR_AF\texacAllFreq\texacAFRFreq\texacAMRFreq\texacEASFreq\texacSASFreq\texacNFEFreq\texacFINFreq\texacOTHFreq\texacAFRHmz\texacAMRHmz\texacEASHmz\texacSASHmz\texacNFEHmz\texacFINHmz\texacOTHHmz\texacFilter\texacReadDepth\tespRefEACount\tespRefAACount\tespRefAllCount\tespAltEACount\tespAltAACount\tespAltAllCount\tespEAMAF\tespAAMAF\tespAllMAF\tespEAAAF\tespAAAAF\tespAllAAF\tespAvgReadDepth\tclinVarIds\tclinVarOrigins\tclinVarMethods\tclinVarClinSignifs\tclinVarReviewStatus\tclinVarPhenotypes\thgmdId\thgmdPhenotype\thgmdPubMedId\thgmdSubCategory\tcosmicIds\tcosmicTissues\tcosmicFreqs\tcosmicSampleCounts\tinsNucs\tdelNucs\tsubstType\twtNuc\tvarNuc\tnucChange\tphastCons\tphyloP\twtAA_1\twtAA_3\twtCodon\twtCodonFreq\tvarAA_1\tvarAA_3\tvarCodon\tvarCodonFreq\tposAA\tnOrthos\tconservedOrthos\tconservedDistSpecies\tBLOSUM45\tBLOSUM62\tBLOSUM80\twtAAcomposition\tvarAAcomposition\twtAApolarity\tvarAApolarity\twtAAvolume\tvarAAvolume\tgranthamDist\tAGVGDclass\tAGVGDgv\tAGVGDgd\tSIFTprediction\tSIFTweight\tSIFTmedian\tMAPPprediction\tMAPPpValue\tMAPPpValueMedian\tAdjacentSNPflag"+"\n")
            line = fusions.readline()
            for line in fusions:
                #strip any new line char + split on ,
                line = line.strip()
                line = line.split(',')


                #only extract results which have known fusions. 
                if line[2] == 'FALSE':
                    continue
                
                #extract 5' and 3' fusions
                five_prime = line[12]
                five_prime = five_prime.split('|')
                five_prime_transcript = five_prime[0]

                five_prime_exon = five_prime[1]
                five_prime_exon = five_prime_exon.replace("exon:","")
                
                
                
                three_prime = line[13]
                three_prime = three_prime.split('|')
                three_prime_transcript = three_prime[0]

                three_prime_exon = three_prime[1]
                three_prime_exon = three_prime_exon.replace("exon:","")
                
                
                location ='Fusion'

                sample_name = line[0]
                sample_name = sample_name.replace("_R1", "")
                
                gsp2_gene = line[15]
                gsp2_gene = gsp2_gene.split('_')
                gene = gsp2_gene[0]
                
                chrom_line = line[16]
                chrom_line = chrom_line.split(':')
                chrom = chrom_line[0]
                start_pos = chrom_line[1]
                
                
                fusion = line[1] 
                # total RNA Reads   
                depth = line[18]
                supportive_reads = line[21]
                unique_starts = line [22]
                percent_rna = line[23]
                




                master_var.write(sample_name + '\t' + chrom + '\t' + start_pos + '\t' + gene + '\t' + fusion + '\t' + "NA" + '\t' + 'NA' + '\t' + "" + '\t' +
                             depth + '\t' + 'NA' + '\t' + "" + '\t' + gene + '\t' + 'NA' + '\t' + "" + '\t' + 'NA' + '\t' + "" +
                             '\t' + fusion + '\t' + 'NA' + '\t' + "" + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + "" + '\t' + three_prime_transcript + '\t' + three_prime_exon + '\t' + location + '\t' + gene + '\t' + "Fusion" + '\t' + "" +
                             '\t' + five_prime_transcript + '\t' + five_prime_exon + '\t' + chrom + '\t' + start_pos + '\t' + "" + '\t' + "" + '\t' + "" + '\t' + "" +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + "" + '\t' + "" + '\t' + "" + '\t' + "" + '\t' +
                             "" + '\t' + "" + '\t' + "" + '\t' + "" + '\t' + depth + '\t' + supportive_reads + '\t' + unique_starts + '\t' + percent_rna + '\t'
                             + "" + '\t' + "" + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + "" + '\t' + "" + '\t' + "" + '\t' + "" + '\t' + "" + '\t' + "" +
                             '\t' + "" + '\t' + 'NA' + '\t' + "" + '\t' + "" + '\t' + "" + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + "" + '\t' + "" + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + "" +
                            '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\n')
                

# This method takes an ArcherDX VM v 5.06 summary file as input and generates
# CPD Filemaker consumable files. 
def sum_stats(run_dir, run_name):
    #create file list of all files to concatenate
    summary_stats = run_dir + '/summaries/Summary_sample.csv'
    metrics = run_dir + '/'+ run_name + '_RunStatsFinal.txt'
    with open(summary_stats, 'r') as stats:
        with open (metrics, 'w+') as metrics:
            metrics.write('SampleName\tTotalStartingReads\tTotalReadsInputForAlignment\tPercantageLost\tReadsMapped\tPercentageReadMapping\tTotalOnTarget\tPercentageOnTarget\t' +
                         'TotalOnTargetFilterReads\tPercentageOnTargetFilter\tPercentageUsable\tMeanCoverage\tPercentBases_above_0\tPercentBases_above_1\tPercentBases_above_250\t' +
                         'PercentBases_above_1000\tNo._of_amplicons(exonic_part)_below_250x\tClipCount\tNo._of_amplicons(exonic_part)_below_150x\n')
            #skip header
            line = stats.readline()
            for line in stats:
                   
                line = line.strip()
                line = line.split(',')
                
                sample_name = line[0]
                sample_name = sample_name.replace("_R1", "")
                
                # dna + rna metrics
                start_reads = line[2]
                aligned_reads = line[11]
                percent_aligned = float(aligned_reads)/float(start_reads) * 100
                percent_lost = 100.00 - percent_aligned

                #rna_reads = reads on target
                rna_reads = line[3]
                # percent rna = % on target
                percent_rna = float(line[5])
                #unique rna = reads after quality filtering 
                unique_rna = line[4]
                # perecent_unique_rna
                percent_unique_rna = float(line[8])
                #percent_unique_target = % usable
                percent_unique_target = float(line[16])
                
                #mean_rna_frag
                mean_rna_frag = float(line[18])
                
                metrics.write(sample_name + '\t' + start_reads + '\t' + aligned_reads + '\t' + "{0:.2f}".format(percent_lost) + '\t' + rna_reads 
                             + '\t' + "{0:.2f}".format(percent_aligned) + '\t' + rna_reads + '\t' + "{0:.2f}".format(percent_rna) + '\t' + unique_rna + '\t' + "{0:.2f}".format(percent_unique_rna)
                             + '\t' + "{0:.2f}".format(percent_unique_target) + '\t' + "{0:.2f}".format(mean_rna_frag) + '\t\t\t\t\n') 

                



def main():
    run_dir = sys.argv[1]
    run_name = sys.argv[2]
    archer_variants(run_dir, run_name)
    sum_stats(run_dir, run_name)

main()
