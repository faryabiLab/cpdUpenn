import sys 
import traceback


# onco parse
def onco_parse(anno,anno_out,sample):
    with open(anno, 'r') as anno_in:
        with open(anno_out, 'w+') as parsed:
            parsed.write("#Samplename\tChrom\tPos\tRef\tAlt\tVariant Type(SnpEff)\tHomozygous\tQuality\tCoverage\tWarnings\tGene_ID\tGene_name\tBio_type\tTrancript_ID\tExon_ID\tExon_Rank\tEffect\told_AA/new_AA\tOld_codon/New_codon\tCodon_Num(CDS)\tCodon_Degeneracy\tCDS_size\tCodons_around\tAAs_around\tCustom_interval_ID\tc_Change(SnpEff)\tp_Change(SnpEff)\tLocation(annovar)\tGene\tConsequence\tTranscript\tcDNA change\tProtien Change\tCHROM(vcf)\tPOS\tdbsnp137_ID\tRef\tAlt\tQUAL\tFILTER\tINFO\tGenotypeFormat\tGenotype\tAC(1000g)\tAN(1000g)\tAF(1000g)\tAMR_AF\tASN_AF\tEUR_AF\tAFR_AF\tVariantType\tFDP\tFRD\tFAD\tFAF\tgene\tgeneId\ttranscript\tstrand\ttransLen\tprotein\tUniprot\tvarType\tcodingEffect\tvarLocation\tassembly\tgDNAstart\tgDNAend\tgNomen\tcDNAstart\tcDNAend\tcNomen\tpNomen\talt_pNomen\texon\tintron\tomimId\tdistNearestSS\tnearestSSType\twtSSFScore\twtMaxEntScore\twtNNSScore\twtGSScore\twtHSFScore\tvarSSFScore\tvarMaxEntScore\tvarNNSScore\tvarGSScore\tvarHSFScore\tnearestSSChange\tlocalSpliceEffect\tproteinDomain1\tproteinDomain2\tproteinDomain3\tproteinDomain4\trsId\trsValidated\trsSuspect\trsValidations\trsValidationNumber\trsAncestralAllele\trsHeterozygosity\trsClinicalSignificance\trsMAF\trsMAFAllele\trsMAFCount\t1000g_AF\t1000g_AFR_AF\t1000g_SAS_AF\t1000g_EAS_AF\t1000g_EUR_AF\t1000g_AMR_AF\texacAllFreq\texacAFRFreq\texacAMRFreq\texacEASFreq\texacSASFreq\texacNFEFreq\texacFINFreq\texacOTHFreq\texacAFRHmz\texacAMRHmz\texacEASHmz\texacSASHmz\texacNFEHmz\texacFINHmz\texacOTHHmz\texacFilter\texacReadDepth\tespRefEACount\tespRefAACount\tespRefAllCount\tespAltEACount\tespAltAACount\tespAltAllCount\tespEAMAF\tespAAMAF\tespAllMAF\tespEAAAF\tespAAAAF\tespAllAAF\tespAvgReadDepth\tclinVarIds\tclinVarOrigins\tclinVarMethods\tclinVarClinSignifs\tclinVarReviewStatus\tclinVarPhenotypes\thgmdId\thgmdPhenotype\thgmdPubMedId\thgmdSubCategory\tcosmicIds\tcosmicTissues\tcosmicFreqs\tcosmicSampleCounts\tinsNucs\tdelNucs\tsubstType\twtNuc\tvarNuc\tnucChange\tphastCons\tphyloP\twtAA_1\twtAA_3\twtCodon\twtCodonFreq\tvarAA_1\tvarAA_3\tvarCodon\tvarCodonFreq\tposAA\tnOrthos\tconservedOrthos\tconservedDistSpecies\tBLOSUM45\tBLOSUM62\tBLOSUM80\twtAAcomposition\tvarAAcomposition\twtAApolarity\tvarAApolarity\twtAAvolume\tvarAAvolume\tgranthamDist\tAGVGDclass\tAGVGDgv\tAGVGDgd\tSIFTprediction\tSIFTweight\tSIFTmedian\tMAPPprediction\tMAPPpValue\tMAPPpValueMedian\tAdjacentSNPflag"+"\n")
            count = 1
            for line in anno_in:
                line = line.strip()
                count = count + 1
                #Skip Oncotator header lines
                if line.startswith('##'):
                    continue;
                if line.startswith('Hugo_Symbol'):
                    continue;
                if line.startswith('#'):
                    continue;
                #split tabs into array
                values = line.split('\t')

                #begin assigning values
                chrom = 'chr' + values[4] #Chromosome
                start_pos = values[5] # Start_position
                end_pos = values[6] # End_position
                ref = values[10] # Reference_Allele
                alt = values[12] # Tumor_Seq_Allele2
                gene = values[0] # Hugo_Symbol
                geneid = values[1] # Entrez_Gene_ID
                onc_effect = values[9] # Variant_Type
                transcript_loc = values[37] # Transcript_Exon

		# Added 5/24/2018 - now reporting strand orientation (forward/reverse) in Solid
		transcript_strand = values[36] # Transcript_Strand
		if transcript_strand == "+":
		    strand="1"
		elif transcript_strand == "-":
		    strand="-1"

                onc_gchange = values[34] # Genome_Change
                onc_type = values[9] # Variant_Type
                sample_name = sample # change to ${prefix} input
                rs_id = values[378] # dbsnp_147_RS << Updated
                rs_val = values[390] # dbsnp_147_VLD

                onc_c_change = values[39] # cDNA_Change
                onc_p_change = values[41] # Protein_Change
                onc_codon_change = values[40] # Codon_Change
                c_dna_start = ""
                c_dna_stop = ""
                #if there is a codon change, extract the cD'NA' start and stop
                if (len(onc_codon_change) > 1):
                    try:
                        start_stop = onc_codon_change[onc_codon_change.index("(") + 1:onc_codon_change.rindex(")")].split('-')
                        c_dna_start = start_stop[0]
                        c_dna_stop = start_stop[1]
                    except:
                        c_dna_start =""
                        c_dna_stop = ""
                    
                onc_gchange = values[34] # Genome_Change

                hgvs_c_change = values[252] # HGVS_coding_DNA_change
                hgvs_p_change = values[254] # HGVS_protein_change
                hgvs_g_change = values[253] # HGVS_genomic_change

                sift_pred = values[315] # dbNSFP_SIFT_pred
                clinvar_id = values[124] # ClinVar_rs
                hgnc_omimid = values[239] # HGNC_OMIM ID(supplied by NCBI)

                onekg_ac = values[80] # 1000gp3_AC
                onekg_an = values[84] # 1000gp3_AN
                onekg_af = values[81] # 1000gp3_AF
                onekg_amr = values[83] # 1000gp3_AMR_AF
                onekg_eur = values[91] # 1000gp3_EUR_AF
                onekg_afr = values[82] # 1000gp3_AFR_AF
                onekg_asn = values[89] # 1000gp3_EAS_AF
                
                #parse snpeff annotations from onco
                #snpeff annotations are separated by a | in either annotation sets of
                # 15 or 30. effect transcript, p and c changes are consistently located in between pipes
                snpeff_ann = values[104] # ANN
                snpeff_split = snpeff_ann.split('|')
                snpeff_effect= snpeff_split[1]
                snpeff_transcript = snpeff_split[6]
                snpeff_gene = snpeff_split[3]
                snpeff_coding = snpeff_split[7]
                if (len(snpeff_split[8]) > 0):
                    snpeff_rank = snpeff_split[8]
                    snpeff_rank = snpeff_rank.split('/')
                    snpeff_rank = snpeff_rank[0]
                else:
                    snpeff_rank = ""
                snpeff_c_change = snpeff_split[9]
                snpeff_p_change = snpeff_split[10]
                
            
                variant_class = values[8] # Variant_Classification
                

                if variant_class == "Intron":
                    location = "Intronic"
                elif variant_class == "3'UTR":
                    location = "3'UTR"
                elif variant_class == "5'UTR":
                    location = "5'UTR"
                elif variant_class == "5'FLANK":
                    location = "5'FLANK"
                elif variant_class == "3'FLANK":
                    location = "3'FLANK"
                elif variant_class == "IGR":
                    location = "Intergenic"
                elif variant_class == "Splice_Site":
                    location = "Splice Site"
                else:
                    location = "Exonic"
                     
                #additional conditionals to ensure appropriate location
                if ((location is "Intron") and (snpeff_coding is "Coding")):
                    location = "Exonic"
                if snpeff_p_change.startswith("n."):
                    location = "Intron"
                #provide intron if necessary
                intron = ""
                if "intron" in snpeff_effect:
                    intron = transcript_loc

                depth = values[395] # depth_across_samples
                DP4 = values[125].split(',')
                alt_depth = int(DP4[2]) + int(DP4[3]) # DP4 column [2]+[3]
                ref_depth = int(DP4[0]) + int(DP4[1]) # DP4 column [0]+[1]

                vaf = float(values[261])*100 # allele_frequency
                #if temp_vaf == 50 or temp_vaf == 100:
                #    vaf = values[79]/float(values[79]+values[80])
                #else:
                #    vaf = temp_vaf

                qual = values[407] # qual

                #---------------------------------------------------------------

                #perform checks to determine which vcf generated the annotation. Unnecessary for SWIFT!
                #checking for freebayes

                #if (len(values[350]) >= 1): # SRF
                 #   if (len(values[340]) < 1): # SAF
                  #      alt_r1 = 0
                   # else:
                    #    alt_r1 = int(values[340]) # SAF
                   # alt_r2 = int(values[343]) # SAR
                   # ref_r1 = int (values[350]) # SRF
                   # ref_r2 = int(values[352]) # SRR

                   # alt_depth = alt_r1 + alt_r2
                   # ref_depth = ref_r1 + ref_r2
                   # vaf = vaf_percent(alt_depth, int(depth))

                   # alt_qual = values[321] # QA
                   # ref_qual = values[323] # QR
                   # qual = avg_qual(ref_qual,alt_qual, ref_depth, alt_depth) 

                   # if float(qual) < 20.0:
                   #     print str(start_pos)
                   #     continue

                #check for scalpel
                #if (len(values[362]) > 1): # allelic_depth
                #    allelic_depth = values[362].split(',')
                #    alt_depth = allelic_depth[1]
                #    ref_depth = allelic_depth[0]
                #    depth = int(alt_depth) + int(ref_depth)
                #    vaf = vaf_percent(alt_depth, depth)
                #    qual = "> 25"
                #check for pindel
                #if (len(values[358]) > 0): # VD__FORMAT__
                #    alt_depth = values[358] # VD__FORMAT__
                #    ref_depth = values[329] # RD
                #    depth = int(alt_depth) + int(ref_depth)
                #    vaf = vaf_percent(alt_depth, depth)
                #    qual = "NA"

		# Quality checks
                if int(depth) < 50 :
                    print str(start_pos + " << Insufficient depth")
                    continue
                if alt_depth < 5 :
                    print str(start_pos+ " << Insufficient alt_depth")
                    continue
                if float(vaf) < 0.1:
                    # replace with equivalent from GATK VCF... somehow...
                    print str(start_pos + " << Variant frequency too low to pass threshold")
                    continue
                parsed.write(sample_name + '\t' + chrom + '\t' + start_pos + '\t' + ref + '\t' + alt + '\t' + onc_type + '\t' + 'NA' + '\t' + qual + '\t' +
                             str(depth) + '\t' + 'NA' + '\t' + geneid + '\t' + snpeff_gene + '\t' + 'NA' + '\t' + snpeff_transcript + '\t' + 'NA' + '\t' + snpeff_rank +
                             '\t' + snpeff_effect + '\t' + 'NA' + '\t' + onc_codon_change + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + rs_id + '\t' + snpeff_c_change + '\t' + snpeff_p_change + '\t' + location + '\t' + gene + '\t' + snpeff_effect + '\t' + snpeff_transcript +
                             '\t' + onc_c_change + '\t' + onc_p_change + '\t' + chrom + '\t' + start_pos + '\t' + rs_id + '\t' + ref + '\t' + alt + '\t' + qual +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + onekg_ac + '\t' + onekg_an + '\t' + onekg_af + '\t' + onekg_amr + '\t' +
                             onekg_asn + '\t' + onekg_eur + '\t' + onekg_afr + '\t' + onc_type + '\t' + str(depth) + '\t' + str(ref_depth) + '\t' + str(alt_depth) + '\t' + str(vaf) + '\t'
                             + snpeff_gene + '\t' + geneid + '\t' + 'NA' + '\t' + strand + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + start_pos + '\t' + end_pos + '\t' + hgvs_g_change + '\t' + c_dna_start + '\t' + c_dna_stop + '\t' + hgvs_c_change +
                             '\t' + hgvs_p_change + '\t' + 'NA' + '\t' + snpeff_rank + '\t' + intron + '\t' + hgnc_omimid + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + rs_id + '\t' + rs_val + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + sift_pred +
                            '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\n')


def vaf_percent(numer,denom):
    numer =float(numer)
    denom = float(denom)
    try:
        percent = (numer/(denom)) * 100
        return str("%.2f" % percent)
    except:
        return "0"

def avg_qual(ref_qual, alt_qual, ref_depth, alt_depth):
    try:
        ref_qual = float(ref_qual)
        alt_qual = float(alt_qual)
        ref_depth = float(ref_depth)
        alt_depth = float(alt_depth)
      
        try:
            avg_ref_qual = ref_qual/ref_depth
        except:
            avg_ref_qual = 0.0
        try:
            avg_alt_qual = alt_qual/alt_depth
        except:
            avg_alt_qual = 0.0
        avg_qual = ((ref_depth * avg_ref_qual) + (alt_depth * avg_alt_qual))/ (alt_depth + ref_depth)
        return str("%.2f" % avg_qual)
    except:
        traceback.print_exc()

def main():
    anno = sys.argv[1]
    anno_out=sys.argv[2]
    sample=sys.argv[3]
    onco_parse(anno,anno_out,sample)

main()

