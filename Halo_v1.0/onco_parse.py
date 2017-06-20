import sys 
import traceback


# onco parse
def onco_parse(anno,anno_out):
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
                chrom = 'chr' + values[4]
                start_pos = values[5]
                end_pos = values[6]
                ref = values[10]
                alt = values[12]
                gene = values[0]
                geneid = values[1]
                onc_effect = values[9]
                transcript_loc = values[37]

                onc_gchange = values[34]
                onc_type = values[9]
                sample_name = values[16]
                rs_id = values[484]
                rs_val = values[496]

                onc_c_change = values[39]
                onc_p_change = values[41]
                onc_codon_change = values[40]
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
                    
                onc_gchange = values[34]

                hgvs_c_change = values[279]
                hgvs_p_change = values[281]
                hgvs_g_change = values[280]

                sift_pred = values[420]
                clinvar_id = values[140]
                hgnc_omimid = values[266]

                onekg_ac = values[82]
                onekg_an = values[86]
                onekg_af = values[83]
                onekg_amr = values[85]
                onekg_eur = values[93]
                onekg_afr = values[84]
                onekg_asn = values[91]

                #parse snpeff annotations from onco
                #snpeff annotations are sperated by a | in either annotation sets of
                # 15 or 30. effect transcript, p and c changes are consistently located in between pipes
                snpeff_ann = values[116]
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
                
            
                variant_class = values[8]
                

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

                depth = values[532]
                #perform checks to determine which vcf generated the annotation
                #checking for freebayes
                if (len(values[350]) >= 1):
                    if (len(values[340]) < 1):
                        alt_r1 = 0
                    else:
                        alt_r1 = int(values[340])
                    alt_r2 = int(values[343])
                    ref_r1 = int (values[350])
                    ref_r2 = int(values[352])

                    alt_depth = alt_r1 + alt_r2
                    ref_depth = ref_r1 + ref_r2
                    vaf = vaf_percent(alt_depth, ref_depth)

                    alt_qual = values[321]
                    ref_qual = values[323]
                    qual = avg_qual(ref_qual,alt_qual, ref_depth, alt_depth) 
                    if float(qual) < 20.0:
                        continue

                if (len(values[362]) > 1):
                    allelic_depth = values[362].split(',')
                    alt_depth = allelic_depth[1]
                    ref_depth = allelic_depth[0]
                    depth = int(alt_depth) + int(ref_depth)
                    vaf = vaf_percent(alt_depth, ref_depth)
                    qual = "> 25"

                if int(depth) < 50:
                    continue
                if alt_depth < 5:
                    continue
                if float(vaf) < 2.0:
                    continue
                parsed.write(sample_name + '\t' + chrom + '\t' + start_pos + '\t' + ref + '\t' + alt + '\t' + onc_type + '\t' + 'NA' + '\t' + qual + '\t' +
                             str(depth) + '\t' + 'NA' + '\t' + geneid + '\t' + snpeff_gene + '\t' + 'NA' + '\t' + snpeff_transcript + '\t' + 'NA' + '\t' + snpeff_rank +
                             '\t' + snpeff_effect + '\t' + 'NA' + '\t' + onc_codon_change + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
                             '\t' + rs_id + '\t' + snpeff_c_change + '\t' + snpeff_p_change + '\t' + location + '\t' + gene + '\t' + snpeff_effect + '\t' + snpeff_transcript +
                             '\t' + snpeff_c_change + '\t' + snpeff_p_change + '\t' + chrom + '\t' + start_pos + '\t' + rs_id + '\t' + ref + '\t' + alt + '\t' + qual +
                             '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + onekg_ac + '\t' + onekg_an + '\t' + onekg_af + '\t' + onekg_amr + '\t' +
                             onekg_asn + '\t' + onekg_eur + '\t' + onekg_afr + '\t' + onc_type + '\t' + str(depth) + '\t' + str(ref_depth) + '\t' + str(alt_depth) + '\t' + vaf + '\t'
                             + snpeff_gene + '\t' + geneid + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' + '\t' + 'NA' +
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
        percent = (numer/(numer+denom)) * 100
        return str("%.2f" % percent)
    except:
        return "0"

def avg_qual (ref_qual, alt_qual, ref_depth, alt_depth):
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
    onco_parse(anno,anno_out)

main()

