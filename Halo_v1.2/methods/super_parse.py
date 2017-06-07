#f_parser.py
#
# Ashkan Bigdeli - 2/26/2017
#
#
#
# Takes in a VCF file and standerizes
# the output for cpd relevant metrics:
#                                     GT = Genotype
#                                     DP = Total Depth
#                                     RD = Reference Depth
#                                     VD = Variant Depth
#                                     VAF = Variant Allele Frequency
#                                     VC = Variant Caller
#                                     RPQ = Average Phred Quality Of Reference
#                                     VPQ = Average Phred Quality Of Alternate                                     
#
#
# Modifications to any input VCF should alter only the Format and Sample Column 
# of VCF 4.0 or greater. Output format is VCF 4.1


import sys, time



# mutect2_parse
#
# @param1 = MuTect2 VCF
# 
#
# This method will auto write an output file with a fixed extension, the parser
# works for MuTect2, Scalpel, and Pindel as run through the CPD HALO_V1.0 pipeline
def indel_parse(caller, input_vcf, out_vcf):

    #write cpd values, and format to vcf
    with open(out_vcf, 'w+') as parsed:
        # write CPD header lines, currently in vcf format 4.1
        parsed.write('##fileformat=VCFv4.1\n')
        parsed.write('##cpd_parsed=' + time.strftime("%m/%d/%Y--%H:%M:%S\n"))
        parsed.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        parsed.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n')
        parsed.write('##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Reference Depth">\n')
        parsed.write('##FORMAT=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">\n')
        parsed.write('##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Frequency">\n')
        parsed.write('##FORMAT=<ID=VC,Number=1,Type=String,Description="Variant Caller">\n')
        parsed.write('##FORMAT=<ID=PS,Number=1,Type=Float,Description="Average Phred Score Of Ref/Alt Alignment">\n')
        with open(input_vcf, 'r') as input_vcf:            
        # take care of the header lines
            for line in input_vcf:
                line = line.strip()

                if line.startswith(('##fileformat=VCF','##FORMAT','##GATK','##contig','##source')):
                    continue
                elif line.startswith(('##INFO','##FILTER','#CHROM','##reference','##fileDate','##source','##bcftools')):
                    parsed.write(line + '\n')
                else:
                    #extract and calculate values relevant to the cpd
                    line_array = line.split('\t')
                    info_col = line_array[9]
                    
                    info_array = info_col.split(':')
                    
                    genotype = info_array[0]
                    
                    depth_array = info_array[1].split(',')       
                    ref_depth = int(depth_array[0])
                    alt_depth = int(depth_array[1])                    
                    total_depth = ref_depth + alt_depth                   
                    vaf = vaf_percent(alt_depth, total_depth)
                    

                    if caller == "MuTect2":
                        #check is variant is phased, pull qual from different location if needed
                        if ("PID" in line_array[8]) or ("PTG" in line_array[8]):
    
                            phred = info_array[8]
                        else: 
                            phred = info_array[6]
 
                        phred_split = phred.split(',')
                        phred = int(phred_split[0]) + int(phred_split[1])
                    
                        avg_phred_score = avg_phred(phred,total_depth)
                    else:
                        avg_phred_score = '86'
                    #print untouched attributes + format field
                    for i in range (0,8):
                        parsed.write(line_array[i] + '\t')
                        
                    #print cpd format and values
                    parsed.write('GT:DP:RD:VD:VAF:VC:PS\t'+  genotype + ':' + str(total_depth) + ':' + str(ref_depth) + ':' + str(alt_depth) + ':' + vaf + ':' + caller + ':' + avg_phred_score + '\n')
                    

                    


# avg_phred
#
# @param1 = Summed phred quality at a given position
# @param2 = Total depth at position
#
# @return = String value, either Numeric or NA if error has occured
#
# This method is meant to calculate the summed Phred Quality into an average
def avg_phred (phred,depth):
    phred = float(phred)
    depth = float(depth)
    try:
        avg_phred = phred/depth
    except:
        return 'NA'
    return str("%.2f" % avg_phred)

# vaf_percent
#
# @param1 = Variant allele depth
# @param2 = Total depth
#
# @return = String value, percentage or NA if error has occured                    
def vaf_percent(numer,denom):
    numer =float(numer)
    denom = float(denom)
    try:
        percent = (numer/(denom)) * 100
    except:
        return "0"
    return str("%.2f" % percent)
         

def main():
    
      
      caller = sys.argv[1]
      input_vcf = sys.argv[2]
      out_vcf = sys.argv[3]

      indel_parse(caller, input_vcf, out_vcf)     
          
main()
                   


