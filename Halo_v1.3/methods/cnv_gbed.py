# Ashkan Bigdeli
# 2/21/2018
#
# A help program to collapse cnv output on chromosome for graphing purposes

import sys
from collections import OrderedDict

# takes cnv input and writes a bed file for cnvkit.py's graphing function
def collapse_chrom(cnv_in, bed_out):
    with open(cnv_in, 'r') as cnv_in:
        with open(bed_out, 'w') as bed_out:
            cnv_in.readline() # skip header
            bd = OrderedDict()
            for line in cnv_in:
                values = line.split("\t")
                chrom = values[2]
                start = values[3]
                stop  = values[4]
                
                if chrom not in bd:
                    bd[chrom] = [start,stop]
                else:
                    bd[chrom][1] = stop

            for key, value in bd.items():
                bed_out.write(key + "\t" + value[0] + "\t" + value[1] + "\n")
            
def main():
    cnv_in = sys.argv[1]
    bed_out=sys.argv[2]
    collapse_chrom(cnv_in,bed_out)

main()

