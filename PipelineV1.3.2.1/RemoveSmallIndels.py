##V1.2
## Run this script as: python RemoveSmallIndels.py <sample.GATK.pre_filtered.vcf> <sample.GATK.pre_filtered.greaterThan6bpIndels.vcf>
## This script removes the Indels of size smaller than 6 bp.

import sys
import re
InputFile=open(sys.argv[1],'r')
OutputFile=open(sys.argv[2],'w')

InputLines=InputFile.readlines()

for eachline in InputLines:
	eachline=eachline.rstrip()
	if(re.match('#',eachline)):
		OutputFile.write(eachline+"\n")
	else:
		eachlineComponents=eachline.split("\t")
		if(len(eachlineComponents[3])>6 or len(eachlineComponents[4])>6):
			OutputFile.write(eachline+"\n")
		else:
			continue
OutputFile.close()
InputFile.close()

