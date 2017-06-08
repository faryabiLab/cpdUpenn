##V1.0
## Run this script as: python MultiAllelic.py <input vcf file> <output vcf file>

## This script splits the multiallelic positions in a vcf file into single allele positions.

import sys,re
FHandle=open(sys.argv[1],'r')
OHandle=open(sys.argv[2],'w')
Lines=FHandle.readlines()
for eachline in Lines:
	eachline=eachline.rstrip()
	if(re.match("#",eachline)):
		OHandle.write(eachline+"\n")
	else:
		Components=eachline.split("\t")
		if(re.search(",",Components[4])):
			Alleles=Components[4].split(",")
			NumAlleles=len(Alleles)
			for i in range(0,NumAlleles):
				OHandle.write(Components[0]+"\t"+Components[1]+"\t"+Components[2]+"\t"+Components[3]+"\t"+Alleles[i]+"\t"+Components[5]+"\t"+Components[6]+"\t"+Components[7]+"\t"+Components[8]+"\t"+Components[9]+"\n")
		else:
			OHandle.write(Components[0]+"\t"+Components[1]+"\t"+Components[2]+"\t"+Components[3]+"\t"+Components[4]+"\t"+Components[5]+"\t"+Components[6]+"\t"+Components[7]+"\t"+Components[8]+"\t"+Components[9]+"\n")


OHandle.close()
FHandle.close()		
