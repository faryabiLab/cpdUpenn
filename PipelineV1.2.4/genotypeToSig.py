##V1.2.4
## Run this script as: python genotypeToSig.py <Sample.signatureGenotype.vcf> <Sample.SignatureFile> <Description>
## This Script collects the Sample's genotype info at all the 'bi-allelic' sites (common dbsnp137 polymorphic sites) from Sample.signatureGenotype.vcf
## and then prints a signature in an output file (Sample.SignatureFile) which would represent a sample's identity.
## Changes in V1.2.4 : 1) 'Description' parameter was added to the commandline .
##		       2) An alphabet (F for FFPE and H for HEME) which is extracted from the 'Description' string is added at the start of a SignatureString.



import sys
import re

Input=open(sys.argv[1],'r')
Output=open(sys.argv[2],'w')
Panel=sys.argv[3]
SignatureString=Panel[0]

while 1 :
	InputLine=Input.readline()
	if(InputLine):
		if(re.match('#',InputLine)):
			continue
		else:
			InputLine=InputLine.rstrip()
			InputCompo=InputLine.split("\t")
			Genotype=InputCompo[9].split(":")[0]
			## Eliminate those positions which were not genotyped with the help of 'if' condition below. Such positions have genotype './.'
			if(re.match('\.',Genotype)):
				continue
			else:
				GenoNum=int(Genotype.split("/")[0])+int(Genotype.split("/")[1]) ## Creating a customized Genotype number(GenoNum) by adding up the Genotype Components. For Example: for genotype '0/1' the GenoNum=0+1=1, likewise for genotype '1/1' the GenoNum=2 and for '0/0' GenoNum=0
				SignatureString=SignatureString+str(GenoNum)
	else:
		Output.write(SignatureString+"\n")
		break

Output.close()
Input.close()
		
			
			
	
