##V2.0
##Run this script as: python ampliconSummary.py <SampleName>
##This script takes in the sample.Depth.below250X(and 150X).amplicons files as an input and adds the sample name to the each line of that file to give the output named as sample.Depth.below250X.amplicons.summary file.
##Changes in V2.0 : This script now also generates Depth.below150X.amplicons.summary. The description of the script is changed accordingly.

import sys
import re

InputFileName=sys.argv[1]+"/"+sys.argv[1]+".Depth.below250X.amplicons"
OutputFileName=sys.argv[1]+"/"+sys.argv[1]+".Depth.below250X.amplicons.summary"
InputFileName2=sys.argv[1]+"/"+sys.argv[1]+".Depth.below150X.amplicons"
OutputFileName2=sys.argv[1]+"/"+sys.argv[1]+".Depth.below150X.amplicons.summary"
InputFile=open(InputFileName,'r')
OutputFile=open(OutputFileName,'w')
InputFile2=open(InputFileName2,'r')
OutputFile2=open(OutputFileName2,'w')

InputFileLines=InputFile.readlines()
InputFileLines2=InputFile2.readlines()

for eachline in InputFileLines:
	eachline=eachline.rstrip()
	newline=sys.argv[1]+"\t"+eachline
	OutputFile.write(newline+"\n")
for eachline2 in InputFileLines2:
	eachline2=eachline2.rstrip()
	newline2=sys.argv[1]+"\t"+eachline2
	OutputFile2.write(newline2+"\n")
OutputFile.close()
InputFile.close()
OutputFile2.close()
InputFile2.close()



