##V1.2
##Run this script as: python ampliconSummary.py <SampleName>
##This script takes in the sample.Depth.below250X.amplicons file as an input and adds the sample name to the each line of that file to give the output named as sample.Depth.below250X.amplicons.summary file.


import sys
import re

InputFileName=sys.argv[1]+"/"+sys.argv[1]+".Depth.below250X.amplicons"
OutputFileName=sys.argv[1]+"/"+sys.argv[1]+".Depth.below250X.amplicons.summary"
InputFile=open(InputFileName,'r')
OutputFile=open(OutputFileName,'w')

InputFileLines=InputFile.readlines()

for eachline in InputFileLines:
	eachline=eachline.rstrip()
	newline=sys.argv[1]+"\t"+eachline
	OutputFile.write(newline+"\n")
OutputFile.close()
InputFile.close()


