##V1.4
##Run this script as: python ampliconSummary_forMergeData.py <SampleName>
##This script takes in the sample.Depth.below250X.amplicons_sorted file as an input and adds the sample name to the each line of that file to give the output named as sample.Depth.below250X.amplicons_sorted.summary file.


import sys
import re

InputFileName=sys.argv[1]+"/"+sys.argv[1]+".Depth.below250X.amplicons_sorted"
OutputFileName=sys.argv[1]+"/"+sys.argv[1]+".Depth.below250X.amplicons_sorted.summary"
InputFile=open(InputFileName,'r')
OutputFile=open(OutputFileName,'w')

InputFileLines=InputFile.readlines()

for eachline in InputFileLines:
	eachline=eachline.rstrip()
	newline=sys.argv[1]+"\t"+eachline
	OutputFile.write(newline+"\n")
OutputFile.close()
InputFile.close()


