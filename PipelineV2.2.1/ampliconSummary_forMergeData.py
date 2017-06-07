##V2.0
##Run this script as: python ampliconSummary_forMergeData.py <MergeDirectoryPath>
##This script takes in the sample.Depth.below250X.amplicons_sorted file as an input and adds the sample name to the each line of that file to give the output named as sample.Depth.below250X.amplicons_sorted.summary file.
## Changes in V2.0 : Added lines to make Depth.below.150X.amplicons_sorted.summary


import sys
import re

MergeDirectoryPath=sys.argv[1]
MergeDirectoryName=sys.argv[1].split("/")[5]

InputFileName=MergeDirectoryPath+"/"+MergeDirectoryName+".Depth.below250X.amplicons_sorted"
OutputFileName=MergeDirectoryPath+"/"+MergeDirectoryName+".Depth.below250X.amplicons_sorted.summary"
InputFileName2=MergeDirectoryPath+"/"+MergeDirectoryName+".Depth.below150X.amplicons_sorted"
OutputFileName2=MergeDirectoryPath+"/"+MergeDirectoryName+".Depth.below150X.amplicons_sorted.summary"
InputFile=open(InputFileName,'r')
OutputFile=open(OutputFileName,'w')
InputFile2=open(InputFileName2,'r')
OutputFile2=open(OutputFileName2,'w')


InputFileLines=InputFile.readlines()
InputFileLines2=InputFile2.readlines()

for eachline in InputFileLines:
	eachline=eachline.rstrip()
	newline=MergeDirectoryName+"\t"+eachline
	OutputFile.write(newline+"\n")
for eachline2 in InputFileLines2:
	eachline2=eachline2.rstrip()
	newline2=MergeDirectoryName+"\t"+eachline2
	OutputFile2.write(newline2+"\n")
OutputFile.close()
InputFile.close()
OutputFile2.close()
InputFile2.close()



