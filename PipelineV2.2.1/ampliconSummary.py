##V2.1.1_hpc
##Run this script as: python ampliconSummary.py <SampleFolderPath>
##This script takes in the sample.Depth.below250X(and 150X).amplicons files as an input and adds the sample name to the each line of that file to give the output named as sample.Depth.below250X.amplicons.summary file.
##Changes in V2.0 : This script now also generates Depth.below150X.amplicons.summary. The description of the script is changed accordingly.
##Changes in V2.1.1_hpc: 1) Changed the argument on command line from <SampleName> to <SampleFolderPath>. Added variables "SampleFolderPath" and "Sample_name" which are assigned values as <SampleFolderPath> and "Sample_name extracted from <SampleFolderPath> respectively.

import sys
import re
SampleFolderPath=sys.argv[1]
Sample_name=SampleFolderPath.split("/")[5]

InputFileName=SampleFolderPath+"/"+Sample_name+".Depth.below250X.amplicons"
OutputFileName=SampleFolderPath+"/"+Sample_name+".Depth.below250X.amplicons.summary"
InputFileName2=SampleFolderPath+"/"+Sample_name+".Depth.below150X.amplicons"
OutputFileName2=SampleFolderPath+"/"+Sample_name+".Depth.below150X.amplicons.summary"
InputFile=open(InputFileName,'r')
OutputFile=open(OutputFileName,'w')
InputFile2=open(InputFileName2,'r')
OutputFile2=open(OutputFileName2,'w')

InputFileLines=InputFile.readlines()
InputFileLines2=InputFile2.readlines()

for eachline in InputFileLines:
	eachline=eachline.rstrip()
	newline=Sample_name+"\t"+eachline
	OutputFile.write(newline+"\n")
for eachline2 in InputFileLines2:
	eachline2=eachline2.rstrip()
	newline2=Sample_name+"\t"+eachline2
	OutputFile2.write(newline2+"\n")
OutputFile.close()
InputFile.close()
OutputFile2.close()
InputFile2.close()



