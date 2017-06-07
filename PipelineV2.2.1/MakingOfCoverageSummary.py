##V2.1.1_hpc
## Run this script as : python MakingOfCoverageSummary.py <SampleFolderPath> <Description>
## This script converts the table (or matrix) format of the Depth.sample_interval_summary into a transpose format where the chromosome intervals become the column names (in the first row) followed by the mean coverage values for the corresponding intervals in the second row.
## Changed in V1.3 : 1) This script now generates panel specific (HEME or FFPE) MeanCovByPosition files.This was made possible by adding a '<Description>' argument to the command line and changing the 'outfileName' variable in the script accordingly.
##		       2) Explaination is little bit modified in the starting of the script to give more clear idea to the readers.
##Changes in V2.1.1_hpc  : 1) The argument <SampleName> is changed to <SampleFolderPath> on command line. Added variables "SampleFolderPath" and "Sample_name" and assigned them the values of <SampleFolderPath> argument and "Sample_name extracted from <SampleFolderPath>" respectively.
##			   2) Added variable "Descript" to which the value of <Description> argument.

import sys
SampleFolderPath=sys.argv[1]
Descript=sys.argv[2]
Sample_name=SampleFolderPath.split("/")[5]
infileName=SampleFolderPath+"/"+Sample_name+".Depth.sample_interval_summary" ## Opens up Sample's Depth.sample_interval_summary file.
outfileName=SampleFolderPath+"/"+Sample_name+"."+Descript+"."+"MeanCovByPosition"  ## Specifying the output file.
infile=open(infileName,'r')
outfile=open(outfileName,'w')
infileLines=infile.readlines()
matrix=[]    ## Inializing the matrix.


for eachline in infileLines:
	eachline=eachline.rstrip()
	Components=eachline.split("\t")
	Position=Components[0]
	MeanCov=Components[4]
	matrix.append((Position,MeanCov))    ## Storing the lines of Sample's Depth.sample_interval_summary file in the form of matrix of rows and columns.
transpose=zip(*matrix)			     ## Changing rows to columns and columns to rows for the data in matrix so that the output looks better. This is called 'transpose of a matrix'.	
NumRowsInTranspose=len(transpose)            ## Determing and storing of number of rows from the new transpose matrix for running the loop on each row further.

for i in range(0,NumRowsInTranspose):		## Running for loop to print the final transpose format in the output file.
	for j in range(0,len(matrix)):
		outfile.write(str(transpose[i][j])+"\t")
        outfile.write("\n")

outfile.close()
infile.close()
	

	

