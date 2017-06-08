##V1.0
## Run this script as : python MakingOfCoverageSummary.py <SampleName>
## This script converts the format of the Depth.sample_interval_summary to display the rows (chromosome positions) as columns and the column (Sample's Mean Coverage at each position) as a row.  

import sys
infileName=sys.argv[1]+"/"+sys.argv[1]+".Depth.sample_interval_summary" ## Opens up Sample's Depth.sample_interval_summary file.
outfileName=sys.argv[1]+"/"+sys.argv[1]+".MeanCovByPosition"  ## Specifying the output file.
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
	

	

