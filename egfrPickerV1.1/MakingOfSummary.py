##V1.0
## Run this script as : python MakingOfSummary.py <SampleName>
## This script stores the format of the Sample.AmpCount.awked.out file in a matrix, makes a transpose of that matrix and prints that transpose to the Sample.summary file.  

import sys
SampleName=sys.argv[1]
infileName=SampleName+".AmpCount.awked.out" ## Opens up Sample's AmpCount.awked.out file.
outfileName=SampleName+".summary"  ## Specifying the output file.
infile=open(infileName,'r')
outfile=open(outfileName,'w')
infileLines=infile.readlines()
matrix=[]    ## Inializing the matrix.

matrix.append(("SampleName",SampleName))
for eachline in infileLines:
	eachline=eachline.rstrip()
	Components=eachline.split(" ")
	Amplicon=Components[0]
	Count=Components[1]
	matrix.append((Amplicon,Count))    ## Storing the lines of Sample's AmpCount.awked.out file in the form of matrix of rows and columns.
transpose=zip(*matrix)			     ## Changing rows to columns and columns to rows for the data in matrix so that the output looks better. This is called 'transpose of a matrix'.	
NumRowsInTranspose=len(transpose)            ## Determing and storing of number of rows from the new transpose matrix for running the loop on each row further.

for i in range(0,NumRowsInTranspose):		## Running for loop to print the final transpose format in the output file.
	for j in range(0,len(matrix)):
		outfile.write(str(transpose[i][j])+"\t")
        outfile.write("\n")

outfile.close()
infile.close()
	

	

