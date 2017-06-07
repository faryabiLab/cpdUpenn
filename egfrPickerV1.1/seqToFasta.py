##V1.1
## Run this script as : python seqToFasta.py <input_sequence_text_file> <output_fasta_sequence_file>
import sys
infile=open(sys.argv[1],'r')
outfile=open(sys.argv[2],'w')
infileLines=infile.readlines()
count=1
for eachline in infileLines:
	eachline=eachline.rstrip()
	outfile.write(">"+str(count)+"\n"+eachline+"\n")
	count=count+1
outfile.close()
infile.close()

