##V1.0
## Run this script as : python FilterAlignScore.py <sam-file to be filtered using minimum alignment score cut-off of 95> <output sam-file>



import re, sys
FileHandle=open(sys.argv[1],'r')  ## Opens the Input File for reading.
OutFileHandle=open(sys.argv[2],'w')  ## Opens the Output File for writing.
FileLines=FileHandle.readlines()   ## Reads Lines from the Input File Handle and stores them in a list form.
for eachline in FileLines:         ## Starting a for loop to go through and process each line stored in the list.
	eachline=eachline.rstrip() ## removes all non-whitespace characters such as "\n" at the end of the line.
	if(re.match('@',eachline)):   ## Checks for character '@' at the start of a line. This is done to detect each header line.
		OutFileHandle.write(eachline + "\n")     ## Write an header to the Output File. 
	else:
		LineComponents=eachline.split("\t")  ## This is for non-header lines,i.e. main sam file lines.It splits each such line into its fields based on a tab-delimiter and stores each field into a List.
		CIGAR=LineComponents[5]              ## Stores the CIGAR string for the read-alignment in a variable
		CIGAR=re.sub('\d+S','',CIGAR)        ## Substitutes all Ss(Soft-Clipped Bases), Hs (Hard-Clipped Bases),Ns (skipped reference bases) and Ps (padded reference bases) in CIGAR alongwith their preceeding digits, with an empty string. This is done as we dont want to include them for the calculation of the alignment length.
		CIGAR=re.sub('\d+H','',CIGAR)
		CIGAR=re.sub('\d+N','',CIGAR)
		CIGAR=re.sub('\d+P','',CIGAR)
		CIGARComponents=re.split("\D",CIGAR)   ## Splits the CIGAR string at the points where it finds characters other than digits. This is done to obtain the number of matches, mismatches and indels for the calculation of the alignment length.
		TotalAlignLength=0                     ## For each read-alignment, at the start of alignment length calculation, it initializes the Total Alignment Lenght to be equal to zero.
		for eachcompo in CIGARComponents:      ## For each CIGAR string digits stored , it checks whether its actually a digit and then adds that digit into the calculation of alignment length.
			if(eachcompo.isdigit()):
				TotalAlignLength=TotalAlignLength+int(eachcompo)
		
		for component in LineComponents:       ## Runs the for loop on the list which contains its elements as the fields of the non-header line.
			if(re.match('MD:Z:',component)): ## Detects the MD tag (field) while parsing through the list and processes that tag further.
				TagContents=component.split(":")  ## Spits the MD tag based on semi-colon delimiter to access the real MD tag (containing no. of matches, mismatches and Indels).
				MatchInfoNumbers=re.split("\D",TagContents[2]) ## Splits the real MD tag futher on the non-digit characters to access only the no. of matches. Only MD tag contains the no. of matches separate from no. of mismatches and that is not the same case with CIGAR string. CIGAR string does not differentiate between matches and mismatches in its representation.
		                MatchSum=0                                     ## Intializes the sum of matches to be equal to zero.
				for eachvalue in MatchInfoNumbers:             ## All the matches are added to give the sum.
					if(eachvalue.isdigit()):
						MatchSum=MatchSum+int(eachvalue)
		AlignScore=100 * MatchSum/TotalAlignLength                   ## Alignment Score i.e. Match Percentage is calculated.
		if(AlignScore>=95):                                          ## Only those read-alignments that have alignment score greater than or equal to 95, are printed in the output file.
			OutFileHandle.write(eachline + "\n")
		
	
FileHandle.close()                                                          ## Closing the Input and Output Files.
OutFileHandle.close()
		
		

