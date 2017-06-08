##V2.0
##Run this script as: python CombineStats.py <Run_HEME(or FFPE or mPPP or bPPP or CEBPa)_MeanCovByPositionFinal.txt> <RunStatsFinal.txt> <Run_HEME(or FFPE or mPPP or bPPP or CEBPa)_NormalizedMeanCovByPosition.txt>
##This script, normalizes each chromosome interval's mean coverage (present in file Run_HEME( or FFPE or mPPP or bPPP or CEBPa)_MeanCovByPositionFinal.txt) by dividing it by the mean coverage over all positions (present in file RunStatsFinal.txt), for each and every sample.
##The output file generated is named as Run_NormalizedMeanCovByPosition.txt.
##Changes in V1.3 : 1) The 'jlineCompo[9]' has been changed to 'jlineCompo[11]' to indicate the correct column ,i.e., MeanCoverage in the RunStatsFinal.txt. 
##                  2) We now separate the SampleName from '_mean_cvg' string while parsing the lines of MeanCovByPositionFinal.txt so that the SampleName string in this file could match perfectly with the same in RunStatsFinal File when "for loops" are executed.
##                     This is because we are now sure that the SampleNames would contain '-' characters instead of '_' characters everytime within them.The '_' character separates 'mean_cvg' string from the original SampleName.
##		    3) This script now handles the the zero overall mean coverage entries. It performs no division when StoredMeanCov=0. Instead it outputs the normalized values as '0s' in the final NormalizedMeanCovByPosition.txt file.
##                  4) This script now uses "flag" variable to flag the following situations: a) When a sample entry is not present in RunStatsFinal.txt, the flag value remains to be zero. b) When a match is found for the
##                     SampleName in RunStatsFinal.txt, the flag value equates to '1'. Further, for flag value '0' the StoredMeanCov is equated to be '0'. And if StoredMeanCov is '0', the script does not perfrom division and
##                     outputs the normalized values as '0s' in the final NormalizedMeanCovByPosition.txt file.
##Changes in V1.3.3 :1) Replaced integer 0 with string "0" in the statements 'StoredMeanCov=0' and in the if statement 'if(StoredMeanCov==0)'.
##Changes in V1.4.2 :1) include "or StoredMeanCov==0.00" in the 'if(StoredMeanCov=="0")' statement to check for an extra possibility.
##Changes in V2.0:   1) Changed "StoredMeanCov=0.00" to "StoredMeanCov="0.00"".

import sys
infile1=open(sys.argv[1],'r')
infile2=open(sys.argv[2],'r')
outfile=open(sys.argv[3],'w')

infile1Lines=infile1.readlines()
infile2Lines=infile2.readlines()

##Starting loop over lines in the file Run_HEME( or FFPE or mPPP or bPPP or CEBPa)_MeanCovByPositionFinal.txt and for each such line another inner-for-loop is started to find the matching SampleName in RunStatsFinal.txt file.
##This is done to extract the sample's mean coverage from RunStatsFinal file and using the same to normalize the interval coverages in respective MeanCovByPositionFinal.txt file.
##Normalization is nothing but dividing a interval's mean coverage by a sample's mean coverage (i.e. mean coverage over all positions).All Normalized Interval Mean Coverages are then written to the final output.


for iLine in infile1Lines:
	flag=0
	iLine=iLine.rstrip()
	iLineCompo=iLine.split("\t")
        iLineSample=iLineCompo[0].split("_")[0]
	if(iLineSample=="Target"):
		outfile.write(iLine+"\n")
	else:
		for jLine in infile2Lines:
			jLine=jLine.rstrip()
			jLineCompo=jLine.split("\t")
			jLineSample=jLineCompo[0]
			jLineMeanCov=jLineCompo[11]
			if(iLineSample==jLineSample):
				StoredMeanCov=jLineMeanCov
				flag=1
				break
		if(flag==0):
			StoredMeanCov="0"
		if(StoredMeanCov=="0" or StoredMeanCov=="0.00"):
			for i in range(1,len(iLineCompo)):
				iLineCompo[i]="0"
		else:
			for i in range(1,len(iLineCompo)):
				iLineCompo[i]=str(float(iLineCompo[i])/float(StoredMeanCov))
			
		iLine="\t".join(iLineCompo[0:])
		outfile.write(iLine+"\n")

outfile.close()
infile2.close()
infile1.close()
	
		
	
			
					
			
