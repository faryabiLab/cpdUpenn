##V1.2.4
##Run this script as: python CombineStats.py <Run_HEME(or FFPE)_MeanCovByPositionFinal.txt> <RunStatsFinal.txt> <Run_HEME(or FFPE)_NormalizedMeanCovByPosition.txt>
##This script, normalizes each chromosome interval's mean coverage (present in file Run_HEME( or FFPE)_MeanCovByPositionFinal.txt) by dividing it by the mean coverage over all positions (present in file RunStatsFinal.txt), for each and every sample.
##The output file generated is named as Run_NormalizedMeanCovByPosition.txt.
##Changes in V1.2.4 : 1) The 'jlineCompo[9]' has been changed to 'jlineCompo[11]' to indicate the correct column ,i.e., MeanCoverage in the RunStatsFinal.txt. 
##                  2) We now separate the SampleName from '_mean_cvg' string while parsing the lines of Run_HEME(or FFPE)_MeanCovByPositionFinal.txt so that the SampleName string in this file could match perfectly with the same in RunStats File when "for loops" are executed.
##                     This is because we are now sure that the SampleNames would contain '-' characters instead of '_' characters everytime within them.
##		    3) This script now handles the the zero overall mean coverage entries too.

import sys
infile1=open(sys.argv[1],'r')
infile2=open(sys.argv[2],'r')
outfile=open(sys.argv[3],'w')

infile1Lines=infile1.readlines()
infile2Lines=infile2.readlines()

##Starting loop over lines in the file Run_HEME( or FFPE)_MeanCovByPositionFinal.txt and for each such line there is another for-loop started to find the matching SampleName in RunStatsFinal.txt file.
##This is done to extract the sample's mean coverage from RunStatsFinal file and using the same to normalize the interval coverages in Run_HEME( or FFPE)_MeanCovByPositionFinal.txt for the respective sample.
##Normalization is nothing but dividing a interval's mean coverage by a sample's mean coverage (i.e. mean coverage over all positions).All Normalized Interval Mean Coverages are then written to the final output.


for iLine in infile1Lines:
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
				break
		if(StoredMeanCov==0):
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
	
		
	
			
					
			
