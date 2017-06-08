##V1.2.3.a
##Run this script as: python CombineStats.py <Run_MeanCovByPositionFinal.txt> <RunStatsFinal.txt> <Run_NormalizedMeanCovByPosition.txt>
##This script normalizes each position's mean coverage (present in file Run_MeanCovByPositionFinal.txt) by dividing it by the mean coverage over all positions (present in file RunStatsFinal.txt) for each and every sample.
##The output file generated is named as Run_NormalizedMeanCovByPosition.txt.   
##This script now handles the the zero overall mean coverage entries too.
##Changes in V1.2.3.a : 1) The 'jlineCompo[9]' has been changed to 'jlineCompo[11]' to indicate the correct column ,i.e., MeanCoverage in the RunStatsFinal.txt.
##                      2) We now separate the SampleName from '_mean_cvg' string while parsing the lines of Run_HEME(or FFPE)_MeanCovByPositionFinal.txt so that the SampleName string in this file could match perfectly with the same in RunStats File when "for loops" are executed.
##                         This is because we are now sure that the SampleNames would contain '-' characters instead of '_' characters everytime within them.

import sys
infile1=open(sys.argv[1],'r')
infile2=open(sys.argv[2],'r')
outfile=open(sys.argv[3],'w')

infile1Lines=infile1.readlines()
infile2Lines=infile2.readlines()

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
	
		
	
			
					
			
