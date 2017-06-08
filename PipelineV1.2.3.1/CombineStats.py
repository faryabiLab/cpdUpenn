##V1.0
##Run this script as: python CombineStats.py <Run_MeanCovByPositionFinal.txt> <RunStatsFinal.txt> <Run_NormalizedMeanCovByPosition.txt>
##This script normalizes each position's mean coverage (present in file Run_MeanCovByPositionFinal.txt) by dividing it by the mean coverage over all positions (present in file RunStatsFinal.txt) for each and every sample.
##The output file generated is named as Run_NormalizedMeanCovByPosition.txt.   

import sys
infile1=open(sys.argv[1],'r')
infile2=open(sys.argv[2],'r')
outfile=open(sys.argv[3],'w')

infile1Lines=infile1.readlines()
infile2Lines=infile2.readlines()

for iLine in infile1Lines:
	iLine=iLine.rstrip()
	iLineCompo=iLine.split("\t")
        iLineSample=iLineCompo[0]
	if(iLineSample=="Target"):
		outfile.write(iLine+"\n")
	else:
		for jLine in infile2Lines:
			jLine=jLine.rstrip()
			jLineCompo=jLine.split("\t")
			jLineSample=jLineCompo[0]+"_mean_cvg"
			jLineMeanCov=jLineCompo[9]
			if(iLineSample==jLineSample):
				StoredMeanCov=jLineMeanCov
				break
		for i in range(1,len(iLineCompo)):
			iLineCompo[i]=str(float(iLineCompo[i])/float(StoredMeanCov))
			
		iLine="\t".join(iLineCompo[0:])
		outfile.write(iLine+"\n")
	
		
	
			
					
			
