##V2.0
## Run this script as python RecordsProcessing.py <pileupAndDepthJoined file> <GenotypeGivenAlleles file to be created> <AmpliconsWithDepthBelow250X file to be created> <AmpliconsDepthBelow100X file to be created> <DepthForAnnotation file to be created>
## v1.2 of this pipe now does not report variants at locations with less than 100 reads.
## Changes in V2.0:   1) 4th argument on the command line is now <AmpliconsDepthBelow150X file to be created>. Hence, the file handle for this file is stated in the script as 'OutFile3'. The file handle for <DepthForAnnotation file> would now be 'OutFile4' (Changes are mades accordingly in the write statements meant for this file handle).
##		      2) "If" statement is added to the script to help generate 'AmpliconsDepthBelow150X file'.
##                    3) Now the script prints out depth values for the regions reported in 'AmpliconsDepthBelow250X file' and 'AmpliconsDepthBelow150X file'. 



import sys,re,string
InFile=open(sys.argv[1],'r')   ## Opens up the input file for reading and Output files for writing.
OutFile=open(sys.argv[2],'w')
OutFile2=open(sys.argv[3],'w')
OutFile3=open(sys.argv[4],'w')
OutFile4=open(sys.argv[5],'w')
OutFile.write("##fileformat=VCFv4.1"+"\n")  ## Writes headers for each Ouptut vcf file.
OutFile.write("#CHROM"+"\t"+"POS"+"\t"+"ID"+"\t"+"REF"+"\t"+"ALT"+"\t"+"QUAL"+"\t"+"FILTER"+"\t"+"INFO"+"\t"+"FORMAT"+"\t"+"SAMPLE"+"\n")
OutFile2.write("##fileformat=VCFv4.1"+"\n")
OutFile2.write("#CHROM"+"\t"+"POS"+"\t"+"ID"+"\t"+"REF"+"\t"+"ALT"+"\t"+"QUAL"+"\t"+"FILTER"+"\t"+"INFO"+"\t"+"FORMAT"+"\t"+"SAMPLE"+"\n")
OutFile3.write("##fileformat=VCFv4.1"+"\n")
OutFile3.write("#CHROM"+"\t"+"POS"+"\t"+"ID"+"\t"+"REF"+"\t"+"ALT"+"\t"+"QUAL"+"\t"+"FILTER"+"\t"+"INFO"+"\t"+"FORMAT"+"\t"+"SAMPLE"+"\n")
OutFile4.write("##fileformat=VCFv4.1"+"\n")
OutFile4.write("#CHROM"+"\t"+"POS"+"\t"+"ID"+"\t"+"REF"+"\t"+"ALT"+"\t"+"QUAL"+"\t"+"FILTER"+"\t"+"INFO"+"\t"+"FORMAT"+"\t"+"SAMPLE"+"\n")
Lines=InFile.readlines()     ## Reads the lines from input file and store them in a List.
for eachline in Lines:       ## Runs the for loop on List of input lines.
	eachline=eachline.rstrip()
	Components=eachline.split(",")  ## Each Line gets splitted based on comma delimiter.
	RefBase=Components[1]     ## Indices of the list start from 0 and go upto (length_of_List-1). Here, the elements of the "Components" list such as 'refbase', 'total-depth at a position' and 'chromosome position' gets extracted.
	AltBaseList=list()        ## Other lists which would be used later in script gets initialized.
	AltBaseDepthList=list()
	AltBaseDepthPercentList=list()
	AltBaseFinal=""          ## Initializes Variables. 
	AltBaseDepthFinal=""
	AltBaseDepthPercentFinal=""
	TotalDepthAtPosition=Components[2]
	ChromPos=Components[0].split(":")
	if(int(TotalDepthAtPosition)<250):   ## The chromosome positions having total read depth below 250 are extracted.
		OutFile2.write(ChromPos[0]+"\t"+ChromPos[1]+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"FDP="+TotalDepthAtPosition+"\t"+"."+"\t"+"."+"\n")
	if(int(TotalDepthAtPosition)<150):
		OutFile3.write(ChromPos[0]+"\t"+ChromPos[1]+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"."+"\t"+"FDP="+TotalDepthAtPosition+"\t"+"."+"\t"+"."+"\n")
	##if(int(TotalDepthAtPosition)>=50):
	BaseDepths=Components[3].split(" ")  ## Reads aligned under a specific chromosome location may contain A,C,G or T alleles. Therefore, extracting each base alongwith its depth at that particular position.
	FirstBase_DepthParts=BaseDepths[0].split(":")
	SecondBase_DepthParts=BaseDepths[1].split(":")
	ThirdBase_DepthParts=BaseDepths[2].split(":")
	FourthBase_DepthParts=BaseDepths[3].split(":")
	##FifthBase_DepthParts=BaseDepths[4].split(":")
	##now extract exact depth number for all bases
	if(RefBase == FirstBase_DepthParts[0]):    ## Checking if Reference Base is equal to one of the A,G,C and T so that we could get depth for the ref base at that particular chromosome position.
		RefBaseDepth=FirstBase_DepthParts[1]
		AltBase1=SecondBase_DepthParts[0]  ## Selecting the bases other than reference to be the alternate ones and extracting their respective depths too.
		AltBase1Depth=SecondBase_DepthParts[1]
		AltBase2=ThirdBase_DepthParts[0]
		AltBase2Depth=ThirdBase_DepthParts[1]
		AltBase3=FourthBase_DepthParts[0]
		AltBase3Depth=FourthBase_DepthParts[1]
		
	elif(RefBase == SecondBase_DepthParts[0]):
		RefBaseDepth=SecondBase_DepthParts[1]
		AltBase1=FirstBase_DepthParts[0]
                AltBase1Depth=FirstBase_DepthParts[1]
                AltBase2=ThirdBase_DepthParts[0]
               	AltBase2Depth=ThirdBase_DepthParts[1]
                AltBase3=FourthBase_DepthParts[0]
               	AltBase3Depth=FourthBase_DepthParts[1]
		
	elif(RefBase == ThirdBase_DepthParts[0]):
		RefBaseDepth=ThirdBase_DepthParts[1]
		AltBase1=FirstBase_DepthParts[0]
               	AltBase1Depth=FirstBase_DepthParts[1]
               	AltBase2=SecondBase_DepthParts[0]
               	AltBase2Depth=SecondBase_DepthParts[1]
               	AltBase3=FourthBase_DepthParts[0]
               	AltBase3Depth=FourthBase_DepthParts[1]
		
	elif(RefBase == FourthBase_DepthParts[0]):
		RefBaseDepth=FourthBase_DepthParts[1]
		AltBase1=FirstBase_DepthParts[0]
               	AltBase1Depth=FirstBase_DepthParts[1]
                AltBase2=SecondBase_DepthParts[0]
               	AltBase2Depth=SecondBase_DepthParts[1]
                AltBase3=ThirdBase_DepthParts[0]
                AltBase3Depth=ThirdBase_DepthParts[1]
			
	##elif(Components[1] == FifthBase_DepthParts[0]):
	##RefBaseDepthPercent=-1.00
	if(int(TotalDepthAtPosition)>0): ## Checks if there are some reads aligned at a particular position (TotalDepth>0). If yes,then it calculates the Percentage of Depth for ref base and all other possible alternate alleles.Finally prints such positions alongwith their ref and possible alternate alleles info in the DepthForAnnotation File.
		RefBaseDepthPercent=100 * float(RefBaseDepth) / float(TotalDepthAtPosition)
		AltBase1DepthPercent=100 * float(AltBase1Depth) / float(TotalDepthAtPosition)
		AltBase2DepthPercent=100 * float(AltBase2Depth) / float(TotalDepthAtPosition)
                AltBase3DepthPercent=100 * float(AltBase3Depth) / float(TotalDepthAtPosition)
		OutFile4.write(ChromPos[0]+"\t"+ChromPos[1]+"\t"+"."+"\t"+RefBase+"\t"+AltBase1+"\t"+"."+"\t"+"."+"\t"+"FDP="+TotalDepthAtPosition+";FRD="+RefBaseDepth+";FAD="+AltBase1Depth+";FAF="+"{0:.2f}".format(AltBase1DepthPercent)+"\t"+"."+"\t"+"."+"\n")
		OutFile4.write(ChromPos[0]+"\t"+ChromPos[1]+"\t"+"."+"\t"+RefBase+"\t"+AltBase2+"\t"+"."+"\t"+"."+"\t"+"FDP="+TotalDepthAtPosition+";FRD="+RefBaseDepth+";FAD="+AltBase2Depth+";FAF="+"{0:.2f}".format(AltBase2DepthPercent)+"\t"+"."+"\t"+"."+"\n")
		OutFile4.write(ChromPos[0]+"\t"+ChromPos[1]+"\t"+"."+"\t"+RefBase+"\t"+AltBase3+"\t"+"."+"\t"+"."+"\t"+"FDP="+TotalDepthAtPosition+";FRD="+RefBaseDepth+";FAD="+AltBase3Depth+";FAF="+"{0:.2f}".format(AltBase3DepthPercent)+"\t"+"."+"\t"+"."+"\n")
	if(int(TotalDepthAtPosition)>=100):   ##GenotypeGivenAlleles file is prepared. All those chromosome positions with the TotalDepth>=100 (this has changed since v 1.0-1.1 which called 50 bp or more), RefBaseDepthPercent <=96 and AltBasesDepthPercent >=4.00 are printed. Here, the Percentage Depth is checked for each alternate base (to be >= 4.00) and only those alternate bases are printed in the output file.
		if(RefBaseDepthPercent<=96.00):
			if(AltBase1DepthPercent>=4.00):
				AltBaseList.append(AltBase1)
				AltBaseDepthList.append(AltBase1Depth)
				AltBaseDepthPercentList.append("{0:.2f}".format(AltBase1DepthPercent))
			if(AltBase2DepthPercent>=4.00):
				AltBaseList.append(AltBase2)
				AltBaseDepthList.append(AltBase2Depth)
				AltBaseDepthPercentList.append("{0:.2f}".format(AltBase2DepthPercent))
			if(AltBase3DepthPercent>=4.00):
				AltBaseList.append(AltBase3)
				AltBaseDepthList.append(AltBase3Depth)
				AltBaseDepthPercentList.append("{0:.2f}".format(AltBase3DepthPercent))
			AltBaseFinal=string.join(AltBaseList,",")
			AltBaseDepthFinal=string.join(AltBaseDepthList,",")
			AltBaseDepthPercentFinal=string.join(AltBaseDepthPercentList,",")
			if(len(AltBaseFinal)!=0):
				OutFile.write(ChromPos[0]+"\t"+ChromPos[1]+"\t"+"."+"\t"+RefBase+"\t"+AltBaseFinal+"\t"+"."+"\t"+"."+"\t"+"FDP="+TotalDepthAtPosition+";FRD="+RefBaseDepth+";FAD="+AltBaseDepthFinal+";FAF="+AltBaseDepthPercentFinal+"\t"+"."+"\t"+"."+"\n")
        
OutFile.close()
OutFile2.close()
OutFile3.close()
InFile.close()		
	
	
	
