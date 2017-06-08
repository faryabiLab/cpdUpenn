##V1.0
## Run this script as: python SummaryStats.py <SampleName>
## This script collects the statistics from different files and displays them together in a single file.

import subprocess,sys,re
OutFileName=sys.argv[1]+"/"+sys.argv[1]+".StatSummary"  ## Specifies the name of the Output File
OutFileHandle=open(OutFileName,'w')                                 ## Opens Output File for writing
OutFileHandle.write("SampleName" + "\t" + "TotalStartingReads"+"\t"+"TotalReadsInputForAlignment" + "\t" + "PercentageLost"+ "\t" + "ReadsMapped" +"\t" + "PercentageReadMapping"+ "\t" +"TotalOnTargetReads" + "\t"+ "PercentageOnTarget" + "\t"+ "TotalOnTargetFilterReads"+ "\t"+ "PercentageOnTargetAfterFilter"+ "\t"+ "PercentageUsable"+"\t"+"MeanCoverage" + "\t" + "PercentBases_above_0" + "\t" + "PercentBases_above_1" + "\t" + "PercentBases_above_250" + "\t" + "PercentBases_above_1000" + "\t"+ "CovBelow250XAmpliconNum"+ "\t"+ "ClipCount" +"\n") ## Writes the header to the Output File.
FileName=sys.argv[1]+"/"+sys.argv[1]+".novo.sorted.flagstat"              ## Specifies the name of sample's initial-sam flagstat file.
FileName2=sys.argv[1]+"/"+sys.argv[1]+".novo.sorted.fixed.ontarget.flagstat"    ## Specifies the name of sample's ontarget-sam flagstat file.
FileName3=sys.argv[1]+"/"+sys.argv[1]+".novo.sorted.fixed.ontarget.afterFilter.flagstat"   ## Specifies the name of sample's ontarget.filtered-sam flagstat file.
FileName4=sys.argv[1]+"/"+sys.argv[1]+".Depth.sample_summary"      ## Specifies the name of sample's Depth Summary file.
AmpliconFile=sys.argv[1]+"/"+sys.argv[1]+".Depth.below250X.amplicons"
ClipStatFile=sys.argv[1]+"/"+sys.argv[1]+".clippingstats.txt"
BtrimStatsR1=sys.argv[1]+"/"+sys.argv[1]+"_R1.btrim.stats.txt"
BtrimStatsR2=sys.argv[1]+"/"+sys.argv[1]+"_R2.btrim.stats.txt"




BtrimStatsR1Handle=open(BtrimStatsR1,'r')       ## Btrim Statistics for both end reads gets extracted.
BtrimStatsR2Handle=open(BtrimStatsR2,'r')
while 1:
	BtrimStatsR1Line=BtrimStatsR1Handle.readline()
	BtrimStatsR2Line=BtrimStatsR2Handle.readline()
	if(re.search("Total sequences",BtrimStatsR1Line) and re.search("Total sequences",BtrimStatsR2Line)):
		NumR1FastqReads=BtrimStatsR1Line.split(" ")[2]
		NumR2FastqReads=BtrimStatsR2Line.split(" ")[2]
		TotalStartingReads=int(NumR1FastqReads)+int(NumR2FastqReads)
		OutFileHandle.write(sys.argv[1] + "\t"+str(TotalStartingReads)+"\t")
		break

FileHandle=open(FileName,'r')                 ## Opens the initial-sam flagstat file and runs for loop on its lines.
FileLines=FileHandle.readlines()
for line in FileLines:
	line=line.rstrip()
	if(re.search('QC-passed reads',line)):   ## Detects the line containing the term "QC-passed reads" and processes it further to obtain the information on total number of reads present in the initial-sam(alignment) file.
		lineSections=line.split('+')     
		TotalReads=lineSections[0]
		TotalReads=TotalReads.rstrip()
	if(re.search('(?<!mate )mapped',line)):  ## Detects the line containing the term "mate mapped" and processes it further to obtain the information on total number of reads mapped to the reference genome.
		lineParts=line.split('+')
		ReadsMapped=lineParts[0]
		ReadsMapped=ReadsMapped.rstrip()
PercentageReadMapping="{0:.2f}".format(100 * float(ReadsMapped)/float(TotalReads))  ## Calculating and printing the percentage of reads mapped.
PercentageLost="{0:.2f}".format(100 * ( (float(TotalStartingReads)-float(TotalReads))/float(TotalStartingReads)))
OutFileHandle.write(TotalReads + "\t"+ PercentageLost + "\t" + ReadsMapped + "\t" + PercentageReadMapping + "\t")



FileHandle2=open(FileName2,'r')                ## Opens up the sample's ontarget-sam flagstat file and processes it to obtain the total number of reads-on-target. Percentage of mapped reads that are on target is also calculated.
FileLines2=FileHandle2.readlines()
for line2 in FileLines2:
	line2=line2.rstrip()
	if(re.search('QC-passed reads',line2)):
		lineSections2=line2.split('+')
		TotalOnTargetReads=lineSections2[0]
		TotalOnTargetReads=TotalOnTargetReads.rstrip()
PercentageOnTarget="{0:.2f}".format(100 * float(TotalOnTargetReads)/float(ReadsMapped))
OutFileHandle.write(TotalOnTargetReads + "\t"+ PercentageOnTarget + "\t")


FileHandle3=open(FileName3,'r')               ## Opens up the sample's ontarget.filtered-sam flagstat file and processes it to obtain the percentage of mapped reads on target in the filtered sam file. Percentage of usable reads (i.e reads on target after filter) out of Total reads is also obtained. 
FileLines3=FileHandle3.readlines()
for line3 in FileLines3:
	line3=line3.rstrip()
        if(re.search('QC-passed reads',line3)):
		lineSections3=line3.split('+')
                TotalOnTargetReadsAfterFilter=lineSections3[0]
               	TotalOnTargetReadsAfterFilter=TotalOnTargetReadsAfterFilter.rstrip()
PercentageOnTargetAfterFilter = "{0:.2f}".format(100 * float(TotalOnTargetReadsAfterFilter)/float(ReadsMapped))
PercentageUsable="{0:.2f}".format(100 * float(TotalOnTargetReadsAfterFilter)/float(TotalStartingReads))
OutFileHandle.write(TotalOnTargetReadsAfterFilter + "\t" + PercentageOnTargetAfterFilter + "\t" + PercentageUsable + "\t")



FileHandle4=open(FileName4,'r')              ## Sample's Depth Summary file is opened and processed further to obtain MeanCoverage, PercentBases_above_0, PercentBases_above_1, PercentBases_above_250 and PercentBases_above_1000
FileLines4=FileHandle4.readlines()
for line4 in FileLines4:
	line4=line4.rstrip()
	if(re.search(sys.argv[1],line4)):
		lineSections4=line4.split('\t')
		MeanCoverage=lineSections4[2]
		PercentBases_above_0=lineSections4[6]
		PercentBases_above_1=lineSections4[7]
		PercentBases_above_250=lineSections4[8]
		PercentBases_above_1000=lineSections4[9]
OutFileHandle.write(MeanCoverage + "\t" + PercentBases_above_0 + "\t" + PercentBases_above_1 + "\t" + PercentBases_above_250 + "\t" + PercentBases_above_1000 +"\t")

AmpliconFileHandle=open(AmpliconFile,'r')     ## The total number of amplicons having reads with depth below 250X is counted here.
AmpliconFileAllLines=AmpliconFileHandle.readlines()
AmpliconCount=len(AmpliconFileAllLines)
OutFileHandle.write(str(AmpliconCount) +"\t")



ClipStatFileHandle=open(ClipStatFile,'r')    ## Extracting total number of quality-score clipped bases from the Sample's Clipping Stats file.
ClipStatLines=ClipStatFileHandle.readlines()
for eachline in ClipStatLines:
	eachline=eachline.rstrip()
	if(re.match('Number of quality-score clipped bases',eachline)):
		eachlineComponents=eachline.split(" ")
		clipCount=eachlineComponents[5]
		break
OutFileHandle.write(clipCount+"\n")

		


OutFileHandle.close()
FileHandle.close()
FileHandle2.close()
FileHandle3.close()
FileHandle4.close()
AmpliconFileHandle.close()
ClipStatFileHandle.close()
BtrimStatsR1Handle.close()
BtrimStatsR2Handle.close()


        
				
				
		


