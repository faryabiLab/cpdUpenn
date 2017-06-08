##V1.3
## Run this script as: python SummaryStats.py <SampleName> <Library> <SampleCount>
## This script collects the statistics from different files and displays them together in a single file.
## Changes in V1.3 : 1) Division-by-Zero errors are now handled appropriately.
##                   2) No more usage of Btrim Stats to pull out total no. of starting reads. The script now uses the fastqs for that count.
##                   3) "Library" and "SampleCount" were added as command-line arguments.

import subprocess,sys,re,os

FileName=sys.argv[1]+"/"+sys.argv[1]+".novo.sorted.flagstat"              ## Specifies the name of sample's initial-sam flagstat file.
FileName2=sys.argv[1]+"/"+sys.argv[1]+".novo.sorted.fixed.ontarget.flagstat"    ## Specifies the name of sample's ontarget-sam flagstat file.
FileName3=sys.argv[1]+"/"+sys.argv[1]+".novo.sorted.fixed.ontarget.afterFilter.flagstat"   ## Specifies the name of sample's ontarget.filtered-sam flagstat file.
FileName4=sys.argv[1]+"/"+sys.argv[1]+".Depth.sample_summary"      ## Specifies the name of sample's Depth Summary file.
AmpliconFile=sys.argv[1]+"/"+sys.argv[1]+".Depth.below250X.amplicons"
ClipStatFile=sys.argv[1]+"/"+sys.argv[1]+".clippingstats.txt"
Sample_lib=sys.argv[2]
Sample_count=sys.argv[3]


if(Sample_lib=="TSCA"):
	R1_fastq=sys.argv[1]+"/"+sys.argv[1]+"_"+Sample_count+"_L001_R1_001.fastq"
	R2_fastq=sys.argv[1]+"/"+sys.argv[1]+"_"+Sample_count+"_L001_R2_001.fastq"
elif(Sample_lib=="PPP" or Sample_lib=="CEBPa"):
	R1_fastq=sys.argv[1]+"/"+sys.argv[1]+"_"+Sample_count+"_L001_R1_001.cutAd.fastq"
	R2_fastq=sys.argv[1]+"/"+sys.argv[1]+"_"+Sample_count+"_L001_R2_001.cutAd.fastq"

if(os.path.exists(FileName) and os.path.exists(FileName2) and os.path.exists(FileName3) and os.path.exists(FileName4) and os.path.exists(AmpliconFile) and os.path.exists(ClipStatFile) and os.path.exists(R1_fastq) and os.path.exists(R2_fastq)):
	OutFileName=sys.argv[1]+"/"+sys.argv[1]+".StatSummary"  ## Specifies the name of the Output File
	OutFileHandle=open(OutFileName,'w')                                 ## Opens Output File for writing
	OutFileHandle.write("SampleName" + "\t" + "TotalStartingReads"+"\t"+"TotalReadsInputForAlignment" + "\t" + "PercentageLost"+ "\t" + "ReadsMapped" +"\t" + "PercentageReadMapping"+ "\t" +"TotalOnTargetReads" + "\t"+ "PercentageOnTarget" + "\t"+ "TotalOnTargetFilterReads"+ "\t"+ "PercentageOnTargetAfterFilter"+ "\t"+ "PercentageUsable"+"\t"+"MeanCoverage" + "\t" + "PercentBases_above_0" + "\t" + "PercentBases_above_1" + "\t" + "PercentBases_above_250" + "\t" + "PercentBases_above_1000" + "\t"+ "CovBelow250XAmpliconNum"+ "\t"+ "ClipCount"+"\n") ## Writes the header to the Output File.
	
	## Calculating total no. of reads from R1 fastq.
	SubProcessR1=subprocess.Popen(['wc','-l',R1_fastq],stdout=subprocess.PIPE).communicate()[0].split(" ")[0]
	NumR1FastqReads=int(SubProcessR1)/4
	
	## Calculating total no. of reads from R2 fastq.
	SubProcessR2=subprocess.Popen(['wc','-l',R2_fastq],stdout=subprocess.PIPE).communicate()[0].split(" ")[0]
	NumR2FastqReads=int(SubProcessR2)/4
	
	TotalStartingReads=NumR1FastqReads+NumR2FastqReads
	
	if(TotalStartingReads!=0):
		OutFileHandle.write(sys.argv[1] + "\t"+str(TotalStartingReads)+"\t")
	else:
		OutString=("\t"+"0")*17      ## Assigning value "0" to all fields in SummaryStats file EXCEPT the SampleName field
		OutFileHandle.write(sys.argv[1]+OutString+"\n")
		OutFileHandle.close()
		sys.exit("No reads in original fastq. Made StatSummary File accordingly\n")
	
	
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
	
	PercentageLost="{0:.2f}".format(100 * ( (float(TotalStartingReads)-float(TotalReads))/float(TotalStartingReads)))
	
	if(TotalReads!=0):
		if(ReadsMapped!=0):
			PercentageReadMapping="{0:.2f}".format(100 * float(ReadsMapped)/float(TotalReads))  ## Calculating and printing the percentage of reads mapped.
			OutFileHandle.write(TotalReads + "\t"+ PercentageLost + "\t" + ReadsMapped + "\t" + PercentageReadMapping + "\t")
		else:
			PercentageReadMapping=0
			OutString1=("\t"+"0")*12
			OutFileHandle.write(TotalReads + "\t"+ PercentageLost + "\t" + ReadsMapped + "\t" + PercentageReadMapping + OutString1 + "\n")
			OutFileHandle.close()
			sys.exit("Reads Mapped are zero. Made StatSummary File accordingly\n")
	else:
		PercentageReadMapping=0
		OutString2=("\t"+"0")*12
		OutFileHandle.write(TotalReads + "\t"+ PercentageLost + "\t" + ReadsMapped + "\t" + PercentageReadMapping + OutString2 + "\n")
		OutFileHandle.close()
		sys.exit("No Reads Remaining After Btrim Step. Made StatSummary File accordingly\n")
	
	
	
	
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
	
else:
	sys.exit("Error: Not enough files to generate SummaryStats File\n")

        
				
				
		


