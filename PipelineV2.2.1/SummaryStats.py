##V2.1.1_hpc
## Run this script as: python SummaryStats.py <SampleFolderPath> <Library> <SampleCount>
## This script collects the statistics from different files and displays them together in a single file.
## Changes in V1.3 : 1) Division-by-Zero errors are now handled appropriately.
##                   2) No more usage of Btrim Stats to pull out total no. of starting reads. The script now uses the fastqs for that count.
##                   3) "Library" and "SampleCount" were added as command-line arguments.
## Changes in V1.3.3:1) Now R1_fastq and R2_fastq variables are set to the respective cutAd.fastq filenames for TSCA samples. Therefore, the if condition related to defining R1_fastq and R2_fastq has been removed from the script below.
##                   2) Replaced integer 0 with string "0" in the if statements 'if(TotalReads!=0)' and 'if(ReadsMapped!=0)'.
## Changes in V2.0:  1) Amplicon2File variable added to denote the file name of 'Depth.below150X.amplicons' file. os.path.exists(Amplicon2File) statement added in the main "if" statement of the script. The logic to count and print the total no. of lines in Depth.below150X.amplicons file also added. Hence, the output,i.e. Sample.SummaryStats file would now have an extra column "No._of_amplicons(exonic_part)_below_150x" at the end of the file.
##                   2) The "CovBelow250XAmpliconNum" name in header is changed to "No._of_amplicons(exonic_part)_below_250x".
## Changes in V2.1.1_hpc :1) Replaced the argument <SampleName> with <SampleFolderPath> on the command-line. Added variables SampleFolderPath and Sample_name which were assigned the values of <SampleFolderPath> and "Sample_name extracted from <SampleFolderPath>" respectively. Made the changes accordingly downstream.
##			  2) Changed the value of PercentageReadMapping from 0 to "0" for the conditions when "TotalReads is NOT not-equal to "0"" and when "TotalReads!="0" but ReadsMapped="0"".

import subprocess,sys,re,os
SampleFolderPath=sys.argv[1]
Sample_name=SampleFolderPath.split("/")[5]

FileName=SampleFolderPath+"/"+Sample_name+".novo.sorted.flagstat"              ## Specifies the name of sample's initial-sam flagstat file.
FileName2=SampleFolderPath+"/"+Sample_name+".novo.sorted.fixed.ontarget.flagstat"    ## Specifies the name of sample's ontarget-sam flagstat file.
FileName3=SampleFolderPath+"/"+Sample_name+".novo.sorted.fixed.ontarget.afterFilter.flagstat"   ## Specifies the name of sample's ontarget.filtered-sam flagstat file.
FileName4=SampleFolderPath+"/"+Sample_name+".Depth.sample_summary"      ## Specifies the name of sample's Depth Summary file.
AmpliconFile=SampleFolderPath+"/"+Sample_name+".Depth.below250X.amplicons"
Amplicon2File=SampleFolderPath+"/"+Sample_name+".Depth.below150X.amplicons"
ClipStatFile=SampleFolderPath+"/"+Sample_name+".clippingstats.txt"
Sample_lib=sys.argv[2]
Sample_count=sys.argv[3]

R1_fastq=SampleFolderPath+"/"+Sample_name+"_"+Sample_count+"_L001_R1_001.cutAd.fastq"
R2_fastq=SampleFolderPath+"/"+Sample_name+"_"+Sample_count+"_L001_R2_001.cutAd.fastq"

if(os.path.exists(FileName) and os.path.exists(FileName2) and os.path.exists(FileName3) and os.path.exists(FileName4) and os.path.exists(AmpliconFile) and os.path.exists(Amplicon2File) and os.path.exists(ClipStatFile) and os.path.exists(R1_fastq) and os.path.exists(R2_fastq)):
	OutFileName=SampleFolderPath+"/"+Sample_name+".StatSummary"  ## Specifies the name of the Output File
	OutFileHandle=open(OutFileName,'w')                                 ## Opens Output File for writing
	OutFileHandle.write("SampleName" + "\t" + "TotalStartingReads"+"\t"+"TotalReadsInputForAlignment" + "\t" + "PercentageLost"+ "\t" + "ReadsMapped" +"\t" + "PercentageReadMapping"+ "\t" +"TotalOnTargetReads" + "\t"+ "PercentageOnTarget" + "\t"+ "TotalOnTargetFilterReads"+ "\t"+ "PercentageOnTargetAfterFilter"+ "\t"+ "PercentageUsable"+"\t"+"MeanCoverage" + "\t" + "PercentBases_above_0" + "\t" + "PercentBases_above_1" + "\t" + "PercentBases_above_250" + "\t" + "PercentBases_above_1000" + "\t"+ "No._of_amplicons(exonic_part)_below_250x"+ "\t"+ "ClipCount"+"\t"+"No._of_amplicons(exonic_part)_below_150x"+"\n") ## Writes the header to the Output File.
	
	## Calculating total no. of reads from R1 fastq.
	SubProcessR1=subprocess.Popen(['wc','-l',R1_fastq],stdout=subprocess.PIPE).communicate()[0].split(" ")[0]
	NumR1FastqReads=int(SubProcessR1)/4
	
	## Calculating total no. of reads from R2 fastq.
	SubProcessR2=subprocess.Popen(['wc','-l',R2_fastq],stdout=subprocess.PIPE).communicate()[0].split(" ")[0]
	NumR2FastqReads=int(SubProcessR2)/4
	
	TotalStartingReads=NumR1FastqReads+NumR2FastqReads
	
	if(TotalStartingReads!=0):
		OutFileHandle.write(Sample_name + "\t"+str(TotalStartingReads)+"\t")
	else:
		OutString=("\t"+"0")*17      ## Assigning value "0" to all fields in SummaryStats file EXCEPT the SampleName field
		OutFileHandle.write(Sample_name+OutString+"\n")
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
	
	if(TotalReads!="0"):
		if(ReadsMapped!="0"):
			PercentageReadMapping="{0:.2f}".format(100 * float(ReadsMapped)/float(TotalReads))  ## Calculating and printing the percentage of reads mapped.
			OutFileHandle.write(TotalReads + "\t"+ PercentageLost + "\t" + ReadsMapped + "\t" + PercentageReadMapping + "\t")
		else:
			PercentageReadMapping="0"
			OutString1=("\t"+"0")*12
			OutFileHandle.write(TotalReads + "\t"+ PercentageLost + "\t" + ReadsMapped + "\t" + PercentageReadMapping + OutString1 + "\n")
			OutFileHandle.close()
			sys.exit("Reads Mapped are zero. Made StatSummary File accordingly\n")
	else:
		PercentageReadMapping="0"
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
		if(re.search(Sample_name,line4)):
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
	OutFileHandle.write(clipCount+"\t")
	
	Amplicon2FileHandle=open(Amplicon2File,'r')     ## The total number of amplicons having reads with depth below 150X is counted here.
	Amplicon2FileAllLines=Amplicon2FileHandle.readlines()
	Amplicon2Count=len(Amplicon2FileAllLines)
	OutFileHandle.write(str(Amplicon2Count) +"\n")
	
	OutFileHandle.close()
	FileHandle.close()
	FileHandle2.close()
	FileHandle3.close()
	FileHandle4.close()
	AmpliconFileHandle.close()
	ClipStatFileHandle.close()
	Amplicon2FileHandle.close()
	
else:
	sys.exit("Error: Not enough files to generate SummaryStats File\n")

        
				
				
		


