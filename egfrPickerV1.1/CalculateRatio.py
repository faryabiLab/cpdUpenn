##V1.1
##Run this script as python CalculateRatio.py <RunSummaryFinal_file> <RunSummaryFinalCalcuations_file>
##Changes in V1.1: 1) Removed HPRT1,GUCB and ACTB amplicons from the calculations.
import sys
import subprocess

summary=sys.argv[1]
summaryHandle=open(summary,'r')
outFile=sys.argv[2]
outFileHandle=open(outFile,'w')
while 1:
	summaryLine=summaryHandle.readline()
	if(summaryLine):
		summaryLine=summaryLine.rstrip()
		summaryComp=summaryLine.split("\t") ##Extracting Components of summaryLine
		if(summaryComp[0]=="SampleName"):
			summaryComp.append("FDP")
			summaryComp.append("PercentageUsable")
			summaryComp.append("EGFRvIII_ratio")
			summaryComp.append("EGFR_expression_level")
			summaryComp.append("wt_EGFR_expression_level")
			summaryComp.append("EGFRvIII_expression_level")
			summaryComp.append("Degradation_9_12/9_10")
			summaryComp.append("Degradation_9_9/9_10")
			for i in summaryComp:
				outFileHandle.write(i+"\t")
			outFileHandle.write("\n")
		else:
			sampleName,EGFR_exon1_intron1,EGFR_exon1_exon2,EGFR_exon1_exon8,HPRT4,SDHA,RPL13A,EGFR_exon9_exon12,EGFR_exon9_exon10,EGFR_exon9_exon9=summaryComp[0],int(summaryComp[1]),int(summaryComp[2]),int(summaryComp[3]),int(summaryComp[4]),int(summaryComp[5]),int(summaryComp[6]),int(summaryComp[7]),int(summaryComp[8]),int(summaryComp[9])
			#Calculating FDP (i.e. filtered depth)
			joinedFile=summaryComp[0]+"/"+summaryComp[0]+".btrim.summary.joined"
			fdp=int(subprocess.Popen(['wc','-l',joinedFile],stdout=subprocess.PIPE).communicate()[0].split(" ")[0])
			summaryComp.append(fdp)
			#Calculating percentUsable,i.e. observed_filtered_depth(OFPD) / total_depth(FDP)
			ofdp=EGFR_exon1_intron1+EGFR_exon1_exon2+EGFR_exon1_exon8+HPRT4+SDHA+RPL13A+EGFR_exon9_exon12+EGFR_exon9_exon10+EGFR_exon9_exon9
			if(fdp!=0):
				percentUsable="{0:.2f}".format((float(ofdp)/float(fdp))*100)
			else:
				percentUsable="N/A"
			summaryComp.append(percentUsable)
			#Calculating other ratios
			EGFR_wt_plus_vIII=EGFR_exon1_exon2+EGFR_exon1_exon8
			HouseKeepingControl=HPRT4+SDHA+RPL13A
			
			if(EGFR_wt_plus_vIII!=0):
				EGFRvIII_ratio="{0:.2f}".format(float(EGFR_exon1_exon8)/float(EGFR_wt_plus_vIII))
			else:
				EGFRvIII_ratio="N/A"
			summaryComp.append(EGFRvIII_ratio)
			
			if(HouseKeepingControl!=0):
				EGFR_exp_level="{0:.2f}".format(float(EGFR_wt_plus_vIII)/float(HouseKeepingControl))
				wt_EGFR_exp_level="{0:.2f}".format(float(EGFR_exon1_exon2)/float(HouseKeepingControl))
				EGFRvIII_exp_level="{0:.2f}".format(float(EGFR_exon1_exon8)/float(HouseKeepingControl))
			else:
				EGFR_exp_level="N/A"
				wt_EGFR_exp_level="N/A"
				EGFRvIII_exp_level="N/A"
			summaryComp.append(EGFR_exp_level)
			summaryComp.append(wt_EGFR_exp_level)
			summaryComp.append(EGFRvIII_exp_level)
			
			if(EGFR_exon9_exon10!=0):
				Degradation_9_12_9_10="{0:.2f}".format(float(EGFR_exon9_exon12)/float(EGFR_exon9_exon10))
				Degradation_9_9_9_10="{0:.2f}".format(float(EGFR_exon9_exon9)/float(EGFR_exon9_exon10))
			else:
				Degradation_9_12_9_10="N/A"
				Degradation_9_9_9_10="N/A"
			summaryComp.append(Degradation_9_12_9_10)
			summaryComp.append(Degradation_9_9_9_10)
				
			#Printing all the Components to the outFile
			
			for i in summaryComp:
                                outFileHandle.write(str(i)+"\t")
                        outFileHandle.write("\n")
	else:
		break

outFileHandle.close()
summaryHandle.close()

