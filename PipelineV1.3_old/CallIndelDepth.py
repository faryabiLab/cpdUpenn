##V1.3
##Run this script as : python CallIndelDepth.py <Sample.indels.BED> <Sample.novo.sorted.fixed.bam.mpileup> <Sample.novo.sorted.fixed.ontarget.afterFilter.clipped.bam.mpileup> <Sample.IndelDepthForAnnotation.vcf>
##This script takes in the indels.bed file prepared using the cannon_snpeff file. 
##The indel information in the bed file (i.e. position, indel-type, indel-length and indel-string ) is passed to a function which calculates the indel-depth, reference-depth
##and the total depth at indel-position. These depth values which are returned by the function are then printed in an output vcf file named IndelDepthForAnnotaion.vcf alongwith the other necessary indel information.

import sys
import re
import string

IndelBed=open(sys.argv[1],'r')
OutFile=open(sys.argv[4],'w')
OutFile.write("##fileformat=VCFv4.1"+"\n")
OutFile.write("#CHROM"+"\t"+"POS"+"\t"+"ID"+"\t"+"REF"+"\t"+"ALT"+"\t"+"QUAL"+"\t"+"FILTER"+"\t"+"INFO"+"\t"+"FORMAT"+"\t"+"SAMPLE"+"\n")




def checkInBam(position,inType,inLen,inString):
	pre_position=str(int(position)-1)
	
        if(int(inLen)>6):
                filename=open(sys.argv[2],'r')
        else:
                filename=open(sys.argv[3],'r')
	while 1:
		fileLine=filename.readline()
		if(fileLine):
			fileLine=fileLine.rstrip()
			fileLineCompo=fileLine.split()[0:5]
			if(inType=='INS'):
				if(fileLineCompo[1]==pre_position):
					TotalDepth=int(fileLineCompo[3])
					inString=string.replace(inString,"+","+"+inLen)
					ReadString=fileLineCompo[4].upper()
					IndelDepth=ReadString.count(inString)
					RefDepth=ReadString.count('.')+ReadString.count(',')
					break
			elif(inType=='DEL'):
				if(fileLineCompo[1]==pre_position):
					inString=string.replace(inString,"-","-"+inLen)
					ReadString=fileLineCompo[4].upper()
					IndelDepth=ReadString.count(inString)
					continue
				elif(fileLineCompo[1]==position):
					ReadString=fileLineCompo[4]
					RefDepth=ReadString.count('.')+ReadString.count(',')
					TotalDepth=int(fileLineCompo[3])
					break
		else:
			break
	
	filename.close()
        return(TotalDepth,RefDepth,IndelDepth)




while 1:
	IndelBedLine=IndelBed.readline()
	if(IndelBedLine):
		IndelBedLine=IndelBedLine.rstrip()
		if(re.match('#',IndelBedLine)):
			continue
		else:
			IndelBedLineCompo=IndelBedLine.split("\t")
			IndelPos=IndelBedLineCompo[2]
			IndelType=IndelBedLineCompo[5]
			IndelLength=IndelBedLineCompo[6]
			IndelString=IndelBedLineCompo[4]
			FDP,FRD,FAD=checkInBam(IndelPos,IndelType,IndelLength,IndelString)
			FAF=100 * float(FAD) / float(FDP)
			OutFile.write(IndelBedLineCompo[0]+"\t"+IndelBedLineCompo[2]+"\t"+"."+"\t"+"\t".join(IndelBedLineCompo[3:5])+"\t"+"."+"\t"+"."+"\t"+"FDP="+str(FDP)+";FRD="+str(FRD)+";FAD="+str(FAD)+";FAF="+"{0:.2f}".format(FAF)+"\t"+"."+"\t"+"."+"\n")
	else:
		break


OutFile.close()





	
	
			
		
		
		

			
