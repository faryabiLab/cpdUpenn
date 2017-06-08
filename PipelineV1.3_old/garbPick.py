## Run this script as python garbPick.py <primer_R1.primParts File> <primer_R2.primParts File> <Sample.flagJoin File> <Sample.garbageStats> <Sample.original.novo.sam.flagValues> <primer_manifest File>


import re
import sys
import subprocess

primerR1file=open(sys.argv[1],'r')
primerR2file=open(sys.argv[2],'r')

flagJoinFileName=sys.argv[3]
##flagJoinFile=open(flagJoinFileName,'r')

garbageStatsFile=open(sys.argv[4],'w')
flagValueFile=open(sys.argv[5],'r')

primManifest=open(sys.argv[6],'r')

outFileNameList=[]



while 1:
	primerR1line=primerR1file.readline()
	primerR2line=primerR2file.readline()
	primManifestLine=primManifest.readline()
	if(primerR1line and primerR2line and primManifestLine):
		primManifestLine=primManifestLine.rstrip()
		primManifestCompo=primManifestLine.split("\t")
		outFileName=primManifestCompo[0]
		chromReg=primManifestCompo[1]
		outFileNameList.append(outFileName)
		garbageStatsFile.write("\t"+chromReg)
		
		outFile=open(outFileName,'w')
		primerR1line=primerR1line.rstrip()
		primerR2line=primerR2line.rstrip()
		primerR1=primerR1line.replace("\t","|")
		primerR2=primerR2line.replace("\t","|")
		primer=primerR1+"|"+primerR2
		outFile.write(subprocess.Popen(['grep','-P',primer,flagJoinFileName],stdout=subprocess.PIPE).communicate()[0])
		outFile.close()
	else:
		garbageStatsFile.write("\n")
		break

primerR1file.close()
primerR2file.close()
primManifest.close()




while 1:
	flagValue=flagValueFile.readline()
	if(flagValue):
		flagValue=flagValue.rstrip()
		garbageStatsFile.write(flagValue)
		for i in outFileNameList:
			p1=subprocess.Popen(['awk','{print $4}',i],stdout=subprocess.PIPE)
			p2=subprocess.Popen(['grep',flagValue],stdin=p1.stdout,stdout=subprocess.PIPE)
			p1.stdout.close()
			p3=subprocess.Popen(['wc','-l'],stdin=p2.stdout,stdout=subprocess.PIPE).communicate()[0]
			p2.stdout.close()
			samFlagOccurence=p3.rstrip()
			garbageStatsFile.write("\t"+samFlagOccurence)
		garbageStatsFile.write("\n")
	else:
		break	
		

garbageStatsFile.close()
flagValueFile.close()
