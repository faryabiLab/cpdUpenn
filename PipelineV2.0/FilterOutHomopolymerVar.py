##V2.0
##Run this script as: python FilterOutHomopolymer.py <Sample.canon_annotated_flagged> <Sample.canon_annotated_flagged.HRR>
import sys
import re
Input_File=sys.argv[1]
Output_File=sys.argv[2]

Input_File_Handle=open(Input_File,'r')
Output_File_Handle=open(Output_File,'w')

while 1:
	InLine=Input_File_Handle.readline()
	if(InLine):
		if(re.match('^SampleName',InLine)):
			Output_File_Handle.write(InLine)
		else:
			InLine=InLine.rstrip()
			VarType=InLine.split("\t")[5]
			Tag=InLine.split("\t")[39]
			if(VarType=="SNP" and Tag=="HomopolymerRun"):
				FAF=InLine.split("\t")[54]
				if(float(FAF)>=10.00):
					Output_File_Handle.write(InLine+"\n")
			elif((VarType=="INS" and Tag=="HomopolymerRun") or (VarType=="DEL" and Tag=="HomopolymerRun")):
				GT=InLine.split("\t")[42]
				DP=GT.split(":")[2]
				AD1=GT.split(":")[1].split(",")[1]
				AF1=100*float(AD1)/float(DP)
				try:
					AD2=GT.split(":")[1].split(",")[2]
				except IndexError, indErr:
					AD2=''
				try:
					AD3=GT.split(":")[1].split(",")[3]
				except IndexError, indErr:
					AD3=''
				if(AD2):
					AF2=100*float(AD2)/float(DP)
					if(AD3):
						AF3=100*float(AD3)/float(DP)
						if(AF1>=10.00 or AF2>=10.00 or AF3>=10.00):
							Output_File_Handle.write(InLine+"\n")
					else:
						if(AF1>=10.00 or AF2>=10.00):
							Output_File_Handle.write(InLine+"\n")
				else:
					if(AF1>=10.00):
						Output_File_Handle.write(InLine+"\n")
			else:
				Output_File_Handle.write(InLine+"\n")
	else:
		break;
Input_File_Handle.close()
Output_File_Handle.close()

