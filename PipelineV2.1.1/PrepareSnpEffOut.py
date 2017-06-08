##V2.1.1
## Run this script as: python PrepareSnpEffOut.py <Sample.combined.HRunTagged.canon_snpeff(or snpeff).vT.ex.txt> <Sample.combined.HRunTagged.canon_snpeff(or snpEff).vT.ex.tab>
## Changes in V2.1.1 : 1) Changed 'POS=POS+1' to 'POS=POS+len(Ref)' for insertions to make the position annotation for them more accurate.

import re,string,sys
InputFile=sys.argv[1]
OutputFile=sys.argv[2]

InputFileHandle=open(InputFile,'r')
OutputFileHandle=open(OutputFile,'w')

def PrepFunction(Change,ANN,ChangeType,CHROM,POS,Ref,Quality,rsID,Hom_or_Het,Coverage):
        Warnings=ANN.split("|")[15]
        Gene_ID=ANN.split("|")[4]
        Gene_name=ANN.split("|")[3]
        BioType=ANN.split("|")[7]
	try:
		c_Change=ANN.split("|")[9]
	except IndexError, IndErr:
		c_Change=''
	try:
		p_Change=ANN.split("|")[10]
	except IndexError, IndErr:
		p_Change=''
        Transcript_ID=ANN.split("|")[6]
	try:
        	Exon_Rank=ANN.split("|")[8].split("/")[0]
	except IndexError, IndErr:
		Exon_Rank=''
        Annotation=ANN.split("|")[1]
	try:
        	AA_pos=ANN.split("|")[13].split("/")[0]
	except IndexError, IndErr:
		AA_pos=''
	try:
        	CDS_length=ANN.split("|")[12].split("/")[1]
	except IndexError, IndErr:
		CDS_length=''
	Exon_ID=""
        Old_new_AA=""
        Old_new_Codon=""
        Codon_Degeneracy=""
        Codons_around=""
        AAs_around=""
        if(ChangeType=="INS"):
                POS=int(POS)+len(Ref)
		Diff=len(Ref)-len(Change)
		Change="+"+Change[Diff:]
                Ref="*"
        elif(ChangeType=="DEL"):
		Diff=len(Change)-len(Ref)
                Change="-"+Ref[Diff:]
		POS=int(POS)+len(Ref)+Diff
		Ref="*"
	linestring=CHROM+"\t"+str(POS)+"\t"+Ref+"\t"+Change+"\t"+ChangeType+"\t"+Hom_or_Het+"\t"+Quality+"\t"+Coverage+"\t"+Warnings+"\t"+Gene_ID+"\t"+Gene_name+"\t"+BioType+"\t"+Transcript_ID+"\t"+Exon_ID+"\t"+Exon_Rank+"\t"+Annotation+"\t"+Old_new_AA+"\t"+Old_new_Codon+"\t"+AA_pos+"\t"+Codon_Degeneracy+"\t"+CDS_length+"\t"+Codons_around+"\t"+AAs_around+"\t"+rsID+"\t"+c_Change+"\t"+p_Change+"\n"
        OutputFileHandle.write(linestring)

while 1:
	line=InputFileHandle.readline()
	line=line.rstrip()
	if(line):
		if(re.match('^#', line)):
			OutputFileHandle.write("#Chromo"+"\t"+"Position"+"\t"+"Reference"+"\t"+"Change"+"\t"+"Change_type"+"\t"+"Homozygous"+"\t"+"Quality"+"\t"+"Coverage"+"\t"+"Warnings"+"\t"+"Gene_ID"+"\t"+"Gene_name"+"\t"+"Bio_type"+"\t"+"Trancript_ID"+"\t"+"Exon_ID"+"\t"+"Exon_Rank"+"\t"+"Effect"+"\t"+"old_AA/new_AA"+"\t"+"Old_codon/New_codon"+"\t"+"Codon_Num(CDS)"+"\t"+"Codon_Degeneracy"+"\t"+"CDS_size"+"\t"+"Codons_around"+"\t"+"AAs_around"+"\t"+"Custom_interval_ID"+"\t"+"c_Change"+"\t"+"p_Change"+"\n")
		else:
			lineCompo=line.split("\t")
			CHROM,POS,Ref,Change,Quality,rsID=lineCompo[0],lineCompo[1],lineCompo[3],lineCompo[4],lineCompo[5],lineCompo[2]
			Hom,Het=lineCompo[9],lineCompo[10]
			if(Hom=="true"):
				Hom_or_Het="Hom"
			elif(Het=="true"):
				Hom_or_Het="Het"
			Coverage=lineCompo[7]
			if(string.find(Change,",")!=-1):
				ChangeComp=Change.split(",")
				ANNComp=lineCompo[8].split(",")
				ChangeTypeComp=lineCompo[11].split(",")
				for i in range(len(ChangeComp)):
					PrepFunction(ChangeComp[i],ANNComp[i],ChangeTypeComp[i],CHROM,POS,Ref,Quality,rsID,Hom_or_Het,Coverage)
			else:
				ANN=lineCompo[8]
				ChangeType=lineCompo[11]
				if(string.find(ANN,",")!=-1):
					ANNComp=ANN.split(",")
					for i in range(len(ANNComp)):
						PrepFunction(Change,ANNComp[i],ChangeType,CHROM,POS,Ref,Quality,rsID,Hom_or_Het,Coverage)
				else:
					PrepFunction(Change,ANN,ChangeType,CHROM,POS,Ref,Quality,rsID,Hom_or_Het,Coverage)
	else:
		break
InputFileHandle.close()
OutputFileHandle.close()
