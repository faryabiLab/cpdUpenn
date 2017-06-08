##V1.2.4
## Run this script as: python AnnotationJoin.py <snpeff annotated variant text file> <refseq annotated variant file> <in-house 1000g database> <sample's DepthForAnnotation file> <in-house alalmut database> <Joined_Annotation_Output file>
## This script joins all the annotated files and annotation databases together giving a single Annotation Table for variants (SNVs and Indels).
## Starting from Version 1.2.4, the last column would have the header (name) "AdjacentSNPflag" instead of "dummy_value".


import sys,re, string
File1=open(sys.argv[1],'r')  ## opens up each annotation file and stores their lines in the respective lists.
File2=open(sys.argv[2],'r')
File3=open(sys.argv[3],'r')
File4=open(sys.argv[4],'r')
File5=open(sys.argv[5],'r')
OutFile=open(sys.argv[6],'w')
List1=File1.readlines()
List2=File2.readlines()
List3=File3.readlines()
List4=File4.readlines()
List5=File5.readlines()
SampleNameBreakdown=sys.argv[1].split("/")
## Write the header to the output file
OutFile.write("SampleName\tChrom\tPos\tRef\tAlt\tVariant Type(SnpEff)\tHomozygous\tQuality\tCoverage\tWarnings\tGene_ID\tGene_name\tBio_type\tTrancript_ID\tExon_ID\tExon_Rank\tEffect\told_AA/new_AA\tOld_codon/New_codon\tCodon_Num(CDS)\tCodon_Degeneracy\tCDS_size\tCodons_around\tAAs_around\tCustom_interval_ID\tLocation(annovar)\tGene\tConsequence\tTranscript\tcDNA change\tProtien Change\tCHROM(vcf)\tPOS\tdbsnp137_ID\tRef\tAlt\tQUAL\tFILTER\tINFO\tGenotypeFormat\tGenotype\tAC(1000g)\tAN(1000g)\tAF(1000g)\tAMR_AF\tASN_AF\tEUR_AF\tAFR_AF\tVariantType\tFDP\tFRD\tFAD\tFAF\tgene\tgeneId\tchrom\ttranscript\tstrand\ttransLen\tprotein\tUniprot\tvarType\tcodingEffect\tvarLocation\tassembly\tgDNAstart\tgDNAend\tgNomen\tcDNAstart\tcDNAend\tcNomen\tpNomen\talt_pNomen\texon\tintron\tomimId\tdistNearestSS\tnearestSSType\twtSSFScore\twtMaxEntScore\twtNNSScore\twtGSScore\twtHSFScore\tvarSSFScore\tvarMaxEntScore\tvarNNSScore\tvarGSScore\tvarHSFScore\tnearestSSChange\tlocalSpliceEffect\tproteinDomain1\tproteinDomain2\tproteinDomain3\tproteinDomain4\trsId\trsValidated\trsAncestralAllele\trsHeterozygosity\trsClinicalSignificance\trsMAF\trsMAFAllele\trsMAFCount\tinsNucs\tdelNucs\tsubstType\twtNuc\tvarNuc\tnucChange\tphastCons\tphyloP\twtAA_1\twtAA_3\twtCodon\twtCodonFreq\tvarAA_1\tvarAA_3\tvarCodon\tvarCodonFreq\tposAA\tnOrthos\tconservedOrthos\tconservedDistSpecies\tBLOSUM45\tBLOSUM62\tBLOSUM80\twtAAcomposition\tvarAAcomposition\twtAApolarity\tvarAApolarity\twtAAvolume\tvarAAvolume\tgranthamDist\tAGVGDclass\tAGVGDgv\tAGVGDgd\tSIFTprediction\tSIFTweight\tSIFTmedian\tMAPPprediction\tMAPPpValue\tMAPPpValueMedian\tAdjacentSNPflag"+"\n")

for snpEffLine in List1:  ## swifts through the snpeff file and finds all the entries in the other annotation files as well as databases with the similar "Chrom Pos Ref Alt Transcript"  OR "Chrom Pos Ref Alt" entries.
	refseqFlag=0
	AlalmutFlag=0
	VarFlag=0
	DepthFlag=0
	snpEffLine=snpEffLine.rstrip()
	if((string.find(snpEffLine,"#"))==-1): ## skips the snpeff header
		snpEffComponents=snpEffLine.split("\t")
		if(len(snpEffComponents)<24):  ## due to the empty fields at the end of snpeff line, the splitting on tab gives errors. So this sections of script manages those empty fields and puts them appropriately in the output file.
			emptyFields=24-len(snpEffComponents)
			snpEffAnno1="\t".join(snpEffComponents[0:len(snpEffComponents)])
			snpEffAnno2=("\t"+"")*emptyFields
			snpEffAnno=snpEffAnno1+snpEffAnno2
		else:
			snpEffAnno="\t".join(snpEffComponents)
		snpEffText="chr"+"\t".join(snpEffComponents[0:4])+"\t"+snpEffComponents[12] ## extracting snpeff text --> "Chrom Pos Ref Alt Transcript"
		snpEffTextWithoutTranscript="chr"+"\t".join(snpEffComponents[0:4]) ## extracting snpeff text --> "Chrom Pos Ref Alt"
		snpEffTextWithoutTranscriptWithGene="chr"+"\t".join(snpEffComponents[0:4])+"\t"+snpEffComponents[10] ## extracting snpeff text --> "Chrom Pos Ref Alt GeneName"
		snpEffTextForINS="chr"+"\t".join(snpEffComponents[0:2])+"\t"+snpEffComponents[3]+"\t"+snpEffComponents[12] ## extracting snpeff text for insertions --> "Chrom Pos Alt Transcript"
		snpEffTextForINSWithoutTranscriptWithGene="chr"+"\t".join(snpEffComponents[0:2])+"\t"+snpEffComponents[3] + "\t"+ snpEffComponents[10] ## extracting snpeff text for insertions --> "Chrom Pos Alt Gene"
	        snpEffTextForDEL="chr"+"\t".join(snpEffComponents[0:2])+"\t"+snpEffComponents[3]+"\t"+snpEffComponents[12] ## extracting snpeff text for deletions --> "Chrom Pos Alt Transcript"
		snpEffTextForDELWithoutTranscriptWithGene="chr"+"\t".join(snpEffComponents[0:2])+"\t"+snpEffComponents[3] +"\t"+ snpEffComponents[10] ## extracting snpeff text for deletions --> "Chrom Pos Alt Gene"
	
	        if(snpEffComponents[4]=='SNP'):  ## Joining the data from refseq for SNPs in snpeff
			for refseqLine in List2:
                        	refseqLine=refseqLine.rstrip()
                        	refseqCompo=refseqLine.split("\t")
                        	if(refseqCompo[8]!=""):
                                	refseqText="\t".join(refseqCompo[0:2])+"\t"+"\t".join(refseqCompo[3:5])+"\t"+refseqCompo[8]
                        	else:
                                	refseqText="\t".join(refseqCompo[0:2])+"\t"+"\t".join(refseqCompo[3:5])+"\t"+refseqCompo[6]
                        	if(snpEffText==refseqText or snpEffTextWithoutTranscriptWithGene==refseqText):
                                	refseqFlag=1
                                	refseqAnno="\t"+"\t".join(refseqCompo[5:])
		elif(snpEffComponents[4]=='INS'): ## Joining the data from refseq for Insertions in snpeff
			for refseqLine in List2:
				refseqLine=refseqLine.rstrip()
                                refseqCompo=refseqLine.split("\t")
				if(refseqCompo[8]!=""): ## if refseq line contains transcript info then extract the text "Chrom Pos Alt" with "Transcript" else extract only "Chrom Pos Alt"
					refseqText=refseqCompo[0]+"\t"+str(int(refseqCompo[1])+1)+"\t"+"+"+refseqCompo[4]+"\t"+refseqCompo[8] ## extracting from refseq line --> "Chrom Pos Alt Transcript"
				else:
					refseqText=refseqCompo[0]+"\t"+str(int(refseqCompo[1])+1)+"\t"+"+"+refseqCompo[4] +"\t"+ refseqCompo[6] ## extracting from refseq line --> "Chrom Pos Alt Gene"
				if(snpEffTextForINS==refseqText or snpEffTextForINSWithoutTranscriptWithGene==refseqText): ## comparing the snpeff's "Chrom Pos Alt Transcript" text and "Chrom Pos Alt" with refseq text. Only one would be true.
					refseqFlag=1
					refseqAnno="\t"+"\t".join(refseqCompo[5:]) ## storing the respective refseq annotations.
		elif(snpEffComponents[4]=='DEL'): ## Joining the data from refseq for Deletions in snpeff
			for refseqLine in List2:
                                refseqLine=refseqLine.rstrip()
                                refseqCompo=refseqLine.split("\t")
                                if(refseqCompo[8]!=""):
                                        refseqText=refseqCompo[0]+"\t"+refseqCompo[1]+"\t"+"-"+refseqCompo[3]+"\t"+refseqCompo[8]
                                else:
                                        refseqText=refseqCompo[0]+"\t"+refseqCompo[1]+"\t"+"-"+refseqCompo[3]+"\t"+refseqCompo[6]
                                if(snpEffTextForDEL==refseqText or snpEffTextForDELWithoutTranscriptWithGene==refseqText):
                                        refseqFlag=1
                                        refseqAnno="\t"+"\t".join(refseqCompo[5:])


                if(refseqFlag==0):
                        refseqAnno=("\t"+"")*16


		for VarLine in List3: ## Joins the data from in-house 1000g database
			VarLine=VarLine.rstrip()
			VarLineComponents=VarLine.split("\t")
			VarText="chr"+"\t".join(VarLineComponents[0:2])+"\t"+"\t".join(VarLineComponents[3:5])
			if(snpEffTextWithoutTranscript == VarText):
				VarFlag=1
				VarAnno="\t"+"\t".join(VarLineComponents[12:20])
				break
		if(VarFlag==0):
			VarAnno=("\t"+"")*8

		
		for DepthLine in List4: ## Joins the data from DepthForAnnotation file
                	DepthLine=DepthLine.rstrip()
                	DepthLineComponents=DepthLine.split("\t")
                	DepthLineText="\t".join(DepthLineComponents[0:2])+"\t"+"\t".join(DepthLineComponents[3:5])
                	if(snpEffTextWithoutTranscript==DepthLineText):
                        	DepthFlag=1
                        	DepthList=DepthLineComponents[7].split(";")
				DepthAnno="\t"+"\t".join(DepthList)
				DepthAnno=string.replace(DepthAnno,"FAD=","")
				DepthAnno=string.replace(DepthAnno,"FRD=","")
				DepthAnno=string.replace(DepthAnno,"FAF=","")
				DepthAnno=string.replace(DepthAnno,"FDP=","")
				break
        	if(DepthFlag==0):
                	DepthAnno=("\t"+"")*4


		for AlalmutLine in List5: ## Joins the data from in-house alalmut database
			AlalmutLine=AlalmutLine.rstrip()
			if((string.find(AlalmutLine,"#"))==-1):
				AlalmutComponents=AlalmutLine.split("\t")
				AlalmutText="\t".join(AlalmutComponents[0:2])+"\t"+"\t".join(AlalmutComponents[3:5])+"\t"+AlalmutComponents[9].split(".")[0]
				if(snpEffText==AlalmutText):
					AlalmutFlag=1
					AlalmutAnno="\t"+"\t".join(AlalmutComponents[6:])
					break
		if(AlalmutFlag==0):
			AlalmutAnno=("\t"+"")*88+"\t0"
			

		OutFile.write(SampleNameBreakdown[0]+"\t"+snpEffAnno+refseqAnno+VarAnno+DepthAnno+AlalmutAnno+"\n") ## Writes all the joined annotations together
		

OutFile.close()
File1.close()
File2.close()
File3.close()
File4.close()
File5.close()
