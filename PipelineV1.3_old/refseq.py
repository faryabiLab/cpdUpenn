##V1.3
## Run this script as: python refseq.py <input -- annovar's variant function file> <input -- annovar's exonic variant function file> <output -- merged file containing both single variants and their exonic properties>
##Change in V1.3 : IndexError while parsing the exonicNC_list is now handled. This error occured previously for 'block substitutions' that didnt have any protein change reported in the exonic_variant_function file.
##Hence, now if there is no protein change reported for an exonic variant then the corresponding column is remained empty in the final output (i.e. '.refseq' file)


import sys,re, string
File1=open(sys.argv[1],'r')
File2=open(sys.argv[2],'r')
OutFile=open(sys.argv[3],'w')
List1=File1.readlines()
List2=File2.readlines()
for variantLine in List1:
	exonicFlag=0
	variantLine=variantLine.rstrip()
	variantComponents=variantLine.split("\t")
	variantText="\t".join(variantComponents[2:7]) ## Extracts text (Chrom StartPos EndPos ref alt) from variant_function file.
	variantInfo="\t".join(variantComponents[0:2])
	vcfText="\t".join(variantComponents[7:])
	for exonicVarLine in List2:
		exonicVarLine=exonicVarLine.rstrip()
		exonicVarComponents=exonicVarLine.split("\t")
		exonicText="\t".join(exonicVarComponents[3:8])  ## Extracts the text (Chrom StartPos EndPos ref alt)  from exonic_variant_function file.
		if(variantText == exonicText): ## checks for the similar text (i.e. Chrom StartPos EndPos ref alt) between the two files and joins the lines form both the files containing the common text.
			exonicFlag=1
			exonicType=exonicVarComponents[1]         ## extracts the exonic variant type ,i.e. synonymous, non-synonymous, frameshift or stop-gain.
			exonicNC=exonicVarComponents[2]           ##exonicNC is exonic Nomenclature (containing Transcript, NucleotideChange (cChange) and ProteinChange (pChange) for each possible Transcript)
	if(exonicFlag==0): ## if no common text found, then writes the line from variant_function file as it is without any annotations from exonic file (this is the case for introns, splicing sites and UTRs).
		OutFile.write(variantText+"\t"+variantInfo+"\t"+""+"\t"+""+"\t"+""+"\t"+""+"\t"+vcfText+"\n")
	elif(exonicFlag==1):
		exonicNC_List=exonicNC.split(",") ## Each Transcript entry separated by comma before is splitted now on the same and stored in the list.
		for items in exonicNC_List:
			if(items!=""):
				NC_elements=items.split(":") ## For each Transcript entry, the Transcript Name, cChange and pChange info is obtained and printed alongwith its respective exonic entry.
				Transcript=NC_elements[1]
				cChange=NC_elements[3]
				try:
					pChange=NC_elements[4]
				except IndexError, indErr:
					pChange=""
				OutFile.write(variantText+"\t"+variantInfo+"\t"+exonicType+"\t"+Transcript+"\t"+cChange+"\t"+pChange+"\t"+vcfText+"\n") ## This Output file retains all the information from original vcf stored in vcf-text variable before.
			

OutFile.close
File1.close()
File2.close()
