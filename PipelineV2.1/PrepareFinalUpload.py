#V2.1
## Run this script as : python PrepareFinalUpload.py <Run_masterVarFinal file> <Run_masterVarFinal_HEME_V2.3 file> <Run_masterVarFinal_UP file> <Run_masterVarFinal_HEME_V2.3_UP file> <Run_HEME_SamplePairsInfo file>
## Changes in V2.1 : a) Included header in Run_masterVarFinal_HEME_V2.3_UP file.
import sys
import subprocess
##Defining the variables and fileHandles
Run_masterVarFinal_file=sys.argv[1]
Run_masterVarFinal_HEME_V2_3_file=sys.argv[2]
Run_masterVarFinal_UP_fileHandle=open(sys.argv[3],'w')
Run_masterVarFinal_HEME_V2_3_UP_fileHandle=open(sys.argv[4],'w')
Run_HEME_SamplePairsInfo_fileHandle=open(sys.argv[5],'r')
All_HEMEao_V1_3_Samples=''
Run_masterVarFinal_HEME_V2_3_UP_fileHandle.write("SampleName\tChrom\tPos\tRef\tAlt\tVariant Type(SnpEff)\tHomozygous\tQuality\tCoverage\tWarnings\tGene_ID\tGene_name\tBio_type\tTrancript_ID\tExon_ID\tExon_Rank\tEffect\told_AA/new_AA\tOld_codon/New_codon\tCodon_Num(CDS)\tCodon_Degeneracy\tCDS_size\tCodons_around\tAAs_around\tCustom_interval_ID\tc_Change(SnpEff)\tp_Change(SnpEff)\tLocation(annovar)\tGene\tConsequence\tTranscript\tcDNA change\tProtien Change\tCHROM(vcf)\tPOS\tdbsnp137_ID\tRef\tAlt\tQUAL\tFILTER\tINFO\tGenotypeFormat\tGenotype\tAC(1000g)\tAN(1000g)\tAF(1000g)\tAMR_AF\tASN_AF\tEUR_AF\tAFR_AF\tVariantType\tFDP\tFRD\tFAD\tFAF\tgene\tgeneId\ttranscript\tstrand\ttransLen\tprotein\tUniprot\tvarType\tcodingEffect\tvarLocation\tassembly\tgDNAstart\tgDNAend\tgNomen\tcDNAstart\tcDNAend\tcNomen\tpNomen\talt_pNomen\texon\tintron\tomimId\tdistNearestSS\tnearestSSType\twtSSFScore\twtMaxEntScore\twtNNSScore\twtGSScore\twtHSFScore\tvarSSFScore\tvarMaxEntScore\tvarNNSScore\tvarGSScore\tvarHSFScore\tnearestSSChange\tlocalSpliceEffect\tproteinDomain1\tproteinDomain2\tproteinDomain3\tproteinDomain4\trsId\trsValidated\trsSuspect\trsValidations\trsValidationNumber\trsAncestralAllele\trsHeterozygosity\trsClinicalSignificance\trsMAF\trsMAFAllele\trsMAFCount\t1000g_AF\t1000g_AFR_AF\t1000g_SAS_AF\t1000g_EAS_AF\t1000g_EUR_AF\t1000g_AMR_AF\texacAllFreq\texacAFRFreq\texacAMRFreq\texacEASFreq\texacSASFreq\texacNFEFreq\texacFINFreq\texacOTHFreq\texacAFRHmz\texacAMRHmz\texacEASHmz\texacSASHmz\texacNFEHmz\texacFINHmz\texacOTHHmz\texacFilter\texacReadDepth\tespRefEACount\tespRefAACount\tespRefAllCount\tespAltEACount\tespAltAACount\tespAltAllCount\tespEAMAF\tespAAMAF\tespAllMAF\tespEAAAF\tespAAAAF\tespAllAAF\tespAvgReadDepth\tclinVarIds\tclinVarOrigins\tclinVarMethods\tclinVarClinSignifs\tclinVarReviewStatus\tclinVarPhenotypes\thgmdId\thgmdPhenotype\thgmdPubMedId\thgmdSubCategory\tcosmicIds\tcosmicTissues\tcosmicFreqs\tcosmicSampleCounts\tinsNucs\tdelNucs\tsubstType\twtNuc\tvarNuc\tnucChange\tphastCons\tphyloP\twtAA_1\twtAA_3\twtCodon\twtCodonFreq\tvarAA_1\tvarAA_3\tvarCodon\tvarCodonFreq\tposAA\tnOrthos\tconservedOrthos\tconservedDistSpecies\tBLOSUM45\tBLOSUM62\tBLOSUM80\twtAAcomposition\tvarAAcomposition\twtAApolarity\tvarAApolarity\twtAAvolume\tvarAAvolume\tgranthamDist\tAGVGDclass\tAGVGDgv\tAGVGDgd\tSIFTprediction\tSIFTweight\tSIFTmedian\tMAPPprediction\tMAPPpValue\tMAPPpValueMedian\tAdjacentSNPflag"+"\n")
## Running a while loop over the Run_HEME_SamplePairsInfo file to extract the sample-names from each line and perform subsequent functions further
while 1:
	SamplePairLine=Run_HEME_SamplePairsInfo_fileHandle.readline()
	if(SamplePairLine):
		SamplePairLine=SamplePairLine.rstrip()
		AllSampleNames=SamplePairLine.split(",")
		HEME_V1_2_Sample,HEMEao_V1_3_Sample,HEME_V2_3_Sample=AllSampleNames[0],AllSampleNames[1],AllSampleNames[2]
		All_HEMEao_V1_3_Samples=All_HEMEao_V1_3_Samples+HEMEao_V1_3_Sample+"|"      ##Making a big string of HEMEao_V1.3 sample-names with "|" characters in between them
		proc1=subprocess.Popen(['sed','s/'+HEME_V2_3_Sample+'/'+HEMEao_V1_3_Sample+'/g',Run_masterVarFinal_HEME_V2_3_file],stdout=subprocess.PIPE) ## Replacing HEME_V2.3 sample-names with their respective HEMEao_V1.3 sample-names in the Run_masterVarFinal_HEME_V2.3 file and storing all the variant entries with replaced names in a new file named Run_masterVarFinal_HEME_V2_3_UP file. The original Run_masterVarFinal_HEME_V2.3 file is still retained.
		proc2=subprocess.Popen(['grep',HEMEao_V1_3_Sample],stdin=proc1.stdout,stdout=subprocess.PIPE).communicate()[0]
		Run_masterVarFinal_HEME_V2_3_UP_fileHandle.write(proc2)	
		proc1.stdout.close()
	else:
		break

All_HEMEao_V1_3_Samples=All_HEMEao_V1_3_Samples.strip("|")  ## Eliminating the last "|" character from the string to make work the logic in 'egrep -v' command correctly
proc3=subprocess.Popen(['egrep','-v',All_HEMEao_V1_3_Samples,Run_masterVarFinal_file],stdout=subprocess.PIPE).communicate()[0]   ##Storing all the variants except those for HEMEao_V1.3 samples in Run_masterVarFinal_UP file
Run_masterVarFinal_UP_fileHandle.write(proc3)

Run_HEME_SamplePairsInfo_fileHandle.close()
Run_masterVarFinal_UP_fileHandle.close()
Run_masterVarFinal_HEME_V2_3_UP_fileHandle.close()

