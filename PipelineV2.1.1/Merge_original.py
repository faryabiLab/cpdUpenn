##V2.1
## Run this script as python Merge.py SampleSheet_forMergeAnalysis.txt
## Changes in V1.4:    1) Changed the script name from 'ampliconSummary.py' to 'ampliconSummary_forMergeData.py' in the respective command line.
##                     2) Changed 'SampleInfo[3]' to 'SampleInfo[2]' to make sure the Barcode_string is calculated as 'SampleInfo[2]+"."+SampleInfo[3]'.
## Changes in V1.4.2:  1) This script now also makes a new file named "Run_HEME_SamplePairsInfo.txt" which would contain the names of HEME_V1.2 sample, its respective HEMEao_V1.3 pair and the final merged sample (i.e.HEME_V2.3 sample).
## Changes in V2.0:    1) Added sleep 1m command between variant calling commands.
##                     2) Added the commands to make Depth.below150x.amplicons_sorted file.
##                     3) Variants with HRun>6 and  having allele-frequency>=10% are now retained.
##                     4) samtools, bedtools, snpEff and Annovar tools updated to their new versions.
##		       5) AnnotationJoin.py's command line was changed to replace hg19.alamut database with hg19.alamut_v1.4.2.DescriptionP1_DescriptionP2(i.e. HEME_V1.2_HEMEao_V1.3) to reflect the most recent version of database specific to HEME panel.
##                     6) SamplePairsInfoFile.close() command added at the end of the script.
##                     7) As the new snpeff version yeilds only vcf output, code-lines (involving both snpEff and snpSift {part of snpEff tool}) have been added to process that output vcf to yeild a similar type of tab-delimited text file (same as before but with two additional columns) so that the downstream joining of annotations from other sources/tools, does not get affected much. This is also done to be concordant with the format of variant entries in CPD's KnowledgeBase.
##Changes in V2.1:             1) GATK version upgraded to 3.4-46. HomopolymerRun Tagging not supported anymore. Hence, removed codes lines related to GATK's HRun Filter and FilterOutHomopolymerVar.py. Therefore, filenames would no more contain the word "HRunTagged" in them.
##                             2) Added -Djava.io.tmpdir=tmp parameter to each GATK command line that did not have it earlier. The argument "-maxAlleles" changed to "-maxAltAlleles" on respective GATK command lines to suit the new version of GATK. The arguments '-K CPD_upenn_gatk.key' and '-et NO_ET' added to each command line to prevent uploading the GATK-run-stats-report to the AWS server.
##                             3) Changed '-nt 6' to '-nct 6' and '-Xmx12G' to '-Xmx6G' on the command line of GATK's UnifiedGenotyper in order to improve the performance of the process.
##                             4) -genotypeMergeOptions and -priority arguments added to command line of GATK's CombineVariants tool. The version of GATK now requires that these arguments should be mentioned on commandline. Changed '--variant' to '--variant:variant1 (or variant2 or variant3 or variant4)' on the "CombineVariants" command lines to represent different variant-call-files (i.e. vcfs).

import sys              ## importing sys library to be used further in the script
import subprocess
SAMPLE_FILE=open(sys.argv[1],'r')  ##Defining a file-handle named "SAMPLE_FILE" to handle the file named "SampleSheet_forMergeAnalysis.txt" that is opened using open command

##Reading all the lines one by one from the SAMPLE_FILE handle and storing them into a list named "Samples"
Samples=SAMPLE_FILE.readlines()
Count=0

##Defining all lists to be used later in the code
Counts=[]
SampleShortNames=[]
SampleNames=[]
Descriptions=[]
Barcodes=[]
SamplePairs=[]
BarcodePairs=[]
SampleShortNamePairs=[]
ActualDescriptionPairs=[]

##Starting a for loop over the lines contained in the "Samples" list
for eachsample in Samples:
	Count=Count+1
	Counts.append(Count)
	eachsample=eachsample.rstrip()
	SampleInfo=eachsample.split(" ")
	SampleShortNames.append(SampleInfo[0])
	SampleNames.append(SampleInfo[1])
	Descriptions.append(SampleInfo[4]+","+SampleInfo[0])
	Barcode_string=SampleInfo[2]+"."+SampleInfo[3]
	Barcodes.append(Barcode_string)

##Starting a loop over the Descriptions list to search for HEME_V1.2 samples
for eachDescription in Descriptions:
	if((eachDescription.find("HEME_V1.2"))!=-1):
		Description=eachDescription
		Index=Descriptions.index(Description)       ## Extracting the list-index (from the "Descriptions" list) for the Description that is detected to be HEME_V1.2
		ActualDescription=Description.split(",")[0]
		SampleShortName=SampleShortNames[Index]         ## Extracting the SampleShortName for the Sample annotated with HEME_V1.2 description
		SampleName=SampleNames[Index]                   ## Extracting the SampleName for the Sample annotated with HEME_V1.2 description
		Barcode=Barcodes[Index]                         ## Extracting the Barcode for the Sample annotated with HEME_V1.2 description
		searchText=SampleShortName+"-"+"ao"            ## Defining the searchText to search for the current Sample's HemeAddOn partner
		##Running an inner for loop to search for the current Sample's HemeAddOn Partener in the SampleShortNames list
		for l in SampleShortNames:
			if(l==searchText):
				aoIndex=SampleShortNames.index(l)
				##Confirming the Description annotated for the HemeAddOn Partner
				if((Descriptions[aoIndex].find("HEMEao_V1.3"))!=-1):
					ActualDescription_ao=Descriptions[aoIndex].split(",")[0]
					SampleName_ao=SampleNames[aoIndex]
					Barcode_ao=Barcodes[aoIndex]
					SampleShortName_ao=l
					SamplePairs.append(SampleName+","+SampleName_ao)
					BarcodePairs.append(Barcode+","+Barcode_ao)
					SampleShortNamePairs.append(SampleShortName+","+SampleShortName_ao)
					ActualDescriptionPairs.append(ActualDescription+","+ActualDescription_ao)
				else:
					sys.exit("Error in the SampleSheet_forMergeAnalysis.txt")
			else:
				continue
	else:
		continue

SamplePairsInfoFile=open("Run_HEME_SamplePairsInfo.txt",'w')
for eachSamplePair in SamplePairs:
	SampleP1=eachSamplePair.split(",")[0]
	SampleP2=eachSamplePair.split(",")[1]
	eachSamplePairIndex=SamplePairs.index(eachSamplePair)
	eachActualDescriptionPair=ActualDescriptionPairs[eachSamplePairIndex]
	DescriptionP1=eachActualDescriptionPair.split(",")[0]
	DescriptionP2=eachActualDescriptionPair.split(",")[1]
	if(DescriptionP1=="HEME_V1.2" and DescriptionP2=="HEMEao_V1.3"):
		LibraryP1="TSCA"
		LibraryP2="TSCA"
		minIndelFrac="0.04"
		minIndelCount="3"
	BED_FILE_for_P1=LibraryP1+"_"+"CANCER_PANEL_"+DescriptionP1+"_target.BED"
	BED_FILE_for_P1_PrimersOn=LibraryP1+"_"+"CANCER_PANEL_"+DescriptionP1+"_target_PrimersOn.BED"
	BED_FILE_for_P2=LibraryP2+"_"+"CANCER_PANEL_"+DescriptionP2+"_target.BED"
	BED_FILE_for_P2_PrimersOn=LibraryP2+"_"+"CANCER_PANEL_"+DescriptionP2+"_target_PrimersOn.BED"
	RGLB_String=LibraryP1+"_"+DescriptionP1+"_"+DescriptionP2
	RGPU_string=BarcodePairs[eachSamplePairIndex]
	ShortNameForP1=SampleShortNamePairs[eachSamplePairIndex].split(",")[0]
	MergeDirectoryName=ShortNameForP1+"_merged"
	SamplePairsInfoFile.write(SampleP1+","+SampleP2+","+MergeDirectoryName+"\n")
	process1=subprocess.Popen(['mkdir',MergeDirectoryName],stdout=subprocess.PIPE).communicate()[0]
	MergePipeScriptFileName=MergeDirectoryName+"/"+MergeDirectoryName+"."+"PipeScript"
	PipelineScript_FILE=open(MergePipeScriptFileName,'w')
	PipelineScript_FILE.write("cp "+SampleP1+".novo.sorted.fixed.ba*"+" "+MergeDirectoryName+"\n")
	PipelineScript_FILE.write("cp "+SampleP1+".novo.sorted.fixed.ontarget.afterFilter.clipped.ba*"+" "+MergeDirectoryName+"\n")
	PipelineScript_FILE.write("cp "+SampleP1+".*.ontarget.ba*"+" "+MergeDirectoryName+"\n")
	PipelineScript_FILE.write("cp "+SampleP2+".novo.sorted.fixed.ba*"+" "+MergeDirectoryName+"\n")
	PipelineScript_FILE.write("cp "+SampleP2+".novo.sorted.fixed.ontarget.afterFilter.clipped.ba*"+" "+MergeDirectoryName+"\n")
	PipelineScript_FILE.write("cp "+SampleP2+".*.ontarget.ba*"+" "+MergeDirectoryName+"\n")
	PipelineScript_FILE.write("cp "+SampleP1+"/"+SampleP1+"*.fastq.gz"+" "+MergeDirectoryName+"\n")
	PipelineScript_FILE.write("cp "+SampleP2+"/"+SampleP2+"*.fastq.gz"+" "+MergeDirectoryName+"\n")
	PipelineScript_FILE.write("zcat "+MergeDirectoryName+"/"+"*_L001_R1_001.fastq.gz"+" "+">"+MergeDirectoryName+"/"+MergeDirectoryName+"_L001_R1_001.fastq"+"\n")
	PipelineScript_FILE.write("zcat "+MergeDirectoryName+"/"+"*_L001_R2_001.fastq.gz"+" "+">"+MergeDirectoryName+"/"+MergeDirectoryName+"_L001_R2_001.fastq"+"\n")
	
	##Merging Phase
	PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools merge "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.MERGED.bam "+MergeDirectoryName+"/"+SampleP1+".novo.sorted.fixed.bam "+MergeDirectoryName+"/"+SampleP2+".novo.sorted.fixed.bam "+"\n")
	PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools merge "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.afterFilter.clipped.MERGED.bam "+MergeDirectoryName+"/"+SampleP1+".novo.sorted.fixed.ontarget.afterFilter.clipped.bam "+MergeDirectoryName+"/"+SampleP2+".novo.sorted.fixed.ontarget.afterFilter.clipped.bam "+"\n")
	PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools merge "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.MERGED.bam "+MergeDirectoryName+"/"+SampleP1+".novo.sorted.fixed.ontarget.bam "+MergeDirectoryName+"/"+SampleP2+".novo.sorted.fixed.ontarget.bam "+"\n")
	PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools merge "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.primer.sorted.fixed.ontarget.MERGED.bam "+MergeDirectoryName+"/"+SampleP1+".novo.primer.sorted.fixed.ontarget.bam "+MergeDirectoryName+"/"+SampleP2+".novo.primer.sorted.fixed.ontarget.bam "+"\n")
	PipelineScript_FILE.write("java -Xmx2g -jar /home/bdaber01/picard-tools-1.46/AddOrReplaceReadGroups.jar I= "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.MERGED.bam O= "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.MERGED.fixed.bam RGID=1 RGLB="+RGLB_String+" RGPL=Illumina RGPU="+RGPU_string+" RGSM="+MergeDirectoryName+" RGCN=CPD" + "\n")
	PipelineScript_FILE.write("java -Xmx2g -jar /home/bdaber01/picard-tools-1.46/AddOrReplaceReadGroups.jar I= "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.MERGED.bam O= "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.MERGED.fixed.bam RGID=1 RGLB="+RGLB_String+" RGPL=Illumina RGPU="+RGPU_string+" RGSM="+MergeDirectoryName+" RGCN=CPD" + "\n")
	PipelineScript_FILE.write("java -Xmx2g -jar /home/bdaber01/picard-tools-1.46/AddOrReplaceReadGroups.jar I= "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.primer.sorted.fixed.ontarget.MERGED.bam O= "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.primer.sorted.fixed.ontarget.MERGED.fixed.bam RGID=1 RGLB="+RGLB_String+" RGPL=Illumina RGPU="+RGPU_string+" RGSM="+MergeDirectoryName+" RGCN=CPD" + "\n")
	PipelineScript_FILE.write("java -Xmx2g -jar /home/bdaber01/picard-tools-1.46/AddOrReplaceReadGroups.jar I= "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.afterFilter.clipped.MERGED.bam O= "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.afterFilter.clipped.MERGED.fixed.bam RGID=1 RGLB="+RGLB_String+" RGPL=Illumina RGPU="+RGPU_string+" RGSM="+MergeDirectoryName+" RGCN=CPD" + "\n")
	PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools flagstat "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.MERGED.fixed.bam > "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.MERGED.fixed.flagstat" + "\n")
	PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools flagstat "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.afterFilter.clipped.MERGED.fixed.bam > "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.afterFilter.clipped.MERGED.fixed.flagstat" + "\n")
	PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools flagstat "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.MERGED.fixed.bam > "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.MERGED.fixed.flagstat" + "\n")
	PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools flagstat "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.primer.sorted.fixed.ontarget.MERGED.fixed.bam > "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.primer.sorted.fixed.ontarget.MERGED.fixed.flagstat" + "\n")
	PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools index "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.MERGED.fixed.bam"+"\n")
	PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools index "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.afterFilter.clipped.MERGED.fixed.bam"+"\n")
	PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools index "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.MERGED.fixed.bam"+"\n")
	PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools index "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.primer.sorted.fixed.ontarget.MERGED.fixed.bam"+"\n")
	PipelineScript_FILE.write("echo 'Running GATK DepthOfCoverage on merged bam..\n\n'"+"\n")
	PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -Xmx2g -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -T DepthOfCoverage -K CPD_upenn_gatk.key -et NO_ET -o "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth -I "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.afterFilter.clipped.MERGED.fixed.bam --minBaseQuality 22 -baseCounts -ct 0 -ct 1 -ct 250 -ct 1000 -L "+BED_FILE_for_P1+" -L "+BED_FILE_for_P2+" -R /Drive_D/Databases/genomes/hg19_genome.fa" + "\n")
	PipelineScript_FILE.write("echo 'Running MakingOfGenotypeGivenAlleles.sh script to create the GenotypeGivenAlleles file , the .Depth.below250X.vcf file and the .DepthForAnnotation File..\n\n'"+"\n")
	PipelineScript_FILE.write("./MakingOfGenotypeGivenAlleles.sh "+MergeDirectoryName+" /home/bdaber01/Drive_D/Databases/RefBases/hg19.ref."+DescriptionP1+"_"+DescriptionP2+"\n")
	PipelineScript_FILE.write("echo 'Intersecting the .Depth.below250X.vcf with the target bed file to get the target regions having reads with depth below 250..\n\n'"+"\n")
	PipelineScript_FILE.write("/home/bdaber01/bedtools-2.25.0/bin/intersectBed -a "+BED_FILE_for_P1+"_refseq_intersect.BED -b "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below250X.vcf -u > "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below250X."+DescriptionP1+".amplicons" + "\n")
	PipelineScript_FILE.write("/home/bdaber01/bedtools-2.25.0/bin/intersectBed -a "+BED_FILE_for_P2+"_refseq_intersect.BED -b "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below250X.vcf -u > "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below250X."+DescriptionP2+".amplicons" + "\n")
	PipelineScript_FILE.write("(cat "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below250X."+DescriptionP1+".amplicons"+";"+"cat "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below250X."+DescriptionP2+".amplicons"+")"+" > "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below250X.amplicons"+"\n")
	PipelineScript_FILE.write("sort "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below250X.amplicons"+" > "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below250X.amplicons_sorted"+"\n")
	PipelineScript_FILE.write("/home/bdaber01/bedtools-2.25.0/bin/intersectBed -a "+BED_FILE_for_P1+"_refseq_intersect.BED -b "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below150X.vcf -u > "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below150X."+DescriptionP1+".amplicons" + "\n")
	PipelineScript_FILE.write("/home/bdaber01/bedtools-2.25.0/bin/intersectBed -a "+BED_FILE_for_P2+"_refseq_intersect.BED -b "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below150X.vcf -u > "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below150X."+DescriptionP2+".amplicons" + "\n")
	PipelineScript_FILE.write("(cat "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below150X."+DescriptionP1+".amplicons"+";"+"cat "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below150X."+DescriptionP2+".amplicons"+")"+" > "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below150X.amplicons"+"\n")
	PipelineScript_FILE.write("sort "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below150X.amplicons"+" > "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below150X.amplicons_sorted"+"\n")
	PipelineScript_FILE.write("echo 'Creating Depth.below250X(and 150X).amplicons_sorted.summary file..\n\n'"+"\n")
	PipelineScript_FILE.write("python ampliconSummary_forMergeData.py "+MergeDirectoryName+"\n")
	PipelineScript_FILE.write("python MakingOfCoverageSummary.py "+MergeDirectoryName+" "+DescriptionP1+"_"+DescriptionP2+"\n")
	
	##Variant Calling
	
	PipelineScript_FILE.write("echo 'Running GATK Variant Call in Discovery mode..\n\n'"+"\n")
	PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -Xmx6G -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -K CPD_upenn_gatk.key -et NO_ET -T UnifiedGenotyper -I "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.afterFilter.clipped.MERGED.fixed.bam --dbsnp /home/bdaber01/Drive_D/Databases/dbsnp/hg19_dbsnp137_correct.vcf -o "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.vcf -nct 6 -glm BOTH -mbq 22 -dt NONE -minIndelCnt "+minIndelCount+" -minIndelFrac "+minIndelFrac+" -nda -maxAltAlleles 5 -stand_emit_conf 10.0 -L "+BED_FILE_for_P1+" -L "+BED_FILE_for_P2+"\n")
	PipelineScript_FILE.write("echo 'waiting for 1 minute before executing next command'"+"\n")
	PipelineScript_FILE.write("sleep 1m"+"\n")
	# calling 6bp or larger on a pre-filtered bam file not containing primers
	PipelineScript_FILE.write("echo 'Running GATK Variant Call in Discovery mode on PRE_FILTERED Target Merged Bam file..\n\n'"+"\n")
	PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -Xmx6G -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -K CPD_upenn_gatk.key -et NO_ET -T UnifiedGenotyper -I "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.MERGED.fixed.bam --dbsnp /home/bdaber01/Drive_D/Databases/dbsnp/hg19_dbsnp137_correct.vcf -o "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.pre_filtered.vcf -nct 6 -glm INDEL -mbq 22 -dt NONE -minIndelCnt "+minIndelCount+" -minIndelFrac "+minIndelFrac+" -nda -maxAltAlleles 5 -stand_emit_conf 10.0 -L "+BED_FILE_for_P1+" -L "+BED_FILE_for_P2+"\n")
	PipelineScript_FILE.write("echo 'waiting for 1 minute before executing next command'"+"\n")
	PipelineScript_FILE.write("sleep 1m"+"\n")
	PipelineScript_FILE.write("echo 'Remove Small Indels (i.e. less than 6bp) from the GATK.pre_filtered.vcf file ..\n\n'"+"\n")
	PipelineScript_FILE.write("python RemoveSmallIndels.py "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.pre_filtered.vcf "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.pre_filtered.greaterThan6bpIndels.vcf "+"\n")
	
	# calling 6bp or larger on a pre-filtered bam file containing primers
	PipelineScript_FILE.write("echo 'Running GATK Variant Call in Discovery mode on PRIMERS-ON-PRE_FILTERED Target Bam file..\n\n'"+"\n")
	PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -Xmx6G -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -K CPD_upenn_gatk.key -et NO_ET -T UnifiedGenotyper -I "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.primer.sorted.fixed.ontarget.MERGED.fixed.bam --dbsnp /home/bdaber01/Drive_D/Databases/dbsnp/hg19_dbsnp137_correct.vcf -o "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.primer.pre_filtered.vcf -nct 6 -glm INDEL -mbq 22 -dt NONE -minIndelCnt "+minIndelCount+" -minIndelFrac "+minIndelFrac+" -nda -maxAltAlleles 5 -stand_emit_conf 10.0 -L "+BED_FILE_for_P1_PrimersOn+" -L "+BED_FILE_for_P2_PrimersOn+"\n")
	PipelineScript_FILE.write("echo 'waiting for 1 minute before executing next command'"+"\n")
	PipelineScript_FILE.write("sleep 1m"+"\n")
	PipelineScript_FILE.write("echo 'Remove Small Indels (i.e. less than 6bp) from the GATK.primer.pre_filtered.vcf file ..\n\n'"+"\n")
	PipelineScript_FILE.write("python RemoveSmallIndels.py "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.primer.pre_filtered.vcf "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.primer.pre_filtered.greaterThan6bpIndels.vcf "+"\n")
	
	# calling SNVs in GenotypeGivenAlleles mode
	PipelineScript_FILE.write("echo 'Running GATK Variant Call in GenotypeGivenAlleles mode..\n\n'"+"\n")
	PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -Xmx6G -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -K CPD_upenn_gatk.key -et NO_ET -T UnifiedGenotyper -I "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.afterFilter.clipped.MERGED.fixed.bam --dbsnp /home/bdaber01/Drive_D/Databases/dbsnp/hg19_dbsnp137_correct.vcf -o "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.af4percent.vcf -nct 6 -glm SNP -mbq 22 -dt NONE -alleles:VCF "+MergeDirectoryName+"/"+MergeDirectoryName+".GenotypeGivenAlleles.sorted.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_emit_conf 10.0 -L "+BED_FILE_for_P1+" -L "+BED_FILE_for_P2+"\n")
	PipelineScript_FILE.write("echo 'waiting for 1 minute before executing next command'"+"\n")
	PipelineScript_FILE.write("sleep 1m"+"\n")
	
	# combining all variant calls together
	PipelineScript_FILE.write("echo 'Combining all four GATK Variant Calls to give a combined.vcf file..\n\n'"+"\n")
	PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -Xmx2G -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -R /Drive_D/Databases/genomes/hg19_genome.fa -T CombineVariants -K CPD_upenn_gatk.key -et NO_ET --variant:variant1 "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.vcf --variant:variant2 "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.af4percent.vcf --variant:variant3 "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.pre_filtered.greaterThan6bpIndels.vcf --variant:variant4 "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.primer.pre_filtered.greaterThan6bpIndels.vcf -genotypeMergeOptions PRIORITIZE -priority variant3,variant4,variant2,variant1 -o "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.vcf"+"\n")
	
	# Filtering the variants
	
	#PipelineScript_FILE.write("echo 'Filtering (tagging) Variants called at Homopolymer sites..\n\n'"+"\n")
	#PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -Xmx2G -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -K CPD_upenn_gatk.key -et NO_ET -T VariantFiltration --variant "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.vcf -o "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.vcf -dt NONE --filterName 'HomopolymerRun' --filterExpression 'HRun>6'" + "\n")
	#PipelineScript_FILE.write("echo 'grepping all variants not called at Homoploymer sites for further annotation..\n\n'"+"\n")
	#PipelineScript_FILE.write("grep -v 'HomopolymerRun' "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.filtered.vcf > "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.vcf"+"\n")
	
	# Annotating the variants
	
	PipelineScript_FILE.write("echo 'Running SnpEff Annotation--with canon option activated..\n\n'"+"\n")
	PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -jar /home/bdaber01/snpEff_v_4.1/snpEff.jar -c /home/bdaber01/snpEff_v_4.1/snpEff.config -canon hg19 -i vcf -o vcf -noStats "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.vcf > "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.canon_snpeff.vcf"+"\n")
	PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -jar /home/bdaber01/snpEff_v_4.1/SnpSift.jar varType "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.canon_snpeff.vcf > "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.canon_snpeff.vT.vcf"+"\n")
	PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -jar /home/bdaber01/snpEff_v_4.1/SnpSift.jar extractFields "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.canon_snpeff.vT.vcf CHROM POS ID REF ALT QUAL FILTER DP ANN HET HOM VARTYPE >"+MergeDirectoryName+"/"+MergeDirectoryName+".combined.canon_snpeff.vT.ex.txt"+"\n")
	PipelineScript_FILE.write("python PrepareSnpEffOut.py "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.canon_snpeff.vT.ex.txt "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.canon_snpeff.vT.ex.tab"+"\n")
	PipelineScript_FILE.write("echo 'Running MultiAllelic.py script on the filtered.HRunRemoved.vcf to convert all mutiallelic position to single allelic ones..\n\n'"+"\n")
	PipelineScript_FILE.write("python MultiAllelic.py "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.vcf "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.multiallelic.vcf"+"\n")
	PipelineScript_FILE.write("echo 'Running Annovar Refseq Annotation..\n\n'"+"\n")
	PipelineScript_FILE.write("perl /home/bdaber01/annovar_17thJun2015/convert2annovar.pl  -format vcf4old -includeinfo "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.multiallelic.vcf > "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.multiallelic.avinput"+"\n")
	PipelineScript_FILE.write("perl /home/bdaber01/annovar_17thJun2015/annotate_variation.pl "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.multiallelic.avinput /home/bdaber01/annovar_17thJun2015/humandb/ -buildver hg19"+"\n")
	PipelineScript_FILE.write("python refseq.py "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.multiallelic.avinput.variant_function "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.multiallelic.avinput.exonic_variant_function "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.multiallelic.avinput.refseq"+"\n")
	PipelineScript_FILE.write("echo 'Joining All Annotations..\n\n'"+"\n")
	PipelineScript_FILE.write("python AnnotationJoin.py "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.canon_snpeff.vT.ex.tab "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.multiallelic.avinput.refseq /Drive_D/Databases/1000Genomes/hg19.1000genomes.extractFields2.annovar_input "+MergeDirectoryName+"/"+MergeDirectoryName+".DepthForAnnotation.sorted.vcf /Drive_D/Databases/Alalmut/hg19.alamut_v1.4.2."+DescriptionP1+"_"+DescriptionP2+" "+MergeDirectoryName+"/"+MergeDirectoryName+".canon_annotated"+"\n")
	PipelineScript_FILE.write("echo 'Flagging Adjacent Snps..\n\n'"+"\n")
	PipelineScript_FILE.write("python flagAdjacentSNPs.py "+MergeDirectoryName+"/"+MergeDirectoryName+".canon_annotated "+MergeDirectoryName+"/"+MergeDirectoryName+".canon_annotated_flagged "+"\n")
	#PipelineScript_FILE.write("python FilterOutHomopolymerVar.py "+MergeDirectoryName+"/"+MergeDirectoryName+".canon_annotated_flagged "+MergeDirectoryName+"/"+MergeDirectoryName+".canon_annotated_flagged.HRR "+"\n")
	# Creating SummaryStats File
	
	PipelineScript_FILE.write("echo 'Creating Stats file for the sample..\n\n'"+"\n")
	PipelineScript_FILE.write("python SummaryStats_forMergeData.py "+MergeDirectoryName+"\n")
	##PipelineScript_FILE.write('find '+SampleName+'/'+' -type f -name "*.sam" -exec rm -f {} \;'+"\n") # will remove the intermediate sam files to save space
	##PipelineScript_FILE.write('find '+SampleName+'/'+' -type f -name "*.fastq" -exec rm -f {} \;'+"\n") # will remove the intermediate fastq files after decompression for btrim. The compressed original files are saved
	##PipelineScript_FILE.write('find '+SampleName+'/'+' -type f -name "*.pe" -exec rm -f {} \;'+"\n") # will remove the intermediate .pe files created by btrim
	PipelineScript_FILE.write("mv "+MergeDirectoryName+"/"+MergeDirectoryName+".canon_annotated* "+"."+"\n")
	PipelineScript_FILE.write("mv "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.MERGED.fixed.ba* "+"."+"\n")
	PipelineScript_FILE.write("mv "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.primer.sorted.fixed.ontarget.MERGED.fixed.ba* "+"."+"\n")
	PipelineScript_FILE.write("mv "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.afterFilter.clipped.MERGED.fixed.ba* "+"."+"\n")
	PipelineScript_FILE.write("mv "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.MERGED.fixed.ba* "+"."+"\n")
	##PipelineScript_FILE.write('find '+SampleName+'/'+' -type f -name "*.bam" -exec rm -f {} \;'+"\n")  # will remove the intermediate Bam files
	##PipelineScript_FILE.write('find '+SampleName+'/'+' -type f -name "*.btrim.summary.txt" -exec rm -f {} \;'+"\n")
	PipelineScript_FILE.close()

SAMPLE_FILE.close()
SamplePairsInfoFile.close()
	
    	
  
