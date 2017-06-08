##V1.4.2
## Run this script as python Merge.py SampleSheet_forMergeAnalysis.txt
## Changes in V1.4:    1) Changed the script name from 'ampliconSummary.py' to 'ampliconSummary_forMergeData.py' in the respective command line.
##                     2) Changed 'SampleInfo[3]' to 'SampleInfo[2]' to make sure the Barcode_string is calculated as 'SampleInfo[2]+"."+SampleInfo[3]'.
## Changes in V1.4.2:  1) This script now also makes a new file named "Run_HEME_SamplePairsInfo.txt" which would contain the names of HEME_V1.2 sample, its respective HEMEao_V1.3 pair and the final merged sample (i.e.HEME_V2.3 sample).

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
	PipelineScript_FILE.write("samtools merge "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.MERGED.bam "+MergeDirectoryName+"/"+SampleP1+".novo.sorted.fixed.bam "+MergeDirectoryName+"/"+SampleP2+".novo.sorted.fixed.bam "+"\n")
	PipelineScript_FILE.write("samtools merge "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.afterFilter.clipped.MERGED.bam "+MergeDirectoryName+"/"+SampleP1+".novo.sorted.fixed.ontarget.afterFilter.clipped.bam "+MergeDirectoryName+"/"+SampleP2+".novo.sorted.fixed.ontarget.afterFilter.clipped.bam "+"\n")
	PipelineScript_FILE.write("samtools merge "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.MERGED.bam "+MergeDirectoryName+"/"+SampleP1+".novo.sorted.fixed.ontarget.bam "+MergeDirectoryName+"/"+SampleP2+".novo.sorted.fixed.ontarget.bam "+"\n")
	PipelineScript_FILE.write("samtools merge "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.primer.sorted.fixed.ontarget.MERGED.bam "+MergeDirectoryName+"/"+SampleP1+".novo.primer.sorted.fixed.ontarget.bam "+MergeDirectoryName+"/"+SampleP2+".novo.primer.sorted.fixed.ontarget.bam "+"\n")
	PipelineScript_FILE.write("java -Xmx2g -jar /home/bdaber01/picard-tools-1.46/AddOrReplaceReadGroups.jar I= "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.MERGED.bam O= "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.MERGED.fixed.bam RGID=1 RGLB="+RGLB_String+" RGPL=Illumina RGPU="+RGPU_string+" RGSM="+MergeDirectoryName+" RGCN=CPD" + "\n")
	PipelineScript_FILE.write("java -Xmx2g -jar /home/bdaber01/picard-tools-1.46/AddOrReplaceReadGroups.jar I= "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.MERGED.bam O= "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.MERGED.fixed.bam RGID=1 RGLB="+RGLB_String+" RGPL=Illumina RGPU="+RGPU_string+" RGSM="+MergeDirectoryName+" RGCN=CPD" + "\n")
	PipelineScript_FILE.write("java -Xmx2g -jar /home/bdaber01/picard-tools-1.46/AddOrReplaceReadGroups.jar I= "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.primer.sorted.fixed.ontarget.MERGED.bam O= "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.primer.sorted.fixed.ontarget.MERGED.fixed.bam RGID=1 RGLB="+RGLB_String+" RGPL=Illumina RGPU="+RGPU_string+" RGSM="+MergeDirectoryName+" RGCN=CPD" + "\n")
	PipelineScript_FILE.write("java -Xmx2g -jar /home/bdaber01/picard-tools-1.46/AddOrReplaceReadGroups.jar I= "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.afterFilter.clipped.MERGED.bam O= "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.afterFilter.clipped.MERGED.fixed.bam RGID=1 RGLB="+RGLB_String+" RGPL=Illumina RGPU="+RGPU_string+" RGSM="+MergeDirectoryName+" RGCN=CPD" + "\n")
	PipelineScript_FILE.write("samtools flagstat "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.MERGED.fixed.bam > "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.MERGED.fixed.flagstat" + "\n")
	PipelineScript_FILE.write("samtools flagstat "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.afterFilter.clipped.MERGED.fixed.bam > "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.afterFilter.clipped.MERGED.fixed.flagstat" + "\n")
	PipelineScript_FILE.write("samtools flagstat "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.MERGED.fixed.bam > "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.MERGED.fixed.flagstat" + "\n")
	PipelineScript_FILE.write("samtools flagstat "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.primer.sorted.fixed.ontarget.MERGED.fixed.bam > "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.primer.sorted.fixed.ontarget.MERGED.fixed.flagstat" + "\n")
	PipelineScript_FILE.write("samtools index "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.MERGED.fixed.bam"+"\n")
	PipelineScript_FILE.write("samtools index "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.afterFilter.clipped.MERGED.fixed.bam"+"\n")
	PipelineScript_FILE.write("samtools index "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.MERGED.fixed.bam"+"\n")
	PipelineScript_FILE.write("samtools index "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.primer.sorted.fixed.ontarget.MERGED.fixed.bam"+"\n")
	PipelineScript_FILE.write("echo 'Running GATK DepthOfCoverage on merged bam..\n\n'"+"\n")
	PipelineScript_FILE.write("java -Xmx2g -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -T DepthOfCoverage -o "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth -I "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.afterFilter.clipped.MERGED.fixed.bam --minBaseQuality 22 -baseCounts -ct 0 -ct 1 -ct 250 -ct 1000 -L "+BED_FILE_for_P1+" -L "+BED_FILE_for_P2+" -R /Drive_D/Databases/genomes/hg19_genome.fa" + "\n")
	PipelineScript_FILE.write("echo 'Running MakingOfGenotypeGivenAlleles.sh script to create the GenotypeGivenAlleles file , the .Depth.below250X.vcf file and the .DepthForAnnotation File..\n\n'"+"\n")
	PipelineScript_FILE.write("./MakingOfGenotypeGivenAlleles.sh "+MergeDirectoryName+" /home/bdaber01/Drive_D/Databases/RefBases/hg19.ref."+DescriptionP1+"_"+DescriptionP2+"\n")
	PipelineScript_FILE.write("echo 'Intersecting the .Depth.below250X.vcf with the target bed file to get the target regions having reads with depth below 250..\n\n'"+"\n")
	PipelineScript_FILE.write("intersectBed -a "+BED_FILE_for_P1+"_refseq_intersect.BED -b "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below250X.vcf -u > "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below250X."+DescriptionP1+".amplicons" + "\n")
	PipelineScript_FILE.write("intersectBed -a "+BED_FILE_for_P2+"_refseq_intersect.BED -b "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below250X.vcf -u > "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below250X."+DescriptionP2+".amplicons" + "\n")
	PipelineScript_FILE.write("(cat "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below250X."+DescriptionP1+".amplicons"+";"+"cat "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below250X."+DescriptionP2+".amplicons"+")"+" > "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below250X.amplicons"+"\n")
	PipelineScript_FILE.write("sort "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below250X.amplicons"+" > "+MergeDirectoryName+"/"+MergeDirectoryName+".Depth.below250X.amplicons_sorted"+"\n")
	PipelineScript_FILE.write("echo 'Creating Depth.below250X.amplicons_sorted.summary file..\n\n'"+"\n")
	PipelineScript_FILE.write("python ampliconSummary_forMergeData.py "+MergeDirectoryName+"\n")
	PipelineScript_FILE.write("python MakingOfCoverageSummary.py "+MergeDirectoryName+" "+DescriptionP1+"_"+DescriptionP2+"\n")
	
	##Variant Calling
	
	PipelineScript_FILE.write("echo 'Running GATK Variant Call in Discovery mode..\n\n'"+"\n")
	PipelineScript_FILE.write("java -Xmx12G -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -T UnifiedGenotyper -I "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.afterFilter.clipped.MERGED.fixed.bam --dbsnp /home/bdaber01/Drive_D/Databases/dbsnp/hg19_dbsnp137_correct.vcf -o "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.vcf -nt 6 -glm BOTH -mbq 22 -dt NONE -minIndelCnt "+minIndelCount+" -minIndelFrac "+minIndelFrac+" -nda -maxAlleles 5 -stand_emit_conf 10.0 -L "+BED_FILE_for_P1+" -L "+BED_FILE_for_P2+"\n")
	# calling 6bp or larger on a pre-filtered bam file not containing primers
	PipelineScript_FILE.write("echo 'Running GATK Variant Call in Discovery mode on PRE_FILTERED Target Merged Bam file..\n\n'"+"\n")
	PipelineScript_FILE.write("java -Xmx12G -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -T UnifiedGenotyper -I "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.MERGED.fixed.bam --dbsnp /home/bdaber01/Drive_D/Databases/dbsnp/hg19_dbsnp137_correct.vcf -o "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.pre_filtered.vcf -nt 6 -glm INDEL -mbq 22 -dt NONE -minIndelCnt "+minIndelCount+" -minIndelFrac "+minIndelFrac+" -nda -maxAlleles 5 -stand_emit_conf 10.0 -L "+BED_FILE_for_P1+" -L "+BED_FILE_for_P2+"\n")
	PipelineScript_FILE.write("echo 'Remove Small Indels (i.e. less than 6bp) from the GATK.pre_filtered.vcf file ..\n\n'"+"\n")
	PipelineScript_FILE.write("python RemoveSmallIndels.py "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.pre_filtered.vcf "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.pre_filtered.greaterThan6bpIndels.vcf "+"\n")
	
	# calling 6bp or larger on a pre-filtered bam file containing primers
	PipelineScript_FILE.write("echo 'Running GATK Variant Call in Discovery mode on PRIMERS-ON-PRE_FILTERED Target Bam file..\n\n'"+"\n")
	PipelineScript_FILE.write("java -Xmx12G -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -T UnifiedGenotyper -I "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.primer.sorted.fixed.ontarget.MERGED.fixed.bam --dbsnp /home/bdaber01/Drive_D/Databases/dbsnp/hg19_dbsnp137_correct.vcf -o "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.primer.pre_filtered.vcf -nt 6 -glm INDEL -mbq 22 -dt NONE -minIndelCnt "+minIndelCount+" -minIndelFrac "+minIndelFrac+" -nda -maxAlleles 5 -stand_emit_conf 10.0 -L "+BED_FILE_for_P1_PrimersOn+" -L "+BED_FILE_for_P2_PrimersOn+"\n")
	PipelineScript_FILE.write("echo 'Remove Small Indels (i.e. less than 6bp) from the GATK.primer.pre_filtered.vcf file ..\n\n'"+"\n")
	PipelineScript_FILE.write("python RemoveSmallIndels.py "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.primer.pre_filtered.vcf "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.primer.pre_filtered.greaterThan6bpIndels.vcf "+"\n")
	
	# calling SNVs in GenotypeGivenAlleles mode
	PipelineScript_FILE.write("echo 'Running GATK Variant Call in GenotypeGivenAlleles mode..\n\n'"+"\n")
	PipelineScript_FILE.write("java -Xmx12G -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -T UnifiedGenotyper -I "+MergeDirectoryName+"/"+MergeDirectoryName+".novo.sorted.fixed.ontarget.afterFilter.clipped.MERGED.fixed.bam --dbsnp /home/bdaber01/Drive_D/Databases/dbsnp/hg19_dbsnp137_correct.vcf -o "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.af4percent.vcf -nt 6 -glm SNP -mbq 22 -dt NONE -alleles:VCF "+MergeDirectoryName+"/"+MergeDirectoryName+".GenotypeGivenAlleles.sorted.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_emit_conf 10.0 -L "+BED_FILE_for_P1+" -L "+BED_FILE_for_P2+"\n")
	
	# combining all variant calls together
	PipelineScript_FILE.write("echo 'Combining all four GATK Variant Calls to give a combined.vcf file..\n\n'"+"\n")
	PipelineScript_FILE.write("java -Xmx2G -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -R /Drive_D/Databases/genomes/hg19_genome.fa -T CombineVariants --variant "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.vcf --variant "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.af4percent.vcf --variant "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.pre_filtered.greaterThan6bpIndels.vcf --variant "+MergeDirectoryName+"/"+MergeDirectoryName+".GATK.primer.pre_filtered.greaterThan6bpIndels.vcf -o "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.vcf"+"\n")
	
	# Filtering the variants
	
	PipelineScript_FILE.write("echo 'Filtering (tagging) Variants called at Homopolymer sites..\n\n'"+"\n")
	PipelineScript_FILE.write("java -Xmx2G -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -T VariantFiltration --variant "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.vcf -o "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.filtered.vcf -dt NONE --filterName 'HomopolymerRun' --filterExpression 'HRun>6'" + "\n")
	PipelineScript_FILE.write("echo 'grepping all variants not called at Homoploymer sites for further annotation..\n\n'"+"\n")
	PipelineScript_FILE.write("grep -v 'HomopolymerRun' "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.filtered.vcf > "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.filtered.HRR.vcf"+"\n")
	
	# Annotating the variants
	
	PipelineScript_FILE.write("echo 'Running SnpEff Annotation--with canon option activated..\n\n'"+"\n")
	PipelineScript_FILE.write("java -jar /home/bdaber01/snpEff_3_0/snpEff.jar -c /home/bdaber01/snpEff_3_0/snpEff.config -canon hg19 -i vcf -o txt -noStats "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.filtered.HRR.vcf > "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.filtered.HRR.canon_snpeff.txt"+"\n")
	PipelineScript_FILE.write("echo 'Running MultiAllelic.py script on the filtered.HRunRemoved.vcf to convert all mutiallelic position to single allelic ones..\n\n'"+"\n")
	PipelineScript_FILE.write("python MultiAllelic.py "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.filtered.HRR.vcf "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.filtered.HRR.multiallelic.vcf"+"\n")
	PipelineScript_FILE.write("echo 'Running Annovar Refseq Annotation..\n\n'"+"\n")
	PipelineScript_FILE.write("perl /home/bdaber01/annovar/convert2annovar.pl  -format vcf4old -includeinfo "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.filtered.HRR.multiallelic.vcf > "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.filtered.HRR.multiallelic.avinput"+"\n")
	PipelineScript_FILE.write("perl /home/bdaber01/annovar/annotate_variation.pl "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.filtered.HRR.multiallelic.avinput /home/bdaber01/annovar/humandb/ -buildver hg19"+"\n")
	PipelineScript_FILE.write("python refseq.py "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.filtered.HRR.multiallelic.avinput.variant_function "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.filtered.HRR.multiallelic.avinput.exonic_variant_function "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.filtered.HRR.multiallelic.avinput.refseq"+"\n")
	PipelineScript_FILE.write("echo 'Joining All Annotations..\n\n'"+"\n")
	PipelineScript_FILE.write("python AnnotationJoin.py "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.filtered.HRR.canon_snpeff.txt "+MergeDirectoryName+"/"+MergeDirectoryName+".combined.filtered.HRR.multiallelic.avinput.refseq /Drive_D/Databases/1000Genomes/hg19.1000genomes.extractFields2.annovar_input "+MergeDirectoryName+"/"+MergeDirectoryName+".DepthForAnnotation.sorted.vcf /Drive_D/Databases/Alalmut/hg19.alalmut "+MergeDirectoryName+"/"+MergeDirectoryName+".canon_annotated"+"\n")
	PipelineScript_FILE.write("echo 'Flagging Adjacent Snps..\n\n'"+"\n")
	PipelineScript_FILE.write("python flagAdjacentSNPs.py "+MergeDirectoryName+"/"+MergeDirectoryName+".canon_annotated "+MergeDirectoryName+"/"+MergeDirectoryName+".canon_annotated_flagged "+"\n")
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
	
    	
  
