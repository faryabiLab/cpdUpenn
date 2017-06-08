##V1.2.2.a
## Run this script as: python SampleProcessing.py Samplesheet_Tobeprocessed  
## This script takes the SampleSheet_Tobeprocessed file as an input and generates the Pipeline Script for each sample in the respective sub-directory.
## updates - can handle V1.1 and V1.2 of bed files
## works with 185 bp pair end reads with adjusted primers for downstream primer in FFPE samples
## This version calls insertion and deletions on pre filtered bam files, removes all variants less than 6bp then merges those with the gatk results later. 
## this Version also removes more of the intermediate files:  *.sam, *.bam, *.btrim.summary, *fastq (uncommressed), *.pe (used in btrim)
## Realignment is now done for the original sorted.fixed.bam and the resulting realigned bam is further used for the analysis
## Non-Canon annotated snpeff is not generated and hence the Non-Canon Annotated,i.e. the final joined (.annotated) file is not generated.
## The PipeScript is not generated for Non-CEBPA samples as this version of Pipe also includes the realignment step which is not preferred to be performed for those samples at this point of time.

import sys,shutil
SAMPLE_FILE=open(sys.argv[1],'r')

    
Samples=SAMPLE_FILE.readlines()
Count=0
for eachsample in Samples:
    Count=Count+1
    eachsample=eachsample.rstrip()
    SampleInfo=eachsample.split(" ")
    SampleName=SampleInfo[0]
    Barcode=SampleInfo[1]+"."+SampleInfo[2]
    Description=SampleInfo[3]
    if(Description=="FFPE_V1"):
    	PrimersFile="Primers/PrimersFFPE185bp.txt"
    	PrimersRCFile="Primers/Primers_RC_FFPE185bp.txt"
    	eParam="140"
    	FragSize="140,50"
    elif(Description=="HEME_V1"):
	FragSize="190,50"
	PrimersFile="Primers/PrimersHEME.txt"
    	PrimersRCFile="Primers/Primers_RC_HEME.txt"
    	eParam="190"
    elif(Description=="FFPE_V1.1"):
	PrimersFile="Primers/PrimersFFPE185bp.txt"
	PrimersRCFile="Primers/Primers_RC_FFPE185bp.txt"
	eParam="140"
    	FragSize="140,50"
    elif(Description=="HEME_V1.1"):
    	FragSize="190,50"
    	PrimersFile="Primers/PrimersHEME.txt"
        PrimersRCFile="Primers/Primers_RC_HEME.txt"
    	eParam="190"
    elif(Description=="HEME_V1.2"):
    	FragSize="190,50"
    	PrimersFile="Primers/PrimersHEMEv1.2.txt"
        PrimersRCFile="Primers/Primers_RC_HEMEv1.2.txt"
    	eParam="190"
    elif(Description=="CEBPa"):
    	FragSize="200,150"
    
    SampleCount="S"+ str(Count)
    
    PipeScriptFileName=SampleName+"/"+SampleName+".PipeScript"
    PipelineScript_FILE=open(PipeScriptFileName,'w')
    
# Phase 1 - Mapping, Filtering, Clipping and Stats Generation


    if(Description!="CEBPa"):
    	continue
    	#BED_FILE="TSCA_CANCER_PANEL_"+Description+"_target.BED"
    	#RGLB_String="TSCA_"+Description
    	#PipelineScript_FILE.write("gunzip "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.fastq.gz"+"\n")
    	#PipelineScript_FILE.write("gunzip "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.fastq.gz"+"\n")
    	#PipelineScript_FILE.write("/home/bdaber01/Btrim64 -p "+PrimersFile+" -t "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.fastq -o "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.btrim.fastq -s "+SampleName+"/"+SampleName+"_R2.btrim.summary.txt -w 0 -e "+eParam+" > "+SampleName+"/"+SampleName+"_R2.btrim.stats.txt"+"\n")
    	#PipelineScript_FILE.write("/home/bdaber01/Btrim64 -p "+PrimersRCFile+" -t "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.fastq -o "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.btrim.fastq -s "+SampleName+"/"+SampleName+"_R1.btrim.summary.txt -w 0 -e "+eParam+" > "+SampleName+"/"+SampleName+"_R1.btrim.stats.txt"+"\n")
    	#PipelineScript_FILE.write("perl paired_end_trim.pl "+SampleName+"/"+SampleName+"_R1.btrim.summary.txt "+SampleName+"/"+SampleName+"_R2.btrim.summary.txt "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.btrim.fastq "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.btrim.fastq"+"\n")
    	#PipelineScript_FILE.write("gzip "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.fastq"+"\n")
    	#PipelineScript_FILE.write("gzip "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.fastq"+"\n")
    	#PipelineScript_FILE.write("echo 'Carrying out Alignment of fastq files with hg19 reference genome using Novoalign...\n\n'"+"\n")
    	#PipelineScript_FILE.write("/home/bdaber01/novocraftV3.0.5/novoalign -d /Drive_D/Databases/genomes/hg19_genome.nix -f "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.btrim.fastq.pe "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.btrim.fastq.pe -i PE "+FragSize+" -c 6 -o FullNW -o SAM > "+SampleName+"/"+SampleName+".novo.sam" +"\n")
    elif(Description=="CEBPa"):
    	BED_FILE="NextEra_"+Description+"_target.BED"
    	RGLB_String="NextEra_"+Description
    	Adapters="CTGTCTCTTATACACATCTGACGCTGCCGACG"+" "+"CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
    	PipelineScript_FILE.write("echo 'Carrying out Alignment of fastq files with hg19 reference genome using Novoalign...\n\n'"+"\n")
    	PipelineScript_FILE.write("/home/bdaber01/novocraftV3.0.5/novoalign -d /Drive_D/Databases/genomes/hg19_genome.nix -f "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.fastq.gz "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.fastq.gz -i PE "+FragSize+" -c 6 -o FullNW -o SAM -a "+str(Adapters)+" > "+SampleName+"/"+SampleName+".novo.sam" +"\n")
    
    PipelineScript_FILE.write("echo 'Converting sam to bam using samtools...\n\n'"+"\n")
    PipelineScript_FILE.write("samtools view -bS "+SampleName+"/"+SampleName+".novo.sam > "+SampleName+"/"+SampleName+".novo.bam"+"\n")
    PipelineScript_FILE.write("echo 'Sorting initial bam...\n\n'"+"\n")
    PipelineScript_FILE.write("samtools sort "+SampleName+"/"+SampleName+".novo.bam "+SampleName+"/"+SampleName+".novo.sorted" + "\n")
    PipelineScript_FILE.write("echo 'Fixing the read group for sorted bam using Picard Tools...\n\n'"+"\n")
    PipelineScript_FILE.write("java -Xmx2g -jar /home/bdaber01/picard-tools-1.46/AddOrReplaceReadGroups.jar I= "+SampleName+"/"+SampleName+".novo.sorted.bam O= "+SampleName+"/"+SampleName+".novo.sorted.fixed.bam RGID=1 RGLB="+RGLB_String+" RGPL=Illumina RGPU="+Barcode+" RGSM="+SampleName+" RGCN=CPD" + "\n")
    PipelineScript_FILE.write("echo 'Indexing fixed bam using Samtools..\n\n'"+"\n")
    PipelineScript_FILE.write("samtools index "+SampleName+"/"+SampleName+".novo.sorted.fixed.bam" +"\n")
    PipelineScript_FILE.write("echo 'Performing Realignment on fixed bam..\n\n'"+"\n")
    PipelineScript_FILE.write("java -Xmx2g -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -T RealignerTargetCreator -I "+SampleName+"/"+SampleName+".novo.sorted.fixed.bam -R /Drive_D/Databases/genomes/hg19_genome.fa -o "+SampleName+"/"+SampleName+".forIndelRealigner.Intervals -nt 4" + "\n")
    PipelineScript_FILE.write("java -Xmx2g -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -T IndelRealigner -I "+SampleName+"/"+SampleName+".novo.sorted.fixed.bam -R /Drive_D/Databases/genomes/hg19_genome.fa -targetIntervals "+SampleName+"/"+SampleName+".forIndelRealigner.Intervals -o "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.bam" +"\n") 
    PipelineScript_FILE.write("echo 'Creating Flagstat for fixed.realigned bam...\n\n'"+"\n")
    PipelineScript_FILE.write("samtools flagstat "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.bam > "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.flagstat"+"\n")
    PipelineScript_FILE.write("echo 'Indexing fixed.realigned bam using Samtools..\n\n'"+"\n")
    PipelineScript_FILE.write("samtools index "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.bam" +"\n")
    PipelineScript_FILE.write("echo 'Intersecting target bed with fixed.realigned bam...\n\n'"+"\n")
    PipelineScript_FILE.write("intersectBed -abam "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.bam -b "+BED_FILE+" > "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ontarget.bam" + "\n")
    PipelineScript_FILE.write("echo 'Indexing on target fixed.realigned bam using Samtools..\n\n'"+"\n")
    PipelineScript_FILE.write("samtools index "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ontarget.bam" +"\n")
    PipelineScript_FILE.write("echo 'Creating Flagstat for sorted.fixed.realigned.ontarget.bam...\n\n'"+"\n")
    PipelineScript_FILE.write("samtools flagstat "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ontarget.bam > "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ontarget.flagstat"+"\n")
    PipelineScript_FILE.write("echo 'Samtools filtering of ontarget bam using minimum mapping quality cutoff of 40 for reads...\n\n'"+"\n")
    PipelineScript_FILE.write("samtools view -h -q 40 "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ontarget.bam > "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ontarget.q40.sam"+"\n")
    PipelineScript_FILE.write("echo 'Running FilterAlignScore.py to further filter the reads using minimum alignment score cutoff of 95...\n\n'"+"\n")
    PipelineScript_FILE.write("python FilterAlignScore.py "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ontarget.q40.sam "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ontarget.q40.as95.sam" + "\n")
    PipelineScript_FILE.write("echo 'Converting filtered sam to bam using Samtools...\n\n'"+"\n")
    PipelineScript_FILE.write("samtools view -bS "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ontarget.q40.as95.sam > "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ontarget.afterFilter.bam" + "\n")
    PipelineScript_FILE.write("echo 'Creating Flagstat for filtered bam...\n\n'"+"\n")
    PipelineScript_FILE.write("samtools flagstat "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ontarget.afterFilter.bam > "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ontarget.afterFilter.flagstat" + "\n")
    PipelineScript_FILE.write("echo 'Indexing filtered bam using Samtools..\n\n'"+"\n")
    PipelineScript_FILE.write("samtools index "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ontarget.afterFilter.bam" +"\n")
    PipelineScript_FILE.write("echo 'Soft Clipping the end of reads...\n\n'"+"\n")
    PipelineScript_FILE.write("java -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -T ClipReads -I "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ontarget.afterFilter.bam  -o "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ontarget.afterFilter.clipped.bam -os "+SampleName+"/"+SampleName+".clippingstats.txt -R /Drive_D/Databases/genomes/hg19_genome.fa -CR SOFTCLIP_BASES -QT 22"+"\n")
    PipelineScript_FILE.write("echo 'Indexing clipped bam using Samtools..\n\n'"+"\n")
    PipelineScript_FILE.write("samtools index "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ontarget.afterFilter.clipped.bam" +"\n")
    PipelineScript_FILE.write("echo 'Running GATK DepthOfCoverage..\n\n'"+"\n")  
    PipelineScript_FILE.write("java -Xmx2g -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -T DepthOfCoverage -o "+SampleName+"/"+SampleName+".Depth -I "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ontarget.afterFilter.clipped.bam --minBaseQuality 22 -baseCounts -ct 0 -ct 1 -ct 250 -ct 1000 -L "+BED_FILE+" -R /Drive_D/Databases/genomes/hg19_genome.fa" + "\n")
    PipelineScript_FILE.write("echo 'Running MakingOfGenotypeGivenAlleles.sh script to create the GenotypeGivenAlleles file , the .Depth.below250X.vcf file and the .DepthForAnnotation File..\n\n'"+"\n")
    PipelineScript_FILE.write("./MakingOfGenotypeGivenAlleles.sh "+SampleName+" /home/bdaber01/Drive_D/Databases/RefBases/hg19.ref."+Description+"\n")
    PipelineScript_FILE.write("echo 'Intersecting the .Depth.below250X.vcf with the target bed file to get the target regions having reads with depth below 250..\n\n'"+"\n")
    PipelineScript_FILE.write("intersectBed -a "+BED_FILE+"_refseq_intersect.BED -b "+SampleName+"/"+SampleName+".Depth.below250X.vcf -u > "+SampleName+"/"+SampleName+".Depth.below250X.amplicons" + "\n")
    PipelineScript_FILE.write("echo 'Creating Depth.below250X.amplicons.summary file..\n\n'"+"\n")
    PipelineScript_FILE.write("python ampliconSummary.py "+SampleName+"\n")
    PipelineScript_FILE.write("echo 'Creating Stats file for each sample..\n\n'"+"\n")
    PipelineScript_FILE.write("python SummaryStats.py "+SampleName+" "+Description+"\n")
    PipelineScript_FILE.write("python MakingOfCoverageSummary.py "+SampleName+"\n")
    
    
    
#Phase 2 - Variant Calling

    PipelineScript_FILE.write("echo 'Running GATK Variant Call in Discovery mode..\n\n'"+"\n")
    PipelineScript_FILE.write("java -Xmx2G -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -T UnifiedGenotyper -I "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ontarget.afterFilter.clipped.bam --dbsnp /home/bdaber01/Drive_D/Databases/dbsnp/hg19_dbsnp137.vcf -o "+SampleName+"/"+SampleName+".GATK.vcf -nt 4 -glm BOTH -mbq 22 -dt NONE -minIndelCnt 3 -minIndelFrac 0.04 -nda -maxAlleles 5 -stand_emit_conf 10.0 -L "+BED_FILE + "\n")
    PipelineScript_FILE.write("echo 'Running GATK Variant Call in Discovery mode on PRE_FILTERED Target Bam file..\n\n'"+"\n")
    PipelineScript_FILE.write("java -Xmx2G -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -T UnifiedGenotyper -I "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ontarget.bam --dbsnp /home/bdaber01/Drive_D/Databases/dbsnp/hg19_dbsnp137.vcf -o "+SampleName+"/"+SampleName+".GATK.pre_filtered.vcf -nt 4 -glm INDEL -mbq 22 -dt NONE -minIndelCnt 3 -minIndelFrac 0.04 -nda -maxAlleles 5 -stand_emit_conf 10.0 -L "+BED_FILE + "\n")
    PipelineScript_FILE.write("echo 'Remove Small Indels (i.e. less than 6bp) from the GATK.pre_filtered.vcf file ..\n\n'"+"\n")
    PipelineScript_FILE.write("python RemoveSmallIndels.py "+SampleName+"/"+SampleName+".GATK.pre_filtered.vcf "+SampleName+"/"+SampleName+".GATK.pre_filtered.greaterThan6bpIndels.vcf "+"\n")
    PipelineScript_FILE.write("echo 'Running GATK Variant Call in GenotypeGivenAlleles mode..\n\n'"+"\n")
    PipelineScript_FILE.write("java -Xmx3G -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -T UnifiedGenotyper -I "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ontarget.afterFilter.clipped.bam --dbsnp /home/bdaber01/Drive_D/Databases/dbsnp/hg19_dbsnp137.vcf -o "+SampleName+"/"+SampleName+".GATK.af4percent.vcf -nt 4 -glm SNP -mbq 22 -dt NONE -alleles:VCF "+SampleName+"/"+SampleName+".GenotypeGivenAlleles.sorted.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_emit_conf 10.0 -L "+BED_FILE + "\n")
    PipelineScript_FILE.write("echo 'Combining three GATK Variant Calls to give a combined.vcf file..\n\n'"+"\n")
    PipelineScript_FILE.write("java -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -R /Drive_D/Databases/genomes/hg19_genome.fa -T CombineVariants --variant "+SampleName+"/"+SampleName+".GATK.vcf --variant "+SampleName+"/"+SampleName+".GATK.af4percent.vcf --variant "+SampleName+"/"+SampleName+".GATK.pre_filtered.greaterThan6bpIndels.vcf -o "+SampleName+"/"+SampleName+".combined.vcf"+"\n")
    PipelineScript_FILE.write("echo 'Filtering (tagging) Variants called at Homopolymer sites..\n\n'"+"\n")
    PipelineScript_FILE.write("java -Xmx2G -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -T VariantFiltration --variant "+SampleName+"/"+SampleName+".combined.vcf -o "+SampleName+"/"+SampleName+".combined.filtered.vcf -dt NONE --filterName 'HomopolymerRun' --filterExpression 'HRun>6'" + "\n")
    PipelineScript_FILE.write("echo 'grepping all variants not called at Homoploymer sites for further annotation..\n\n'"+"\n")
    PipelineScript_FILE.write("grep -v 'HomopolymerRun' "+SampleName+"/"+SampleName+".combined.filtered.vcf > "+SampleName+"/"+SampleName+".combined.filtered.HRunRemoved.vcf"+"\n")
    PipelineScript_FILE.write("echo 'Running SnpEff Annotation-- both with canon and without canon option..\n\n'"+"\n")
    PipelineScript_FILE.write("java -jar /home/bdaber01/snpEff_3_0/snpEff.jar -c /home/bdaber01/snpEff_3_0/snpEff.config -canon hg19 -i vcf -o txt -noStats "+SampleName+"/"+SampleName+".combined.filtered.HRunRemoved.vcf > "+SampleName+"/"+SampleName+".combined.filtered.HRunRemoved.canon_snpeff.txt"+"\n")
    #PipelineScript_FILE.write("java -jar /home/bdaber01/snpEff_3_0/snpEff.jar -c /home/bdaber01/snpEff_3_0/snpEff.config hg19 -i vcf -o txt -noStats "+SampleName+"/"+SampleName+".combined.filtered.HRunRemoved.vcf > "+SampleName+"/"+SampleName+".combined.filtered.HRunRemoved.snpeff.txt"+"\n")
    PipelineScript_FILE.write("echo 'Running MultiAllelic.py script on the filtered.HRunRemoved.vcf to convert all mutiallelic position to single allelic ones..\n\n'"+"\n")
    PipelineScript_FILE.write("python MultiAllelic.py "+SampleName+"/"+SampleName+".combined.filtered.HRunRemoved.vcf "+SampleName+"/"+SampleName+".combined.filtered.HRunRemoved.multiallelic.vcf"+"\n")
    PipelineScript_FILE.write("echo 'Running Annovar Refseq Annotation..\n\n'"+"\n")
    PipelineScript_FILE.write("perl /home/bdaber01/annovar/convert2annovar.pl  -format vcf4 -includeinfo "+SampleName+"/"+SampleName+".combined.filtered.HRunRemoved.multiallelic.vcf > "+SampleName+"/"+SampleName+".combined.filtered.HRunRemoved.multiallelic.avinput"+"\n")
    PipelineScript_FILE.write("perl /home/bdaber01/annovar/annotate_variation.pl "+SampleName+"/"+SampleName+".combined.filtered.HRunRemoved.multiallelic.avinput /home/bdaber01/annovar/humandb/ -buildver hg19"+"\n")
    PipelineScript_FILE.write("python refseq.py "+SampleName+"/"+SampleName+".combined.filtered.HRunRemoved.multiallelic.avinput.variant_function "+SampleName+"/"+SampleName+".combined.filtered.HRunRemoved.multiallelic.avinput.exonic_variant_function "+SampleName+"/"+SampleName+".combined.filtered.HRunRemoved.multiallelic.avinput.refseq"+"\n")
    PipelineScript_FILE.write("echo 'Joining All Annotations..\n\n'"+"\n")
    PipelineScript_FILE.write("python AnnotationJoin.py "+SampleName+"/"+SampleName+".combined.filtered.HRunRemoved.canon_snpeff.txt "+SampleName+"/"+SampleName+".combined.filtered.HRunRemoved.multiallelic.avinput.refseq /Drive_D/Databases/1000Genomes/hg19.1000genomes.extractFields2.annovar_input "+SampleName+"/"+SampleName+".DepthForAnnotation.sorted.vcf /Drive_D/Databases/Alalmut/hg19.alalmut "+SampleName+"/"+SampleName+".canon_annotated"+"\n")
    #PipelineScript_FILE.write("python AnnotationJoin.py "+SampleName+"/"+SampleName+".combined.filtered.HRunRemoved.snpeff.txt "+SampleName+"/"+SampleName+".combined.filtered.HRunRemoved.multiallelic.avinput.refseq /Drive_D/Databases/1000Genomes/hg19.1000genomes.extractFields2.annovar_input "+SampleName+"/"+SampleName+".DepthForAnnotation.sorted.vcf /Drive_D/Databases/Alalmut/hg19.alalmut "+SampleName+"/"+SampleName+".annotated"+"\n")
    PipelineScript_FILE.write('find . -type f -name "*.sam" -exec rm -f {} \;'+"\n") # will remove the intermediate sam files to save space
    PipelineScript_FILE.write('find . -type f -name "*.fastq" -exec rm -f {} \;'+"\n") # will remove the intermediate fastq files after decompression for btrim. The compressed original files are saved
    PipelineScript_FILE.write('find . -type f -name "*.pe" -exec rm -f {} \;'+"\n") # will remove the intermediate .pe files created by btrim
    PipelineScript_FILE.write("mv "+SampleName+"/"+SampleName+".canon_annotated "+"."+"\n")
    PipelineScript_FILE.write("mv "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ontarget.afterFilter.clipped.ba* "+"."+"\n")
    PipelineScript_FILE.write("mv "+SampleName+"/"+SampleName+".novo.sorted.fixed.realigned.ba* "+"."+"\n") 
    PipelineScript_FILE.write('find '+SampleName+'/'+' -type f -name "*.bam" -exec rm -f {} \;'+"\n")  # will remove the intermediate Bam files.
    PipelineScript_FILE.write('find . -type f -name "*.btrim.summary.txt" -exec rm -f {} \;'+"\n")
    PipelineScript_FILE.close()

SAMPLE_FILE.close()



    
    
    
    
     
    
    
    
    
    
    
    
    
    
    
    
