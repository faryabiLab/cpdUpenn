##V1.2.4
## Run this script as: python SampleProcessing.py Samplesheet_Tobeprocessed  
## This script takes the SampleSheet_Tobeprocessed file as an input and generates the Pipeline Script for each sample in the respective sub-directory.
## updates - can handle V1.1 and V1.2 of bed files
## works with 185 bp pair end reads with adjusted primers for downstream primer in FFPE samples
## This version calls insertion and deletions on pre filtered bam files, removes all variants less than 6bp then merges those with the gatk results later. 
## this Version also removes more of the intermediate files:  *.sam, *.bam, *.btrim.summary, *fastq (uncommressed), *.pe (used in btrim)
## Changes in Version 1.2.2 : 1) The novoalign software is updated to version 3.00.05. 
##                            2) The quality cut-off value is changed from 15 to 22 on GATK's ClipReads, DepthOfCoverage and UnifiedGenotyper command-lines in response to the change in 
##                               illumina-Miseq's software version which now assigns higher quality values to the bases as compared to the previous version.
## Changes in Version 1.2.3 : 1)The '-Xmx (maximum memory allocated to Java)' parameter on each GATK-UnifiedGenotyper (variant call) command line was change from '2G' or '3G' to '12G' and
##			      the -nt (no. of data threads) parameter on those command lines was changed from '4' to '6' consequently. This means that we would now be running
##			      each GATK-UnifiedGenotyper command on 6 data threads i.e. 6 cpu cores (1 data thread on each core), assigning 2GB memory to each core (hence utilizing total memory = 2GB * 6cores = 12GB ).
##			      This change would speed up the Variant Calling steps and hence would yeild the pipe results faster than before.
## Changes in Version 1.2.4 : 1)Introducing HEMEv2.0 Panel.
##			      2)Replacing "Samtools sort" with "Novosort".
##                            3)'-Xmx' parameter with value of 2G was mentioned on GATK's CombineVariants command line for the first time just to make sure enough memory is made available to java during the execution of the program.
##                            4)No more generation of non-canon snpeff file.
##                            5)The parameter '-Djava.io.tmpdir=tmp' was included on GATK's DepthOfCoverage, ClipReads and UnifiedGenotyper CommandLines to make sure temp files are directed to the 'tmp' directory under current RunFolder under DriveD.
##                              Previously they were directed to the '/tmp' folder on the linux server which could easily run out of space. DriveD has a lot of space to handle the temp files.
##			      6)'Description' was added to the MakingOfCoverageSummary.py command line.
##			      7) Indel Depth Calculation now included in the Pipe.
##			      8) Signature Genotype is now generated for each sample using genotypeSig.py script and is displayed in the Final StatSummary file.
##			      9) Adjacent Snps are now flagged for each sample. The flags are displayed in the last field(or column) of the cannon_annotated_flagged file.
##			      10) Output File names would now contain 'HRR' instead of 'HRunRemoved'.
##			      11) The generated cannon_annotated.flagged is moved to the main folder after analysis.
##                            12) Annovar software was updated. The value for option -format was changed from 'vcf4' to 'vcf4old' on convert2annovar command-line to make sure that none of the variants were lost.
##                                With the new version of convert2annovar.pl script if we keep the value of -format option as 'vcf4' then it would eliminate the variants having genotype '0/0' and we dont want that
##                                as we dont trust the genotyping by GATK for some low frequency variants.
##                            13) Finding and Removing of the files at the end of processing of each sample is now done by using the sample-specific folder name on the command line instead of asking unix to find the files randomly anywhere in the main run folder.
##                                This modification is done to adopt a safer practice.
##			      14) Corrected the hg19.ref.FFPE_V1.1 and hg19.ref.FFPE_V1 files. The new file names are hg19.ref.FFPE_V1.1_correct and hg19.ref.FFPE_V1_correct. An "if" condition has been incorporated in the 
##				  script at "MakingGenotypeGivenAlleles" step to use these corrected files. Correction was done at the chromosome positions chr7:128850461-128850472 to remove the error (i.e. '^M' character)
##                                present in the ref-allele column. This error had no impact on the results reported using the previous pipeline versions as it was found in the intronic region of SMO gene (i.e. chr7:128850461-128850472).
##                            15) The hg19_dbsnp137.vcf is no more used as it had errors at some chromosome positions which had "chrchr" instead of "chr".This error has been rectified and the new file used in the pipe is named as
##                                hg19_dbsnp137_correct.vcf.


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
    elif(Description=="HEME_V2.0"):
    	FragSize="190,50"
    	PrimersFile="Primers/PrimersHEMEv2.0.txt"
    	PrimersRCFile="Primers/Primers_RC_HEMEv2.0.txt"
    	eParam="190"
    else:
    	continue
    BED_FILE="TSCA_CANCER_PANEL_"+Description+"_target.BED"
    RGLB_String="TSCA_"+Description
    SampleCount="S"+ str(Count)
    
    PipeScriptFileName=SampleName+"/"+SampleName+".PipeScript"
    PipelineScript_FILE=open(PipeScriptFileName,'w')
    
# Phase 1 - Mapping, Filtering, Clipping and Stats Generation  
    
    PipelineScript_FILE.write("gunzip "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.fastq.gz"+"\n")
    PipelineScript_FILE.write("gunzip "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.fastq.gz"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/Btrim64 -p "+PrimersFile+" -t "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.fastq -o "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.btrim.fastq -s "+SampleName+"/"+SampleName+"_R2.btrim.summary.txt -w 0 -e "+eParam+" > "+SampleName+"/"+SampleName+"_R2.btrim.stats.txt"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/Btrim64 -p "+PrimersRCFile+" -t "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.fastq -o "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.btrim.fastq -s "+SampleName+"/"+SampleName+"_R1.btrim.summary.txt -w 0 -e "+eParam+" > "+SampleName+"/"+SampleName+"_R1.btrim.stats.txt"+"\n")
    PipelineScript_FILE.write("perl paired_end_trim.pl "+SampleName+"/"+SampleName+"_R1.btrim.summary.txt "+SampleName+"/"+SampleName+"_R2.btrim.summary.txt "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.btrim.fastq "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.btrim.fastq"+"\n")
    PipelineScript_FILE.write("gzip "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.fastq"+"\n")
    PipelineScript_FILE.write("gzip "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.fastq"+"\n")
    
    	
    PipelineScript_FILE.write("echo 'Carrying out Alignment of fastq files with hg19 reference genome using Novoalign...\n\n'"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/novocraftV3.0.5/novoalign -d /Drive_D/Databases/genomes/hg19_genome.nix -f "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.btrim.fastq.pe "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.btrim.fastq.pe -i PE "+FragSize+" -c 6 -o FullNW -o SAM > "+SampleName+"/"+SampleName+".novo.sam" +"\n")
    PipelineScript_FILE.write("echo 'Converting sam to bam using samtools...\n\n'"+"\n")
    PipelineScript_FILE.write("samtools view -bS "+SampleName+"/"+SampleName+".novo.sam > "+SampleName+"/"+SampleName+".novo.bam"+"\n")
    PipelineScript_FILE.write("echo 'Sorting initial bam using Novosort...\n\n'"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/novocraftV3.0.5/novosort -t tmp -c 6 "+SampleName+"/"+SampleName+".novo.bam > "+SampleName+"/"+SampleName+".novo.sorted.bam"+"\n")
    PipelineScript_FILE.write("echo 'Fixing the read group for sorted bam using Picard Tools...\n\n'"+"\n")
    PipelineScript_FILE.write("java -Xmx2g -jar /home/bdaber01/picard-tools-1.46/AddOrReplaceReadGroups.jar I= "+SampleName+"/"+SampleName+".novo.sorted.bam O= "+SampleName+"/"+SampleName+".novo.sorted.fixed.bam RGID=1 RGLB="+RGLB_String+" RGPL=Illumina RGPU="+Barcode+" RGSM="+SampleName+" RGCN=CPD" + "\n")
    PipelineScript_FILE.write("echo 'Creating Flagstat for bam...\n\n'"+"\n")
    PipelineScript_FILE.write("samtools flagstat "+SampleName+"/"+SampleName+".novo.sorted.fixed.bam > "+SampleName+"/"+SampleName+".novo.sorted.flagstat"+"\n")
    PipelineScript_FILE.write("echo 'Indexing fixed bam using Samtools..\n\n'"+"\n")
    PipelineScript_FILE.write("samtools index "+SampleName+"/"+SampleName+".novo.sorted.fixed.bam" +"\n")
    PipelineScript_FILE.write("echo 'Intersecting target bed with bam...\n\n'"+"\n")
    PipelineScript_FILE.write("intersectBed -abam "+SampleName+"/"+SampleName+".novo.sorted.fixed.bam -b "+BED_FILE+" > "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.bam" + "\n")
    PipelineScript_FILE.write("echo 'Indexing on target fixed bam using Samtools..\n\n'"+"\n")
    PipelineScript_FILE.write("samtools index "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.bam" +"\n")
    PipelineScript_FILE.write("echo 'Creating Flagstat for sorted.fixed.ontarget.bam...\n\n'"+"\n")
    PipelineScript_FILE.write("samtools flagstat "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.bam > "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.flagstat"+"\n")
    PipelineScript_FILE.write("echo 'Samtools filtering of ontarget bam using minimum mapping quality cutoff of 40 for reads...\n\n'"+"\n")
    PipelineScript_FILE.write("samtools view -h -q 40 "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.bam > "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.q40.sam"+"\n")
    PipelineScript_FILE.write("echo 'Running FilterAlignScore.py to further filter the reads using minimum alignment score cutoff of 95...\n\n'"+"\n")
    PipelineScript_FILE.write("python FilterAlignScore.py "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.q40.sam "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.q40.as95.sam" + "\n")
    PipelineScript_FILE.write("echo 'Converting filtered sam to bam using Samtools...\n\n'"+"\n")
    PipelineScript_FILE.write("samtools view -bS "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.q40.as95.sam > "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.bam" + "\n")
    PipelineScript_FILE.write("echo 'Creating Flagstat for filtered bam...\n\n'"+"\n")
    PipelineScript_FILE.write("samtools flagstat "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.bam > "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.flagstat" + "\n")
    PipelineScript_FILE.write("echo 'Indexing filtered bam using Samtools..\n\n'"+"\n")
    PipelineScript_FILE.write("samtools index "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.bam" +"\n")
    PipelineScript_FILE.write("echo 'Soft Clipping the end of reads...\n\n'"+"\n")
    PipelineScript_FILE.write("java -Xmx2g -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -T ClipReads -I "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.bam  -o "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.clipped.bam -os "+SampleName+"/"+SampleName+".clippingstats.txt -R /Drive_D/Databases/genomes/hg19_genome.fa -CR SOFTCLIP_BASES -QT 22"+"\n")
    PipelineScript_FILE.write("echo 'Indexing clipped bam using Samtools..\n\n'"+"\n")
    PipelineScript_FILE.write("samtools index "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.clipped.bam" +"\n")
    PipelineScript_FILE.write("echo 'Running GATK DepthOfCoverage..\n\n'"+"\n")  
    PipelineScript_FILE.write("java -Xmx2g -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -T DepthOfCoverage -o "+SampleName+"/"+SampleName+".Depth -I "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.clipped.bam --minBaseQuality 22 -baseCounts -ct 0 -ct 1 -ct 250 -ct 1000 -L "+BED_FILE+" -R /Drive_D/Databases/genomes/hg19_genome.fa" + "\n")
    PipelineScript_FILE.write("echo 'Running MakingOfGenotypeGivenAlleles.sh script to create the GenotypeGivenAlleles file , the .Depth.below250X.vcf file and the .DepthForAnnotation File..\n\n'"+"\n")
    if(Description=="FFPE_V1.1" or Description=="FFPE_V1"):
    	PipelineScript_FILE.write("./MakingOfGenotypeGivenAlleles.sh "+SampleName+" /home/bdaber01/Drive_D/Databases/RefBases/hg19.ref."+Description+"_correct"+"\n")
    else:
    	PipelineScript_FILE.write("./MakingOfGenotypeGivenAlleles.sh "+SampleName+" /home/bdaber01/Drive_D/Databases/RefBases/hg19.ref."+Description+"\n")
    PipelineScript_FILE.write("echo 'Intersecting the .Depth.below250X.vcf with the target bed file to get the target regions having reads with depth below 250..\n\n'"+"\n")
    PipelineScript_FILE.write("intersectBed -a "+BED_FILE+"_refseq_intersect.BED -b "+SampleName+"/"+SampleName+".Depth.below250X.vcf -u > "+SampleName+"/"+SampleName+".Depth.below250X.amplicons" + "\n")
    PipelineScript_FILE.write("echo 'Creating Depth.below250X.amplicons.summary file..\n\n'"+"\n")
    PipelineScript_FILE.write("python ampliconSummary.py "+SampleName+"\n")
    PipelineScript_FILE.write("python MakingOfCoverageSummary.py "+SampleName+" "+Description+"\n")
    
    
    
#Phase 2 - Variant Calling and Merging of All Stats

    PipelineScript_FILE.write("echo 'Running GATK Variant Call in Discovery mode..\n\n'"+"\n")
    PipelineScript_FILE.write("java -Xmx12G -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -T UnifiedGenotyper -I "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.clipped.bam --dbsnp /home/bdaber01/Drive_D/Databases/dbsnp/hg19_dbsnp137_correct.vcf -o "+SampleName+"/"+SampleName+".GATK.vcf -nt 6 -glm BOTH -mbq 22 -dt NONE -minIndelCnt 3 -minIndelFrac 0.04 -nda -maxAlleles 5 -stand_emit_conf 10.0 -L "+BED_FILE + "\n")
    PipelineScript_FILE.write("echo 'Running GATK Variant Call in Discovery mode on PRE_FILTERED Target Bam file..\n\n'"+"\n")
    PipelineScript_FILE.write("java -Xmx12G -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -T UnifiedGenotyper -I "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.bam --dbsnp /home/bdaber01/Drive_D/Databases/dbsnp/hg19_dbsnp137_correct.vcf -o "+SampleName+"/"+SampleName+".GATK.pre_filtered.vcf -nt 6 -glm INDEL -mbq 22 -dt NONE -minIndelCnt 3 -minIndelFrac 0.04 -nda -maxAlleles 5 -stand_emit_conf 10.0 -L "+BED_FILE + "\n")
    PipelineScript_FILE.write("echo 'Remove Small Indels (i.e. less than 6bp) from the GATK.pre_filtered.vcf file ..\n\n'"+"\n")
    PipelineScript_FILE.write("python RemoveSmallIndels.py "+SampleName+"/"+SampleName+".GATK.pre_filtered.vcf "+SampleName+"/"+SampleName+".GATK.pre_filtered.greaterThan6bpIndels.vcf "+"\n")
    PipelineScript_FILE.write("echo 'Running GATK Variant Call in GenotypeGivenAlleles mode..\n\n'"+"\n")
    PipelineScript_FILE.write("java -Xmx12G -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -T UnifiedGenotyper -I "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.clipped.bam --dbsnp /home/bdaber01/Drive_D/Databases/dbsnp/hg19_dbsnp137_correct.vcf -o "+SampleName+"/"+SampleName+".GATK.af4percent.vcf -nt 6 -glm SNP -mbq 22 -dt NONE -alleles:VCF "+SampleName+"/"+SampleName+".GenotypeGivenAlleles.sorted.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_emit_conf 10.0 -L "+BED_FILE + "\n")
    PipelineScript_FILE.write("echo 'Combining three GATK Variant Calls to give a combined.vcf file..\n\n'"+"\n")
    PipelineScript_FILE.write("java -Xmx2G -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -R /Drive_D/Databases/genomes/hg19_genome.fa -T CombineVariants --variant "+SampleName+"/"+SampleName+".GATK.vcf --variant "+SampleName+"/"+SampleName+".GATK.af4percent.vcf --variant "+SampleName+"/"+SampleName+".GATK.pre_filtered.greaterThan6bpIndels.vcf -o "+SampleName+"/"+SampleName+".combined.vcf"+"\n")
    PipelineScript_FILE.write("echo 'Filtering (tagging) Variants called at Homopolymer sites..\n\n'"+"\n")
    PipelineScript_FILE.write("java -Xmx2G -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -T VariantFiltration --variant "+SampleName+"/"+SampleName+".combined.vcf -o "+SampleName+"/"+SampleName+".combined.filtered.vcf -dt NONE --filterName 'HomopolymerRun' --filterExpression 'HRun>6'" + "\n")
    PipelineScript_FILE.write("echo 'grepping all variants not called at Homoploymer sites for further annotation..\n\n'"+"\n")
    PipelineScript_FILE.write("grep -v 'HomopolymerRun' "+SampleName+"/"+SampleName+".combined.filtered.vcf > "+SampleName+"/"+SampleName+".combined.filtered.HRR.vcf"+"\n")
    PipelineScript_FILE.write("echo 'Running SnpEff Annotation-- both with canon and without canon option..\n\n'"+"\n")
    PipelineScript_FILE.write("java -jar /home/bdaber01/snpEff_3_0/snpEff.jar -c /home/bdaber01/snpEff_3_0/snpEff.config -canon hg19 -i vcf -o txt -noStats "+SampleName+"/"+SampleName+".combined.filtered.HRR.vcf > "+SampleName+"/"+SampleName+".combined.filtered.HRR.canon_snpeff.txt"+"\n")
    PipelineScript_FILE.write("echo 'Extracting Indels from the canon_snpeff file and calling their depths using their bams..\n\n'"+"\n")
    PipelineScript_FILE.write("./CallIndelDepthAndMerge.sh "+SampleName+"/"+SampleName+"\n")
    PipelineScript_FILE.write("echo 'Running MultiAllelic.py script on the filtered.HRunRemoved.vcf to convert all mutiallelic position to single allelic ones..\n\n'"+"\n")
    PipelineScript_FILE.write("python MultiAllelic.py "+SampleName+"/"+SampleName+".combined.filtered.HRR.vcf "+SampleName+"/"+SampleName+".combined.filtered.HRR.multiallelic.vcf"+"\n")
    PipelineScript_FILE.write("echo 'Running Annovar Refseq Annotation..\n\n'"+"\n")
    PipelineScript_FILE.write("perl /home/bdaber01/annovar/convert2annovar.pl  -format vcf4old -includeinfo "+SampleName+"/"+SampleName+".combined.filtered.HRR.multiallelic.vcf > "+SampleName+"/"+SampleName+".combined.filtered.HRR.multiallelic.avinput"+"\n")
    PipelineScript_FILE.write("perl /home/bdaber01/annovar/annotate_variation.pl "+SampleName+"/"+SampleName+".combined.filtered.HRR.multiallelic.avinput /home/bdaber01/annovar/humandb/ -buildver hg19"+"\n")
    PipelineScript_FILE.write("python refseq.py "+SampleName+"/"+SampleName+".combined.filtered.HRR.multiallelic.avinput.variant_function "+SampleName+"/"+SampleName+".combined.filtered.HRR.multiallelic.avinput.exonic_variant_function "+SampleName+"/"+SampleName+".combined.filtered.HRR.multiallelic.avinput.refseq"+"\n")
    PipelineScript_FILE.write("echo 'Joining All Annotations..\n\n'"+"\n")
    PipelineScript_FILE.write("python AnnotationJoin.py "+SampleName+"/"+SampleName+".combined.filtered.HRR.canon_snpeff.txt "+SampleName+"/"+SampleName+".combined.filtered.HRR.multiallelic.avinput.refseq /Drive_D/Databases/1000Genomes/hg19.1000genomes.extractFields2.annovar_input "+SampleName+"/"+SampleName+".FinalDepthForAnnotation.vcf /Drive_D/Databases/Alalmut/hg19.alalmut "+SampleName+"/"+SampleName+".canon_annotated"+"\n")
    PipelineScript_FILE.write("echo 'Flagging Adjacent Snps..\n\n'"+"\n")
    PipelineScript_FILE.write("python flagAdjacentSNPs.py "+SampleName+"/"+SampleName+".canon_annotated "+SampleName+"/"+SampleName+".canon_annotated_flagged "+"\n")
    PipelineScript_FILE.write("echo 'Running GATK Variant Call in GenotypeGivenAlleles mode to genotype only the signature panel alleles..\n\n'"+"\n")
    PipelineScript_FILE.write("java -Xmx12G -jar /home/bdaber01/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -T UnifiedGenotyper -I "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.clipped.bam --dbsnp /home/bdaber01/Drive_D/Databases/dbsnp/hg19_dbsnp137_correct.vcf -o "+SampleName+"/"+SampleName+".signatureGenotype.vcf -nt 6 -glm SNP -mbq 15 -dt NONE -alleles:VCF /Drive_D/Databases/sigIDPanel/SignatureIDPanel."+Description+".sorted.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_emit_conf 10.0 -L "+BED_FILE + "\n")
    PipelineScript_FILE.write("echo 'Creating the signatureID file for sample..\n\n'"+"\n")
    PipelineScript_FILE.write("python genotypeToSig.py "+SampleName+"/"+SampleName+".signatureGenotype.vcf "+SampleName+"/"+SampleName+".signatureID "+Description+"\n")
    PipelineScript_FILE.write("echo 'Creating Stats file for the sample..\n\n'"+"\n")
    PipelineScript_FILE.write("python SummaryStats.py "+SampleName+"\n")
    PipelineScript_FILE.write('find '+SampleName+'/'+' -type f -name "*.sam" -exec rm -f {} \;'+"\n") # will remove the intermediate sam files to save space
    PipelineScript_FILE.write('find '+SampleName+'/'+' -type f -name "*.fastq" -exec rm -f {} \;'+"\n") # will remove the intermediate fastq files after decompression for btrim. The compressed original files are saved
    PipelineScript_FILE.write('find '+SampleName+'/'+' -type f -name "*.pe" -exec rm -f {} \;'+"\n") # will remove the intermediate .pe files created by btrim
    PipelineScript_FILE.write("mv "+SampleName+"/"+SampleName+".canon_annotated* "+"."+"\n")
    PipelineScript_FILE.write("mv "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.clipped.ba* "+"."+"\n")
    PipelineScript_FILE.write("mv "+SampleName+"/"+SampleName+".novo.sorted.fixed.ba* "+"."+"\n") 
    PipelineScript_FILE.write('find '+SampleName+'/'+' -type f -name "*.bam" -exec rm -f {} \;'+"\n")  # will remove the intermediate Bam files.
    PipelineScript_FILE.write('find '+SampleName+'/'+' -type f -name "*.btrim.summary.txt" -exec rm -f {} \;'+"\n")
    PipelineScript_FILE.close()

SAMPLE_FILE.close()



    
    
    
    
     
    
    
    
    
    
    
    
    
    
    
    
