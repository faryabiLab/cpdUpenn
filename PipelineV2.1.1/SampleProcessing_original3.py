##V2.1
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
## Changes in Version 1.3 :   1)Replacing "Samtools sort" with "Novosort".
##                            2)'-Xmx' parameter with value of 2G was mentioned on GATK's CombineVariants command line for the first time just to make sure enough memory is made available to java during the execution of the program.
##                            3)No more generation of non-canon snpeff file.
##                            4)The parameter '-Djava.io.tmpdir=tmp' was included on GATK's DepthOfCoverage, ClipReads and UnifiedGenotyper CommandLines to make sure temp files are directed to the 'tmp' directory under current RunFolder under DriveD.
##                              Previously they were directed to the '/tmp' folder on the linux server which could easily run out of space. DriveD has a lot of space to handle the temp files.
##			      5)'Description' was added to the MakingOfCoverageSummary.py command line.
##			      6) Adjacent Snps are now flagged for each sample. The flags are displayed in the last field(or column) of the cannon_annotated_flagged file.The generated cannon_annotated.flagged is moved to the main folder after analysis.
##			      7) Output File names would now contain 'HRR' instead of 'HRunRemoved'.
##                            8) Annovar software was updated. The value for option -format was changed from 'vcf4' to 'vcf4old' on convert2annovar command-line to make sure that none of the variants were lost.
##                               With the new version of convert2annovar.pl script if we keep the value of -format option as 'vcf4' then it would eliminate the variants having genotype '0/0' and we dont want that
##                               as we dont trust the genotyping by GATK for some low frequency variants.
##			      9) Corrected the hg19.ref.FFPE_V1.1 and hg19.ref.FFPE_V1 files. The new file names are hg19.ref.FFPE_V1.1_correct and hg19.ref.FFPE_V1_correct. An "if" condition has been incorporated in the 
##				  script at "MakingGenotypeGivenAlleles" step to use these corrected files. Correction was done at the chromosome positions chr7:128850461-128850472 to remove the error (i.e. '^M' character)
##                                present in the ref-allele column. This error had no impact on the results reported using the previous pipeline versions as it was found in the intronic region of SMO gene (i.e. chr7:128850461-128850472).
##                            10) The hg19_dbsnp137.vcf is no more used as it had errors at some chromosome positions which had "chrchr" instead of "chr".This error has been rectified and the new file used in the pipe is named as
##                                hg19_dbsnp137_correct.vcf.
##                            11) Added bPPP, mPPP and CEBPa Panels.
##                            12) Added an extra variable named "Library" to process samples based on their library type,i.e. "TSCA (HEME or FFPE)", "PPP (mPPP or bPPP)" or "CEBPA". The "BED_FILE" and "RGLB_string" variables now use "Library" variable instead of "TSCA" string.
##                                Also added "Adapter_R1" and "Adapter_R2" variables to indicate adapter strings
##                                present in the sample reads.
##                            13) CutAdapt tool is now used to cut adapters from fastq reads.
##                            14) This pipe version now includes the steps such as : a) Alignment of the reads with primers-on-but-adapters-off (to call 6bp or larger indels from the resulting bams).This is done for both PPP and TSCA samples.
##                                                                                   b) Cut-Adapted-Primer-Trimmed fastqs are aligned to hg19 ref genome for PPP samples to call SNVs and smaller than 6bp INDELS from resulting bam (post-filter) and larger than 6bp indels from bam (pre-filter).
##                                                                                      For TSCA samples, we still use Non-Cut-Adapted-But-Primmer-Trimmed Fastqs for the same purpose.
##                                                                                      All variants calls are merged asusual into single vcf per sample.
##                                                                                   c) For CEBPa samples the primers-on-but-adapters-off bam is also generated.It is then realigned around indels using GATK's Indel-Realigner. Its name is then changed to "novo.sorted.fixed.bam" to match the regular workflow to call large indels (pre-filter) and small indels and SNVs (post-filter).
##			      15) Added "Library" variable as an argument on SummaryStats.py command line to generate stats appropriately for each library type. SampleCount variable was also added to the command-line so that SummaryStats.py could use fastqs
##			          (containing appropriate sample-count values in their names) for the respective samples.
##                            16) Finding and Removing of the files at the end of processing of each sample is now done by using the sample-specific folder name on the command line instead of asking unix to find the files randomly anywhere in the main run folder.
##                               This modification is done to adopt a safer practice. Btrim summary files are now not removed at the end of pipe. The newly generated primers-on-but-adapters-off bams are also moved from respective sample-folders to the main run folder
##                            17) Added "minIndelCount" and "minIndelFrac" variables to the script to process different library samples separately for variant call. These variables are used on GATK command lines to call Indels.
##                            18) Added "minLenThrow" variable to the script. It allows cutAdapt tool to remove all the sequences having length less than "minLenThrow" after cutting adapters on them.
##			      19) Added "novoalignVer" variable to operate different versions of novoalign for different panels.
##			      20) Shifted zipping of fastqs and cutad fastqs to the end of the script so that SummaryStats.py could utilize the unzipped fastqs to get the total starting read-count. Also changed the novoalign command line (used to generate alignment with primers-on-adapters-off) to intake cutad.fastq instead of cutad.fastq.gz.
## Changes in V1.3.2 :         1) Changed Novoalign version to V3.02.07 and Novosort version to V1.03.01 for both TSCA and PPP.
##                             2) Introduced BED_FILE_PRIMERS_ON target file to operate on PrimersOnBam. These files would help detect insertions (mostly duplications) which spans a part into the target sequence and a part into the primers.
##                             3) Novosort command used to sort the Primers-On Bam contained Novoalign version '3.0.5' on its command line. This was reset to 'novoalignVer' parameter. (Note: novosort version is contained within its respective novocraft"novoalignVer" folder).
## Changes in V1.4 :           1) Added HEMEao_V1.3 panel.
##                             2) OnTarget Bams are now also moved from individual sample folders into the main RunFolder.
## Changes in V1.4.1 : 	       1) Introduced the script ConcordanceExtract.sh to extract the R1-primer and R2-primer concordant reads (only for TSCA and PPP samples) and then carry them forward through the analysis. This would eliminate all the indel artifacts seen in the primer region in the R1-primer-R2-primer-discordant reads.
##			       2) Now TSCA-sample's raw-fastqs also gets processed through the cutAdapt step prior to getting their primers trimmed through Btrim.
##			       3) The value of "minLenThrow" variable has been changed from 1 to 31 for all the TSCA samples.
##		               4) Steps to remove the intermediate files generated by ConcordanceExtract.sh script were added towards the end of this script.
##                             5) Steps to gzip the *.tab.sort.cndt.fastq files were towards the end of this script.
##			       6) Included Time Waits (of 1 minute) between the execution of GATK variant call commands to make sure enough memory space gets available to accomodate and parse heavy read-depth bams at each variant call step.
##		               7) The 'short_limit' parameter used on Btrim command line gets defined for each Panel except for CEBPa for which Btrim is not used.
## Changes in V2.0:	       1) Added BRCA and mBRCA panels.
##	                       2) Changed the if statements to include "Library==BRCA".
##	                       3) Included if,else statements to handle the generation of canon_snpeff specially (with all required transcripts) just for "Library=BRCA".
##		               4) Pre-concordant-filter-primmer-trimmed-fastqs would be used during the Primers-off-alignment step. This would be applicable for only BRCA and mBRCA samples (i.e. Library=="BRCA"). If..else conditions are included at Primers-off-reads-alignment step to account for this change. Samples from all other applicable panels (i.e. HEME_V1.2(or ao_V1.3), FFPE_V1.1 and bPPP(or mPPP)) would go through the usual route, i.e. alignment of concordant primers-off reads to the reference genome.
##			       5) The variants that have HRun>6 are retained if they have allele-frequency>=10%.
##                             6) Novoalign, Novosort, samtools, bedtools, SnpEff and Annovar tools updated to their new versions.
##                             7) For Library=BRCA, fastq.pe files are not deleted at the end of pipe execution and are gzipped.
##			       8) hg19.alamut was changed to hg19.alamut_v1.4.2.Description on AnnotationJoin.py's command line to reflect the usage of new updated Alamut databases. This time the files are separate for separate panels.
##                             9) Command added to generate Depth.below.150X.amplicons file.
##                             10) Echo statement for MakingOfGenotypeGivenAlleles.sh was changed to reflect the generation of Depth.below.150X.vcf file.
##                             11) As the new snpeff version yeilds only vcf output, code-lines (involving both snpEff and snpSift {part of snpEff tool}) have been added to process that output vcf to yeild a similar type of tab-delimited text file (same as before but with two additional columns) so that the downstream joining of annotations from other sources/tools, does not get affected much. This is also done to be concordant with the format of variant entries in CPD's KnowledgeBase.
## 			       12) '.tab' files are now retained at the end of pipe-execution.
##Changes in V2.1:             1) GATK version upgraded to 3.4-46. HomopolymerRun Tagging not supported anymore. Hence, removed codes lines related to GATK's HRun Filter and FilterOutHomopolymerVar.py. Therefore, filenames would no more contain the word "HRunTagged" in them.
##                             2) Added -Djava.io.tmpdir=tmp parameter to each GATK command line that did not have it earlier. The argument "-maxAlleles" changed to "-maxAltAlleles" on respective GATK command lines to suit the new version of GATK. The arguments '-K CPD_upenn_gatk.key' and '-et NO_ET' added to each command line to prevent uploading the GATK-run-stats-report to the AWS server.
##		               3) Changed '-nt 6' to '-nct 6' and '-Xmx12G' to '-Xmx6G' on the command line of GATK's UnifiedGenotyper in order to improve the performance of the process.
##                             4) -genotypeMergeOptions and -priority arguments added to command line of GATK's CombineVariants tool. The version of GATK now requires that these arguments should be mentioned on commandline. Changed '--variant' to '--variant:variant1 (or variant2 or variant3 or variant4) on the respective "CombineVariants" command lines to represent different variant-call-files (i.e vcfs).

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
    	Library="TSCA"
    	Adapter_R1="GCGAATTTCGACGATCGTTGCATTAACTCGCGA"
    	Adapter_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    	minIndelFrac="0.04"
    	minIndelCount="3"
    	minLenThrow=31
    	novoalignVer="3.03.02"
    	short_limit="25"
    elif(Description=="HEME_V1"):
	FragSize="190,50"
	PrimersFile="Primers/PrimersHEME.txt"
    	PrimersRCFile="Primers/Primers_RC_HEME.txt"
    	eParam="190"
    	Library="TSCA"
    	Adapter_R1="GCGAATTTCGACGATCGTTGCATTAACTCGCGA"
    	Adapter_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    	minIndelFrac="0.04"
    	minIndelCount="3"
    	minLenThrow=31
    	novoalignVer="3.03.02"
    	short_limit="25"
    elif(Description=="FFPE_V1.1"):
	PrimersFile="Primers/PrimersFFPE185bp.txt"
	PrimersRCFile="Primers/Primers_RC_FFPE185bp.txt"
	eParam="140"
    	FragSize="140,50"
    	Library="TSCA"
    	Adapter_R1="GCGAATTTCGACGATCGTTGCATTAACTCGCGA"
    	Adapter_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    	minIndelFrac="0.04"
    	minIndelCount="3"
    	minLenThrow=31
    	novoalignVer="3.03.02"
    	short_limit="25"
    elif(Description=="HEME_V1.1"):
    	FragSize="190,50"
    	PrimersFile="Primers/PrimersHEME.txt"
        PrimersRCFile="Primers/Primers_RC_HEME.txt"
    	eParam="190"
    	Library="TSCA"
    	Adapter_R1="GCGAATTTCGACGATCGTTGCATTAACTCGCGA"
    	Adapter_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    	minIndelFrac="0.04"
    	minIndelCount="3"
    	minLenThrow=31
    	novoalignVer="3.03.02"
    	short_limit="25"
    elif(Description=="HEME_V1.2"):
    	FragSize="190,50"
    	PrimersFile="Primers/PrimersHEMEv1.2.txt"
        PrimersRCFile="Primers/Primers_RC_HEMEv1.2.txt"
    	eParam="190"
    	Library="TSCA"
    	Adapter_R1="GCGAATTTCGACGATCGTTGCATTAACTCGCGA"
    	Adapter_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    	minIndelFrac="0.04"
    	minIndelCount="3"
    	minLenThrow=31
    	novoalignVer="3.03.02"
    	short_limit="25"
    elif(Description=="bPPP"):
    	FragSize="150,50"
    	PrimersFile="Primers/PrimersToCutR2_bPPP.txt"
    	PrimersRCFile="Primers/PrimersToCutR1_bPPP.txt"
    	eParam="95"
    	Library="PPP"
    	Adapter_R1="GCGAATTTCGACGATCGTTGCATTAACTCGCGA"
    	Adapter_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    	minIndelFrac="0.04"
    	minIndelCount="3"
    	minLenThrow=51
    	novoalignVer="3.03.02"
    	short_limit="45"
    elif(Description=="mPPP"):
    	FragSize="150,50"
    	PrimersFile="Primers/PrimersToCutR2_mPPP.txt"
    	PrimersRCFile="Primers/PrimersToCutR1_mPPP.txt"
    	eParam="95"
    	Library="PPP"
    	Adapter_R1="GCGAATTTCGACGATCGTTGCATTAACTCGCGA"
    	Adapter_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    	minIndelFrac="0.04"
    	minIndelCount="3"
    	minLenThrow=51
    	novoalignVer="3.03.02"
    	short_limit="45"
    elif(Description=="CEBPa"):
    	FragSize="200,150"
    	Library="CEBPa"
    	Adapter_R2="CTGTCTCTTATACACATCTGACGCTGCCGACG"
    	Adapter_R1="CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
    	minIndelFrac="0.10"
    	minIndelCount="9"
    	minLenThrow=1
    	novoalignVer="3.03.02"
    elif(Description=="HEMEao_V1.3"):
        FragSize="190,50"
        PrimersFile="Primers/PrimersHEMEao_v1.3.txt"
        PrimersRCFile="Primers/Primers_RC_HEMEao_v1.3.txt"
        eParam="190"
        Library="TSCA"
        Adapter_R1="GCGAATTTCGACGATCGTTGCATTAACTCGCGA"
        Adapter_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
        minIndelFrac="0.04"
        minIndelCount="3"
        minLenThrow=31
    	novoalignVer="3.03.02"
    	short_limit="25"
    elif(Description=="BRCA"):
    	FragSize="190,50"
    	PrimersFile="Primers/PrimersToCutR2_BRCA.txt"
    	PrimersRCFile="Primers/PrimersToCutR1_BRCA.txt"
    	eParam="95"
    	Library="BRCA"
    	Adapter_R1="GCGAATTTCGACGATCGTTGCATTAACTCGCGA"
    	Adapter_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    	minIndelFrac="0.04"
    	minIndelCount="3"
    	minLenThrow=50
    	novoalignVer="3.03.02"
    	short_limit="25"
    elif(Description=="mBRCA"):
    	FragSize="190,50"
    	PrimersFile="Primers/PrimersToCutR2_mBRCA.txt"
    	PrimersRCFile="Primers/PrimersToCutR1_mBRCA.txt"
    	eParam="95"
    	Library="BRCA"
    	Adapter_R1="GCGAATTTCGACGATCGTTGCATTAACTCGCGA"
    	Adapter_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    	minIndelFrac="0.04"
    	minIndelCount="3"
    	minLenThrow=50
    	novoalignVer="3.03.02"
    	short_limit="25"
    else:
    	continue
    
    BED_FILE=Library+"_"+"CANCER_PANEL_"+Description+"_target.BED"
    BED_FILE_PRIMERS_ON=Library+"_"+"CANCER_PANEL_"+Description+"_target_PrimersOn.BED"
    RGLB_String=Library+"_"+Description
    SampleCount="S"+ str(Count)
    
    PipeScriptFileName=SampleName+"/"+SampleName+".PipeScript"
    PipelineScript_FILE=open(PipeScriptFileName,'w')
    
# Phase 1 - Mapping, Filtering, Clipping and Stats Generation  
    
    PipelineScript_FILE.write("gunzip "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.fastq.gz"+"\n")
    PipelineScript_FILE.write("gunzip "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.fastq.gz"+"\n")
    
    ## Cutting Adapters from the sequences in fastqs.
    PipelineScript_FILE.write("/home/bdaber01/cutadapt-1.3/bin/cutadapt -e 0.15 -m "+str(minLenThrow)+" -a "+Adapter_R1+" --paired-output "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.cutAd.tmp.fastq"+" -o "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.cutAd.tmp.fastq "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.fastq "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.fastq "+"\n")
    PipelineScript_FILE.write("/home/bdaber01/cutadapt-1.3/bin/cutadapt -e 0.15 -m "+str(minLenThrow)+" -a "+Adapter_R2+" --paired-output "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.cutAd.fastq"+" -o "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.cutAd.fastq "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.cutAd.tmp.fastq "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.cutAd.tmp.fastq "+"\n")
    
    ## Performing Btrim on cut-adapt fastqs for TSCA and PPP samples.
    if(Library=="TSCA" or Library=="PPP" or Library=="BRCA"):
    	PipelineScript_FILE.write("/home/bdaber01/Btrim64 -l "+short_limit+" -p "+PrimersFile+" -t "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.cutAd.fastq -o "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.btrim.fastq -s "+SampleName+"/"+SampleName+"_R2.btrim.summary.txt -w 0 -e "+eParam+" > "+SampleName+"/"+SampleName+"_R2.btrim.stats.txt"+"\n")
    	PipelineScript_FILE.write("/home/bdaber01/Btrim64 -l "+short_limit+" -p "+PrimersRCFile+" -t "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.cutAd.fastq -o "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.btrim.fastq -s "+SampleName+"/"+SampleName+"_R1.btrim.summary.txt -w 0 -e "+eParam+" > "+SampleName+"/"+SampleName+"_R1.btrim.stats.txt"+"\n")
    	PipelineScript_FILE.write("perl paired_end_trim.pl "+SampleName+"/"+SampleName+"_R1.btrim.summary.txt "+SampleName+"/"+SampleName+"_R2.btrim.summary.txt "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.btrim.fastq "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.btrim.fastq"+"\n")
    	PipelineScript_FILE.write("./ConcordanceExtract.sh "+SampleName+" "+SampleCount+"\n")
    	#   Mapping- processing of reads containing primer but not adapter sequences
    	PipelineScript_FILE.write("echo 'Carrying out Alignment of primer-but-not-adapter containing fastq files with hg19 reference genome using Novoalign...\n\n'"+"\n")
    	PipelineScript_FILE.write("/home/bdaber01/novocraftV"+novoalignVer+"/novoalign -d /Drive_D/Databases/genomes/hg19_genome.nix -f "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.cutAd.fastq.tab.sort.cndt.fastq "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.cutAd.fastq.tab.sort.cndt.fastq -i PE "+FragSize+" -c 3 -o FullNW -o SAM "+" > "+SampleName+"/"+SampleName+".novo.primer.sam" +"\n")
    elif(Library=="CEBPa"):
    	PipelineScript_FILE.write("echo 'Carrying out Alignment of primer-but-not-adapter containing fastq files with hg19 reference genome using Novoalign...\n\n'"+"\n")
    	PipelineScript_FILE.write("/home/bdaber01/novocraftV"+novoalignVer+"/novoalign -d /Drive_D/Databases/genomes/hg19_genome.nix -f "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.cutAd.fastq "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.cutAd.fastq -i PE "+FragSize+" -c 3 -o FullNW -o SAM "+" > "+SampleName+"/"+SampleName+".novo.primer.sam" +"\n")
    PipelineScript_FILE.write("echo 'Converting sam to bam using /home/bdaber01/samtools-1.2/samtools...\n\n'"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools view -bS "+SampleName+"/"+SampleName+".novo.primer.sam > "+SampleName+"/"+SampleName+".novo.primer.bam"+"\n")
    PipelineScript_FILE.write("echo 'Sorting initial bam using Novosort...\n\n'"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/novocraftV"+novoalignVer+"/novosort -t tmp -c 3 "+SampleName+"/"+SampleName+".novo.primer.bam > "+SampleName+"/"+SampleName+".novo.primer.sorted.bam"+"\n")
    PipelineScript_FILE.write("echo 'Fixing the read group for sorted bam using Picard Tools...\n\n'"+"\n")
    PipelineScript_FILE.write("java -Xmx2g -jar /home/bdaber01/picard-tools-1.46/AddOrReplaceReadGroups.jar I= "+SampleName+"/"+SampleName+".novo.primer.sorted.bam O= "+SampleName+"/"+SampleName+".novo.primer.sorted.fixed.bam RGID=1 RGLB="+RGLB_String+" RGPL=Illumina RGPU="+Barcode+" RGSM="+SampleName+" RGCN=CPD" + "\n")
    PipelineScript_FILE.write("echo 'Creating Flagstat for bam...\n\n'"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools flagstat "+SampleName+"/"+SampleName+".novo.primer.sorted.fixed.bam > "+SampleName+"/"+SampleName+".novo.primer.sorted.flagstat"+"\n")
    PipelineScript_FILE.write("echo 'Indexing fixed bam using Samtools..\n\n'"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools index "+SampleName+"/"+SampleName+".novo.primer.sorted.fixed.bam" +"\n")
    if(Library=="CEBPa"):
    	PipelineScript_FILE.write("echo 'Performing Realignment on primer.sorted.fixed bam..\n\n'"+"\n")
    	PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -Xmx6G -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -T RealignerTargetCreator -K CPD_upenn_gatk.key -et NO_ET -I "+SampleName+"/"+SampleName+".novo.primer.sorted.fixed.bam -R /Drive_D/Databases/genomes/hg19_genome.fa -o "+SampleName+"/"+SampleName+".forIndelRealigner.Intervals -nt 6 -L "+BED_FILE+ "\n")
    	PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -Xmx2g -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -T IndelRealigner -K CPD_upenn_gatk.key -et NO_ET -I "+SampleName+"/"+SampleName+".novo.primer.sorted.fixed.bam -R /Drive_D/Databases/genomes/hg19_genome.fa -targetIntervals "+SampleName+"/"+SampleName+".forIndelRealigner.Intervals -o "+SampleName+"/"+SampleName+".novo.primer.sorted.fixed.realigned.bam" +"\n")
    	PipelineScript_FILE.write("echo 'renaming novo.primer.sorted.fixed.relaigned.bam to novo.sorted.fixed.bam.. \n\n'"+"\n")
    	PipelineScript_FILE.write("mv "+SampleName+"/"+SampleName+".novo.primer.sorted.fixed.realigned.bam "+SampleName+"/"+SampleName+".novo.sorted.fixed.bam "+"\n")
    elif(Library=="TSCA" or Library=="PPP" or Library=="BRCA"):
    	PipelineScript_FILE.write("echo 'Intersecting target bed with bam...\n\n'"+"\n")
    	PipelineScript_FILE.write("/home/bdaber01/bedtools-2.25.0/bin/intersectBed -abam "+SampleName+"/"+SampleName+".novo.primer.sorted.fixed.bam -b "+BED_FILE_PRIMERS_ON+" > "+SampleName+"/"+SampleName+".novo.primer.sorted.fixed.ontarget.bam" + "\n")
    	PipelineScript_FILE.write("echo 'Indexing on target fixed bam using Samtools..\n\n'"+"\n")
    	PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools index "+SampleName+"/"+SampleName+".novo.primer.sorted.fixed.ontarget.bam" +"\n")
    	PipelineScript_FILE.write("echo 'Creating Flagstat for primer.sorted.fixed.ontarget.bam...\n\n'"+"\n")
    	PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools flagstat "+SampleName+"/"+SampleName+".novo.primer.sorted.fixed.ontarget.bam > "+SampleName+"/"+SampleName+".novo.primer.sorted.fixed.ontarget.flagstat"+"\n")
    	
    	#   Mapping- processing of reads without primer sequences
    	PipelineScript_FILE.write("echo 'Carrying out Alignment of fastq files with hg19 reference genome using Novoalign...\n\n'"+"\n")
    	if(Library=="BRCA"):
    		## Pre-Concordance Filter Fastqs would be used in alignment step for BRCA and mBRCA
    		PipelineScript_FILE.write("/home/bdaber01/novocraftV"+novoalignVer+"/novoalign -d /Drive_D/Databases/genomes/hg19_genome.nix -f "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.btrim.fastq.pe "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.btrim.fastq.pe -i PE "+FragSize+" -c 3 -o FullNW -o SAM > "+SampleName+"/"+SampleName+".novo.sam" +"\n")
    	else:
    		PipelineScript_FILE.write("/home/bdaber01/novocraftV"+novoalignVer+"/novoalign -d /Drive_D/Databases/genomes/hg19_genome.nix -f "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.btrim.fastq.pe.tab.sort.cndt.fastq "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.btrim.fastq.pe.tab.sort.cndt.fastq -i PE "+FragSize+" -c 3 -o FullNW -o SAM > "+SampleName+"/"+SampleName+".novo.sam" +"\n")
    	PipelineScript_FILE.write("echo 'Converting sam to bam using /home/bdaber01/samtools-1.2/samtools...\n\n'"+"\n")
    	PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools view -bS "+SampleName+"/"+SampleName+".novo.sam > "+SampleName+"/"+SampleName+".novo.bam"+"\n")
    	PipelineScript_FILE.write("echo 'Sorting initial bam using Novosort...\n\n'"+"\n")
    	PipelineScript_FILE.write("/home/bdaber01/novocraftV"+novoalignVer+"/novosort -t tmp -c 3 "+SampleName+"/"+SampleName+".novo.bam > "+SampleName+"/"+SampleName+".novo.sorted.bam"+"\n")
    	PipelineScript_FILE.write("echo 'Fixing the read group for sorted bam using Picard Tools...\n\n'"+"\n")
    	PipelineScript_FILE.write("java -Xmx2g -jar /home/bdaber01/picard-tools-1.46/AddOrReplaceReadGroups.jar I= "+SampleName+"/"+SampleName+".novo.sorted.bam O= "+SampleName+"/"+SampleName+".novo.sorted.fixed.bam RGID=1 RGLB="+RGLB_String+" RGPL=Illumina RGPU="+Barcode+" RGSM="+SampleName+" RGCN=CPD" + "\n")
    PipelineScript_FILE.write("echo 'Creating Flagstat for bam...\n\n'"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools flagstat "+SampleName+"/"+SampleName+".novo.sorted.fixed.bam > "+SampleName+"/"+SampleName+".novo.sorted.flagstat"+"\n")
    PipelineScript_FILE.write("echo 'Indexing fixed bam using Samtools..\n\n'"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools index "+SampleName+"/"+SampleName+".novo.sorted.fixed.bam" +"\n")
    PipelineScript_FILE.write("echo 'Intersecting target bed with bam...\n\n'"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/bedtools-2.25.0/bin/intersectBed -abam "+SampleName+"/"+SampleName+".novo.sorted.fixed.bam -b "+BED_FILE+" > "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.bam" + "\n")
    PipelineScript_FILE.write("echo 'Indexing on target fixed bam using Samtools..\n\n'"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools index "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.bam" +"\n")
    PipelineScript_FILE.write("echo 'Creating Flagstat for sorted.fixed.ontarget.bam...\n\n'"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools flagstat "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.bam > "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.flagstat"+"\n")
    PipelineScript_FILE.write("echo 'Samtools filtering of ontarget bam using minimum mapping quality cutoff of 40 for reads...\n\n'"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools view -h -q 40 "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.bam > "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.q40.sam"+"\n")
    PipelineScript_FILE.write("echo 'Running FilterAlignScore.py to further filter the reads using minimum alignment score cutoff of 95...\n\n'"+"\n")
    PipelineScript_FILE.write("python FilterAlignScore.py "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.q40.sam "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.q40.as95.sam" + "\n")
    PipelineScript_FILE.write("echo 'Converting filtered sam to bam using Samtools...\n\n'"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools view -bS "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.q40.as95.sam > "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.bam" + "\n")
    PipelineScript_FILE.write("echo 'Creating Flagstat for filtered bam...\n\n'"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools flagstat "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.bam > "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.flagstat" + "\n")
    PipelineScript_FILE.write("echo 'Indexing filtered bam using Samtools..\n\n'"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools index "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.bam" +"\n")
    PipelineScript_FILE.write("echo 'Soft Clipping the end of reads...\n\n'"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -Xmx2g -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -T ClipReads -K CPD_upenn_gatk.key -et NO_ET -I "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.bam  -o "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.clipped.bam -os "+SampleName+"/"+SampleName+".clippingstats.txt -R /Drive_D/Databases/genomes/hg19_genome.fa -CR SOFTCLIP_BASES -QT 22"+"\n")
    PipelineScript_FILE.write("echo 'Indexing clipped bam using Samtools..\n\n'"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/samtools-1.2/samtools index "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.clipped.bam" +"\n")
    PipelineScript_FILE.write("echo 'Running GATK DepthOfCoverage..\n\n'"+"\n")  
    PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -Xmx2g -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -T DepthOfCoverage -K CPD_upenn_gatk.key -et NO_ET -o "+SampleName+"/"+SampleName+".Depth -I "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.clipped.bam --minBaseQuality 22 -baseCounts -ct 0 -ct 1 -ct 250 -ct 1000 -L "+BED_FILE+" -R /Drive_D/Databases/genomes/hg19_genome.fa" + "\n")
    PipelineScript_FILE.write("echo 'Running MakingOfGenotypeGivenAlleles.sh script to create the GenotypeGivenAlleles file , the .Depth.below250X.vcf file, .Depth.below150X.vcf file and the .DepthForAnnotation File..\n\n'"+"\n")
    if(Description=="FFPE_V1.1" or Description=="FFPE_V1"):
    	PipelineScript_FILE.write("./MakingOfGenotypeGivenAlleles.sh "+SampleName+" /home/bdaber01/Drive_D/Databases/RefBases/hg19.ref."+Description+"_correct"+"\n")
    else:
    	PipelineScript_FILE.write("./MakingOfGenotypeGivenAlleles.sh "+SampleName+" /home/bdaber01/Drive_D/Databases/RefBases/hg19.ref."+Description+"\n")
    PipelineScript_FILE.write("echo 'Intersecting the .Depth.below250X.vcf with the target bed file to get the target regions having reads with depth below 250..\n\n'"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/bedtools-2.25.0/bin/intersectBed -a "+BED_FILE+"_refseq_intersect.BED -b "+SampleName+"/"+SampleName+".Depth.below250X.vcf -u > "+SampleName+"/"+SampleName+".Depth.below250X.amplicons" + "\n")
    PipelineScript_FILE.write("/home/bdaber01/bedtools-2.25.0/bin/intersectBed -a "+BED_FILE+"_refseq_intersect.BED -b "+SampleName+"/"+SampleName+".Depth.below150X.vcf -u > "+SampleName+"/"+SampleName+".Depth.below150X.amplicons" + "\n")
    PipelineScript_FILE.write("echo 'Creating Depth.below250X(and 150X).amplicons.summary file..\n\n'"+"\n")
    PipelineScript_FILE.write("python ampliconSummary.py "+SampleName+"\n")
    PipelineScript_FILE.write("python MakingOfCoverageSummary.py "+SampleName+" "+Description+"\n")
    
    
    
#Phase 2 - Variant Calling and Merging of All Stats

    PipelineScript_FILE.write("echo 'Running GATK Variant Call in Discovery mode..\n\n'"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -Xmx6G -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -K CPD_upenn_gatk.key -et NO_ET -T UnifiedGenotyper -I "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.clipped.bam --dbsnp /home/bdaber01/Drive_D/Databases/dbsnp/hg19_dbsnp137_correct.vcf -o "+SampleName+"/"+SampleName+".GATK.vcf -nct 1 -glm BOTH -mbq 22 -dt NONE -minIndelCnt "+minIndelCount+" -minIndelFrac "+minIndelFrac+" -nda -maxAltAlleles 5 -stand_emit_conf 10.0 -L "+BED_FILE + "\n")
    
    PipelineScript_FILE.write("echo 'waiting for 1 minute before executing next command'"+"\n")
    PipelineScript_FILE.write("sleep 1m"+"\n")
    # calling 6bp or larger on a pre-filtered bam file not containing primers for TSCA and PPP samples and on a pre-filtered realgined bam for CEBPa samples
    PipelineScript_FILE.write("echo 'Running GATK Variant Call in Discovery mode on PRE_FILTERED Target Bam file..\n\n'"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -Xmx6G -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -K CPD_upenn_gatk.key -et NO_ET -T UnifiedGenotyper -I "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.bam --dbsnp /home/bdaber01/Drive_D/Databases/dbsnp/hg19_dbsnp137_correct.vcf -o "+SampleName+"/"+SampleName+".GATK.pre_filtered.vcf -nct 1 -glm INDEL -mbq 22 -dt NONE -minIndelCnt "+minIndelCount+" -minIndelFrac "+minIndelFrac+" -nda -maxAltAlleles 5 -stand_emit_conf 10.0 -L "+BED_FILE + "\n")
    PipelineScript_FILE.write("echo 'Remove Small Indels (i.e. less than 6bp) from the GATK.pre_filtered.vcf file ..\n\n'"+"\n")
    PipelineScript_FILE.write("python RemoveSmallIndels.py "+SampleName+"/"+SampleName+".GATK.pre_filtered.vcf "+SampleName+"/"+SampleName+".GATK.pre_filtered.greaterThan6bpIndels.vcf "+"\n")
    if(Library=="TSCA" or Library=="PPP" or Library=="BRCA"):
    	PipelineScript_FILE.write("echo 'waiting for 1 minute before executing next command'"+"\n")
    	PipelineScript_FILE.write("sleep 1m"+"\n")
    	# calling 6bp or larger on a pre-filtered bam file containing primers
    	PipelineScript_FILE.write("echo 'Running GATK Variant Call in Discovery mode on PRIMERS-ON-PRE_FILTERED Target Bam file..\n\n'"+"\n")
    	PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -Xmx6G -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -K CPD_upenn_gatk.key -et NO_ET -T UnifiedGenotyper -I "+SampleName+"/"+SampleName+".novo.primer.sorted.fixed.ontarget.bam --dbsnp /home/bdaber01/Drive_D/Databases/dbsnp/hg19_dbsnp137_correct.vcf -o "+SampleName+"/"+SampleName+".GATK.primer.pre_filtered.vcf -nct 1 -glm INDEL -mbq 22 -dt NONE -minIndelCnt "+minIndelCount+" -minIndelFrac "+minIndelFrac+" -nda -maxAltAlleles 5 -stand_emit_conf 10.0 -L "+BED_FILE_PRIMERS_ON + "\n")
    	PipelineScript_FILE.write("echo 'Remove Small Indels (i.e. less than 6bp) from the GATK.primer.pre_filtered.vcf file ..\n\n'"+"\n")
    	PipelineScript_FILE.write("python RemoveSmallIndels.py "+SampleName+"/"+SampleName+".GATK.primer.pre_filtered.vcf "+SampleName+"/"+SampleName+".GATK.primer.pre_filtered.greaterThan6bpIndels.vcf "+"\n")
    
    PipelineScript_FILE.write("echo 'waiting for 1 minute before executing next command'"+"\n")
    PipelineScript_FILE.write("sleep 1m"+"\n")
    # calling SNVs in GenotypeGivenAlleles mode
    
    PipelineScript_FILE.write("echo 'Running GATK Variant Call in GenotypeGivenAlleles mode..\n\n'"+"\n")
    PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -Xmx6G -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -K CPD_upenn_gatk.key -et NO_ET -T UnifiedGenotyper -I "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.clipped.bam --dbsnp /home/bdaber01/Drive_D/Databases/dbsnp/hg19_dbsnp137_correct.vcf -o "+SampleName+"/"+SampleName+".GATK.af4percent.vcf -nct 1 -glm SNP -mbq 22 -dt NONE -alleles:VCF "+SampleName+"/"+SampleName+".GenotypeGivenAlleles.sorted.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_emit_conf 10.0 -L "+BED_FILE + "\n")
    
    # combining all variant calls together
    if(Library=="TSCA" or Library=="PPP" or Library=="BRCA"):
    	PipelineScript_FILE.write("echo 'Combining all four GATK Variant Calls to give a combined.vcf file..\n\n'"+"\n")
    	PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -Xmx2G -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -R /Drive_D/Databases/genomes/hg19_genome.fa -T CombineVariants -K CPD_upenn_gatk.key -et NO_ET --variant:variant1 "+SampleName+"/"+SampleName+".GATK.vcf --variant:variant2 "+SampleName+"/"+SampleName+".GATK.af4percent.vcf --variant:variant3 "+SampleName+"/"+SampleName+".GATK.pre_filtered.greaterThan6bpIndels.vcf --variant:variant4 "+SampleName+"/"+SampleName+".GATK.primer.pre_filtered.greaterThan6bpIndels.vcf -genotypeMergeOptions PRIORITIZE -priority variant3,variant4,variant2,variant1 -o "+SampleName+"/"+SampleName+".combined.vcf"+"\n")
    elif(Library=="CEBPa"):
    	PipelineScript_FILE.write("echo 'Combining all four GATK Variant Calls to give a combined.vcf file..\n\n'"+"\n")
    	PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -Xmx2G -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -R /Drive_D/Databases/genomes/hg19_genome.fa -T CombineVariants -K CPD_upenn_gatk.key -et NO_ET --variant:variant1 "+SampleName+"/"+SampleName+".GATK.vcf --variant:variant2 "+SampleName+"/"+SampleName+".GATK.af4percent.vcf --variant:variant3 "+SampleName+"/"+SampleName+".GATK.pre_filtered.greaterThan6bpIndels.vcf -genotypeMergeOptions PRIORITIZE -priority variant3,variant2,variant1 -o "+SampleName+"/"+SampleName+".combined.vcf"+"\n")
    
    # Homopolymer Tagging of the variants
    
    #PipelineScript_FILE.write("echo 'Filtering (tagging) Variants called at Homopolymer sites..\n\n'"+"\n")
    #PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -Xmx2G -Djava.io.tmpdir=tmp -jar /home/bdaber01/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -R /home/bdaber01/Drive_D/Databases/genomes/hg19_genome.fa -K CPD_upenn_gatk.key -et NO_ET -T VariantFiltration --variant "+SampleName+"/"+SampleName+".combined.vcf -o "+SampleName+"/"+SampleName+".combined.vcf -dt NONE --filterName 'HomopolymerRun' --filterExpression 'HRun>6'" + "\n")
    #PipelineScript_FILE.write("echo 'grepping all variants not called at Homoploymer sites for further annotation..\n\n'"+"\n")
    #PipelineScript_FILE.write("grep -v 'HomopolymerRun' "+SampleName+"/"+SampleName+".combined.filtered.vcf > "+SampleName+"/"+SampleName+".combined.vcf"+"\n")
    
    # Annotating the variants
    
    PipelineScript_FILE.write("echo 'Running SnpEff Annotation--with canon option activated (BRCA/mBRCA is an exception where canon option is not activated initially, but the unimportant transcripts are filtered out in a later step..\n\n'"+"\n")
    if(Library=="BRCA"):
	PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -jar /home/bdaber01/snpEff_v_4.1/snpEff.jar -c /home/bdaber01/snpEff_v_4.1/snpEff.config hg19 -i vcf -o vcf -noStats "+SampleName+"/"+SampleName+".combined.vcf > "+SampleName+"/"+SampleName+".combined.snpeff.vcf"+"\n")
	PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -jar /home/bdaber01/snpEff_v_4.1/SnpSift.jar varType "+SampleName+"/"+SampleName+".combined.snpeff.vcf > "+SampleName+"/"+SampleName+".combined.snpeff.vT.vcf"+"\n")
	PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -jar /home/bdaber01/snpEff_v_4.1/SnpSift.jar extractFields "+SampleName+"/"+SampleName+".combined.snpeff.vT.vcf CHROM POS ID REF ALT QUAL FILTER DP ANN HET HOM VARTYPE > "+SampleName+"/"+SampleName+".combined.snpeff.vT.ex.txt"+"\n")
	PipelineScript_FILE.write("python PrepareSnpEffOut.py "+SampleName+"/"+SampleName+".combined.snpeff.vT.ex.txt "+SampleName+"/"+SampleName+".combined.snpeff.vT.ex.tab"+"\n") 
    	PipelineScript_FILE.write("egrep -v 'NM_001291241|NM_001122741|NM_001291230|NM_000125|NM_007297|NM_007298|NM_007299|NR_027676|NM_032458|NM_032335|NM_001015877|NM_003140' "+SampleName+"/"+SampleName+".combined.snpeff.vT.ex.tab > "+SampleName+"/"+SampleName+".combined.canon_snpeff.vT.ex.tab"+"\n")
    else:
    	PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -jar /home/bdaber01/snpEff_v_4.1/snpEff.jar -c /home/bdaber01/snpEff_v_4.1/snpEff.config -canon hg19 -i vcf -o vcf -noStats "+SampleName+"/"+SampleName+".combined.vcf > "+SampleName+"/"+SampleName+".combined.canon_snpeff.vcf"+"\n")
	PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -jar /home/bdaber01/snpEff_v_4.1/SnpSift.jar varType "+SampleName+"/"+SampleName+".combined.canon_snpeff.vcf > "+SampleName+"/"+SampleName+".combined.canon_snpeff.vT.vcf"+"\n")
	PipelineScript_FILE.write("/home/bdaber01/jdk1.7.0_80/bin/java -jar /home/bdaber01/snpEff_v_4.1/SnpSift.jar extractFields "+SampleName+"/"+SampleName+".combined.canon_snpeff.vT.vcf CHROM POS ID REF ALT QUAL FILTER DP ANN HET HOM VARTYPE >"+SampleName+"/"+SampleName+".combined.canon_snpeff.vT.ex.txt"+"\n")
	PipelineScript_FILE.write("python PrepareSnpEffOut.py "+SampleName+"/"+SampleName+".combined.canon_snpeff.vT.ex.txt "+SampleName+"/"+SampleName+".combined.canon_snpeff.vT.ex.tab"+"\n")
    PipelineScript_FILE.write("echo 'Running MultiAllelic.py script on the filtered.HRunRemoved.vcf to convert all mutiallelic position to single allelic ones..\n\n'"+"\n")
    PipelineScript_FILE.write("python MultiAllelic.py "+SampleName+"/"+SampleName+".combined.vcf "+SampleName+"/"+SampleName+".combined.multiallelic.vcf"+"\n")
    PipelineScript_FILE.write("echo 'Running Annovar Refseq Annotation..\n\n'"+"\n")
    PipelineScript_FILE.write("perl /home/bdaber01/annovar_17thJun2015/convert2annovar.pl  -format vcf4old -includeinfo "+SampleName+"/"+SampleName+".combined.multiallelic.vcf > "+SampleName+"/"+SampleName+".combined.multiallelic.avinput"+"\n")
    PipelineScript_FILE.write("perl /home/bdaber01/annovar_17thJun2015/annotate_variation.pl "+SampleName+"/"+SampleName+".combined.multiallelic.avinput /home/bdaber01/annovar_17thJun2015/humandb/ -buildver hg19"+"\n")
    PipelineScript_FILE.write("python refseq.py "+SampleName+"/"+SampleName+".combined.multiallelic.avinput.variant_function "+SampleName+"/"+SampleName+".combined.multiallelic.avinput.exonic_variant_function "+SampleName+"/"+SampleName+".combined.multiallelic.avinput.refseq"+"\n")
    PipelineScript_FILE.write("echo 'Joining All Annotations..\n\n'"+"\n")
    PipelineScript_FILE.write("python AnnotationJoin.py "+SampleName+"/"+SampleName+".combined.canon_snpeff.vT.ex.tab "+SampleName+"/"+SampleName+".combined.multiallelic.avinput.refseq /Drive_D/Databases/1000Genomes/hg19.1000genomes.extractFields2.annovar_input "+SampleName+"/"+SampleName+".DepthForAnnotation.sorted.vcf /Drive_D/Databases/Alalmut/hg19.alamut_v1.4.2."+Description+" "+SampleName+"/"+SampleName+".canon_annotated"+"\n")
    PipelineScript_FILE.write("echo 'Flagging Adjacent Snps..\n\n'"+"\n")
    PipelineScript_FILE.write("python flagAdjacentSNPs.py "+SampleName+"/"+SampleName+".canon_annotated "+SampleName+"/"+SampleName+".canon_annotated_flagged "+"\n")
    #PipelineScript_FILE.write("python FilterOutHomopolymerVar.py "+SampleName+"/"+SampleName+".canon_annotated_flagged "+SampleName+"/"+SampleName+".canon_annotated_flagged.HRR"+"\n")
    
    # Creating SummaryStats File
    
    PipelineScript_FILE.write("echo 'Creating Stats file for the sample..\n\n'"+"\n")
    PipelineScript_FILE.write("python SummaryStats.py "+SampleName+" "+Library+" "+SampleCount+"\n")
    PipelineScript_FILE.write("gzip "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.fastq"+"\n")
    PipelineScript_FILE.write("gzip "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.fastq"+"\n")
    PipelineScript_FILE.write("gzip "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.cutAd.fastq"+"\n")
    PipelineScript_FILE.write("gzip "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.cutAd.fastq"+"\n")
    PipelineScript_FILE.write("gzip "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.cutAd.fastq.tab.sort.cndt.fastq"+"\n")
    PipelineScript_FILE.write("gzip "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.cutAd.fastq.tab.sort.cndt.fastq"+"\n")
    PipelineScript_FILE.write("gzip "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R1_001.btrim.fastq.pe.tab.sort.cndt.fastq"+"\n")
    PipelineScript_FILE.write("gzip "+SampleName+"/"+SampleName+"_"+SampleCount+"_L001_R2_001.btrim.fastq.pe.tab.sort.cndt.fastq"+"\n")
    PipelineScript_FILE.write('find '+SampleName+'/'+' -type f -name "*.sam" -exec rm -f {} \;'+"\n") # will remove the intermediate sam files to save space
    PipelineScript_FILE.write('find '+SampleName+'/'+' -type f -name "*.fastq" -exec rm -f {} \;'+"\n") # will remove the intermediate fastq files after decompression for btrim. The compressed original files are saved
    if(Library!="BRCA"):
    	PipelineScript_FILE.write('find '+SampleName+'/'+' -type f -name "*.pe" -exec rm -f {} \;'+"\n") # will remove the intermediate .pe files created by btrim
    else:
    	PipelineScript_FILE.write("gzip "+SampleName+"/"+SampleName+"_"+SampleCount+"_*.pe"+"\n")
    PipelineScript_FILE.write('find '+SampleName+'/'+' -type f -name "*.se" -exec rm -f {} \;'+"\n")
    #PipelineScript_FILE.write('find '+SampleName+'/'+' -type f -name "*.tab" -exec rm -f {} \;'+"\n")
    PipelineScript_FILE.write('find '+SampleName+'/'+' -type f -name "*.cndt" -exec rm -f {} \;'+"\n")
    PipelineScript_FILE.write('find '+SampleName+'/'+' -type f -name "*.sort" -exec rm -f {} \;'+"\n")
    PipelineScript_FILE.write('find '+SampleName+'/'+' -type f -name "*.sort.txt" -exec rm -f {} \;'+"\n")
    PipelineScript_FILE.write("mv "+SampleName+"/"+SampleName+".canon_annotated* "+"."+"\n")
    PipelineScript_FILE.write("mv "+SampleName+"/"+SampleName+".novo.primer.sorted.fixed.ba* "+"."+"\n")
    PipelineScript_FILE.write("mv "+SampleName+"/"+SampleName+".novo.sorted.fixed.ontarget.afterFilter.clipped.ba* "+"."+"\n")
    PipelineScript_FILE.write("mv "+SampleName+"/"+SampleName+".novo.sorted.fixed.ba* "+"."+"\n")
    PipelineScript_FILE.write("mv "+SampleName+"/"+SampleName+".*.ontarget.ba* "+"."+"\n")
    PipelineScript_FILE.write('find '+SampleName+'/'+' -type f -name "*.bam" -exec rm -f {} \;'+"\n")  # will remove the intermediate Bam files.
    #PipelineScript_FILE.write('find '+SampleName+'/'+' -type f -name "*.btrim.summary.txt" -exec rm -f {} \;'+"\n")
    PipelineScript_FILE.close()

SAMPLE_FILE.close()



    
    
    
    
     
    
    
    
    
    
    
    
    
    
    
    
