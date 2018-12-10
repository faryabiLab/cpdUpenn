variants=$2_Run_masterVarFinal.txt
stats=$2_RunStatsFinal.txt
amps=$2_Amplifications.txt
cov=$2_MeanCovByPostion.txt

echo "SampleName	Chrom	Pos	Ref	Alt	Variant Type(SnpEff)	Homozygous	Quality	Coverage	Warnings	Gene_ID	Gene_name	Bio_type	Trancript_ID	Exon_ID	Exon_Rank	Effect	old_AA/new_AA	Old_codon/New_codon	Codon_Num(CDS)	Codon_Degeneracy	CDS_size	Codons_around	AAs_around	Custom_interval_ID	c_Change(SnpEff)	p_Change(SnpEff)	Location(annovar)	Gene	Consequence	Transcript	cDNA change	Protien Change	CHROM(vcf)	POS	dbsnp137_ID	Ref	Alt	QUAL	FILTER	INFO	GenotypeFormat	Genotype	AC(1000g)	AN(1000g)	AF(1000g)	AMR_AF	ASN_AF	EUR_AF	AFR_AF	VariantType	FDP	FRD	FAD	FAF	gene	geneId	transcript	strand	transLen	protein	Uniprot	varType	codingEffect	varLocation	assembly	gDNAstart	gDNAend	gNomen	cDNAstart	cDNAend	cNomen	pNomen	alt_pNomen	exon	intron	omimId	distNearestSS	nearestSSType	wtSSFScore	wtMaxEntScore	wtNNSScore	wtGSScore	wtHSFScore	varSSFScore	varMaxEntScore	varNNSScore	varGSScore	varHSFScore	nearestSSChange	localSpliceEffect	proteinDomain1	proteinDomain2	proteinDomain3	proteinDomain4	rsId	rsValidated	rsSuspect	rsValidations	rsValidationNumber	rsAncestralAllele	rsHeterozygosity	rsClinicalSignificance	rsMAF	rsMAFAllele	rsMAFCount	1000g_AF	1000g_AFR_AF	1000g_SAS_AF	1000g_EAS_AF	1000g_EUR_AF	1000g_AMR_AF	exacAllFreq	exacAFRFreq	exacAMRFreq	exacEASFreq	exacSASFreq	exacNFEFreq	exacFINFreq	exacOTHFreq	exacAFRHmz	exacAMRHmz	exacEASHmz	exacSASHmz	exacNFEHmz	exacFINHmz	exacOTHHmz	exacFilter	exacReadDepth	espRefEACount	espRefAACount	espRefAllCount	espAltEACount	espAltAACount	espAltAllCount	espEAMAF	espAAMAF	espAllMAF	espEAAAF	espAAAAF	espAllAAF	espAvgReadDepth	clinVarIds	clinVarOrigins	clinVarMethods	clinVarClinSignifs	clinVarReviewStatus	clinVarPhenotypes	hgmdId	hgmdPhenotype	hgmdPubMedId	hgmdSubCategory	cosmicIds	cosmicTissues	cosmicFreqs	cosmicSampleCounts	insNucs	delNucs	substType	wtNuc	varNuc	nucChange	phastCons	phyloP	wtAA_1	wtAA_3	wtCodon	wtCodonFreq	varAA_1	varAA_3	varCodon	varCodonFreq	posAA	nOrthos	conservedOrthos	conservedDistSpecies	BLOSUM45	BLOSUM62	BLOSUM80	wtAAcomposition	varAAcomposition	wtAApolarity	varAApolarity	wtAAvolume	varAAvolume	granthamDist	AGVGDclass	AGVGDgv	AGVGDgd	SIFTprediction	SIFTweight	SIFTmedian	MAPPprediction	MAPPpValue	MAPPpValueMedian	AdjacentSNPflag" > $variants

find ./$1 -name "*variants.upload" -type f -exec cat {} \; >> $variants

echo "SampleName	TotalStartingReads	TotalReadsInputForAlignment	PercantageLost	ReadsMapped	PercentageReadMapping	TotalOnTarget	PercentageOnTarget	TotalOnTargetFilterReads	PercentageOnTargetFilter	PercentageUsable	MeanCoverage	PercentBases_above_0	PercentBases_above_1	PercentBases_above_250	PercentBases_above_1000	No._of_amplicons(exonic_part)_below_250x	ClipCount	No._of_amplicons(exonic_part)_below_150x\n" > $stats

find ./$1 -name "*stats.upload" -type f -exec cat {} \; >> $stats


echo "SampleName	Gene	Amp" > ${amps}

find ./$1 -name "*amps.upload" -type f -exec cat {} \; >> $amps


find ./$1 -name "*depths.upload" -type f -exec cat {} \; >> $cov
