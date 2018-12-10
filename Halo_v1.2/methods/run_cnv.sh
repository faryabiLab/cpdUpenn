############################################
# CNVKIT test implementation 
# Ashkan Bigdeli 2/4/2016
#http://cnvkit.readthedocs.io/
name=$1
bam=$2
procs=$3
threshold=$4
outdir=$5
#variants=$3
#procs=$4

#dependencies
#cnvkit='/project/cpdlab/Tools/cnvkit/cnvkit.py'
#norm_ref='/project/cpdlab/Scripts/SolidV2/solidv2_v1.2/files/cnv_normref_solid2.cnn'
#target_bed='/project/cpdlab/Scripts/SolidV2/solidv2_v1.2/files/04818-1457701567_Regions_clean.bed'
#append_variants='/project/cpdlab/Scripts/SolidV2/solidv2_v1.2/methods/append_cnv.py'
#make_gbed='/project/cpdlab/Scripts/SolidV2/solidv2_v1.2/methods/cnv_gbed.py'
#java7='/project/cpdlab/Tools/jdk1.7.0_80/bin/java'
#snpeff='/project/cpdlab/Tools/snpEff/4.1l/snpEff.jar'
#snpeff_conf='/project/cpdlab/Tools/snpEff/4.1l/snpEff.config'


##
cnvkit='/project/cpdlab/Tools/cnvkit/cnvkit-0.9.3/cnvkit.py'
norm_ref='/project/cpdlab/ashkan/dev/cnv_validation/reference_pool/sv2_cnv_pooled_ref.8817.cnn'
split_bed='/project/cpdlab/ashkan/dev/cnv_validation/reference_pool/04818-1457701567_Regions_clean.split.target.bed'
MT='/project/cpdlab/ashkan/dev/cnv_validation/reference_pool/MT'


############# AMPLICON METHOD #################

target_cvg=${name}.targetcoverage.cnn
# calculate coverage in bam file
python ${cnvkit} coverage ${bam} ${split_bed} -p ${procs} -o ${target_cvg}

# correct for bias in regional coverage/GC content, output (.cnr) copy number table of ratios
# run segmentation post correction
corrected_cvg=${name}.cnr
python ${cnvkit} fix ${target_cvg} ${MT} ${norm_ref} --no-edge -o ${corrected_cvg}


#perform segmentation for graphing purposes
segments=${name}.cns
python ${cnvkit} segment ${corrected_cvg} -o ${segments}

# generate scatter plot
# This module attempts to access Qt for visulization, it is unecessary
# and will throw an ignorable warning. An update to CNVkit may
# resolve this issue
#export QT_QPA_PLATFORM=offscreen
#scatter_plot=${name}.broad_view.pdf
#python ${cnvkit} scatter ${corrected_cvg} -s ${segments} --title ${name}.Broad -o ${scatter_plot}

#tab_delim outputs ASSUMPTION OF PURE TUMOR for CNV & Segments
cn=${name}.cn.tsv
gainloss=${name}.gainloss.tsv

python ${cnvkit} call ${corrected_cvg} -m clonal --purity 0.55 -o ${cn}
python ${cnvkit} genemetrics -t ${threshold}  ${cn} --stdev --ci -o ${gainloss}

# method can estimate sex, but lets go ahead and take it out since its not
# clinically validated yet
cn_nosex=${name}.cn.nosex.tsv
cn_nosex_named=${name}.cn.nosex.name.tsv
awk '!/chrX|chrY/' ${gainloss} > ${cn_nosex}
awk -v var=${name} -vOFS='' '{print var"\t"$0}' ${cn_nosex} > ${cn_nosex_named}


#for dev purposes
mv ${name}*cn* ${name}*tsv ${name}*e ${name}*o ${out_dir}


#cn_filtered=${name}.cn.filtered.tsv
#awk '$6 >= 4' ${cn_nosex} > ${cn_filtered}

#cn_final=${name}.cn.final.tsv
#awk -v var=${name} -vOFS='' '{print var"\t"$0}' ${cn_filtered} > ${cn_final}

#cut final tab to provide annotations with snpEff
#anno_bed=${name}.cn.anno.bed
#snpeff_anno=${name}.cn.snpeff_anno.tab
#cut -f2,3,4 ${cn_final} | sed 1d > ${anno_bed}
#${java7} -jar ${snpeff} -c ${snpeff_conf} -canon hg19 -i bed -noStats ${anno_bed} > ${snpeff_anno}

#convert to bed file to graph relevant cnvs
#cn_bed=${name}.cn.bed
#python ${make_gbed} ${cn_final} ${cn_bed}
#${cnvkit} scatter ${corrected_cvg} -s ${segments}  -l ${cn_bed} -o ${name}.cnv.pdf 

# append copy number changes to filemaker variant file
#python ${append_variants} ${variants} ${cn_final} ${snpeff_anno} 
