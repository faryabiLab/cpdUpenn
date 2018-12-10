############################################
# CNVKIT test implementation 
# Ashkan Bigdeli 2/4/2016
#http://cnvkit.readthedocs.io/

out_dir=$1
cd ${out_dir}
name=$2
bam=$3
#variants=$3
#procs=$4

#dependencies
cnvkit='/project/cpdlab/Tools/cnvkit/cnvkit.py'
norm_ref='/project/cpdlab/Scripts/SolidV2/solidv2_v1.2/files/cnv_normref_solid2.cnn'
target_bed='/project/cpdlab/Scripts/SolidV2/solidv2_v1.2/files/04818-1457701567_Regions_clean.bed'
append_variants='/project/cpdlab/Scripts/SolidV2/solidv2_v1.2/methods/append_cnv.py'
make_gbed='/project/cpdlab/Scripts/SolidV2/solidv2_v1.2/methods/cnv_gbed.py'
java7='/project/cpdlab/Tools/jdk1.7.0_80/bin/java'
snpeff='/project/cpdlab/Tools/snpEff/4.1l/snpEff.jar'
snpeff_conf='/project/cpdlab/Tools/snpEff/4.1l/snpEff.config'

############# AMPLICON METHOD #################


split_bed=${name}.split.bed
#analyze target bed

python ${cnvkit} target ${target_bed} --split -o ${split_bed}
touch MT

pool_target_cvg=${name}.targetcoverage.cnn
# calculate coverage in bam file
python ${cnvkit} coverage ${bam} ${split_bed} -p ${procs} -o ${pool_target_cvg}

# correct for bias in regional coverage/GC content, output (.cnr) copy number table of ratios
# run segmentation post correction
corrected_cvg=${name}.cnr
python ${cnvkit} fix ${pool_target_cvg} MT ${norm_ref} -o ${corrected_cvg}


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
gainloss=${name}.gainloss.tab
cn=${name}.cn.tab
python ${cnvkit} gainloss -t 0.585 ${corrected_cvg} -o ${gainloss}
python ${cnvkit} call ${gainloss} -o ${cn}

# method can estimate sex, but lets go ahead and take it out since its not
# clinically validated yet
cn_nosex=${name}.cn.nosex.tab
awk '!/chrX|chrY/' ${cn} > ${cn_nosex}

cn_filtered=${name}.cn.filtered.tsv
awk '$6 >= 4' ${cn_nosex} > ${cn_filtered}

cn_final=${name}.cn.final.tab
awk -v var=${name} -vOFS='' '{print var"\t"$0}' ${cn_filtered} > ${cn_final}

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
