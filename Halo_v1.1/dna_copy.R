#Ashkan Bigdeli
#
# Amps.R takes in the log ratio's of sample reads and segments to 
# determine copy number. 
library(DNAcopy)
library(dplyr)
library(plotrix)
rm(list = ls())
args = commandArgs(trailingOnly=TRUE)

#sample_name="CPDC163141_Solid_16326_B02_UMI_PAL_16234_SEQ_160264"
sample_name=args[1]
sample_ratio=args[2]
outfile=args[3]

# run DNAcopy - https://bioconductor.org/packages/release/bioc/manuals/DNAcopy/man/DNAcopy.pdf 
amp_df = read.table(sample_ratio, sep="\t", header=T, stringsAsFactors = F)
CNA.object = CNA(cbind(amp_df$Amp), amp_df$Chrom, amp_df$Pos, data.type="logratio",sampleid=sample_name, presorted=T)
smoothed.CNA.object = smooth.CNA(CNA.object)
segment.smoothed.CNA.object = segment(smoothed.CNA.object, verbose=1)

#extract information for plotting segments from CNA object
seg_idx = as.data.frame(segment.smoothed.CNA.object$output)
seg_idx$startRow = segment.smoothed.CNA.object$segRows$startRow
seg_idx$endRow = segment.smoothed.CNA.object$segRows$endRow

#generate dataframe of putative positive results (2 fold log2 change)
seg_idx_filt_pos = filter(seg_idx, seg.mean > 1.0 )
#seg_idx_filt_pos$cnv_pos = paste(seg_idx_filt_pos$chrom, ":", seg_idx_filt_pos$loc.start, "-", seg_idx_filt_pos$loc.end, "(", seg_idx_filt_pos$seg.mean, ")")
seg_idx_filt_pos$cnv_pos = paste(seg_idx_filt_pos$chrom, seg_idx_filt_pos$loc.start, sep=":")
seg_idx_filt_pos$cnv_pos = paste(seg_idx_filt_pos$cnv_pos, seg_idx_filt_pos$loc.end, sep="-")
#seg_idx_filt_pos$cnv_pos = paste(seg_idx_filt_pos$cnv_pos, seg_idx_filt_pos$seg.mean, sep="=")

#generate dataframe of putative negative results (2 fold log2 change)
seg_idx_filt_neg = filter(seg_idx, seg.mean < -1.0 )
#seg_idx_filt_neg$cnv_pos = paste(seg_idx_filt_neg$chrom, ":", seg_idx_filt_neg$loc.start, "-", seg_idx_filt_neg$loc.end, "(", seg_idx_filt_neg$seg.mean, ")")
seg_idx_filt_neg$cnv_pos = paste(seg_idx_filt_neg$chrom, seg_idx_filt_neg$loc.start, sep=":")
seg_idx_filt_neg$cnv_pos = paste(seg_idx_filt_neg$cnv_pos, seg_idx_filt_neg$loc.end, sep="\U002D")
#seg_idx_filt_pos$cnv_pos = paste(seg_idx_filt_pos$cnv_pos, seg_idx_filt_neg$seg.mean, sep="=")
#par(oma=c(3,3,3,3))
#par(mai=c(10,10,10,10) + 0.2)
#split on colon, then - to get chr # and start interval

pdf(outfile,width=12, height=12)
#plot sample using DNAcopy plot function
patient_plot = plotSample(segment.smoothed.CNA.object,altcol = T, col=c("green", "black"), xaxt = 'n', xlab='Chrom')
#subset indexes to plot chromosomes on X
x_labeling = subset(seg_idx, !duplicated(seg_idx$chrom))
patient_plot = patient_plot + axis(1, at = x_labeling$startRow, labels =x_labeling$chrom)

#if there are positive values, label them
if (nrow(seg_idx_filt_pos[1]) != 0){
  #seg_idx_filt_pos$cnv_pos = paste(seg_idx_filt_pos$cnv_pos, seg_idx_filt_pos$seg.mean, sep="=")
  patient_plot = patient_plot + spread.labels(seg_idx_filt_pos$endRow, seg_idx_filt_pos$seg.mean, labels = seg_idx_filt_pos$cnv_pos, between = T, 
                                              srt =90, pos = 4, cex = 0.7, col="blue")  
}
#if there are negtaive values, label them
if (nrow(seg_idx_filt_neg[1]) !=0){
  #seg_idx_filt_pos$cnv_pos = paste(seg_idx_filt_pos$cnv_pos, seg_idx_filt_neg$seg.mean, sep="=")
  patient_plot = patient_plot + spread.labels(seg_idx_filt_neg$endRow, seg_idx_filt_neg$seg.mean, labels = seg_idx_filt_neg$cnv_pos, between = T, 
                                              srt =90, pos = 2, cex = 0.7 , col="purple")
}
dev.off()
