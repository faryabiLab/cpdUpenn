# Solid2.py
#
#Ashkan Bigdeli 3/30/2016
#
# Contains the Solid2 object

class Solid2(object):

    adapter1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    adapter2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
    amplicon_bed='/project/cpdlab/cpdUpenn/solid2_dev/bed_files/39795-1430234657_Covered.bed'
    frag_size='190,50'
    min_indel_cnt = '9'
    min_indel_frac = '0.10'
    lib_name = 'Solid2'

    def __init__(self, sample_name, read1, read2, read_index, index2, out_dir):
        self.sample_name = sample_name
        self.read1 = read1
        self.read2 = read2
        self.read_index = read_index
        self.index2 = index2
        self.out_dir = out_dir
