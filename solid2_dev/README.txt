Basic Structure
1) Solid2.py is an object class pertaining to the CPD HaloPlex Solid2 Panel
        Class variables = Adapter1 & 2 sequences, bed file location, fragmentsize, min indel count & fraction, library name.
        Instance variables = Sample name, read 1 path, read 2 path, index path, molecular barcode, output directory
2) Haloplex_Run.py - takes instance variables of Solid2 as arguments and passes them to "sample run" method.
        All method calls are called from CPD_ETR.py

3) Helper files = 
        /methods/CPD_ETR.py library of all used methods
        /util/Paths.py - File with paths for tools/databases
        /bed_files/ folder of all used bed files.
        /panel/ holds all panels, i.e Solid2.py

Usage: 
1) Two scripts are required to call this program, this first is the script to run HaloPlex_RUN.py
   The second is a job submission script for HPC job submission

example: BSUB.SH


bsub -M 24567 -e /project/cpdlab/Solid2/HiSeqSamples/CPDV150986-35/S2_35.e -o /project/cpdlab/Solid2/HiSeqSamples/CPDV150986-35/S2_35.o sh CPDV150986-35.sh

example sample.sh
python -m cProfile -o /project/cpdlab/Solid2/HiSeqSamples/CPDV150986-35/CPDV150986-35.timeStats.profile \
HaloPlex_Run.py CPDV150986-35 /project/cpdlab/Solid2/HiSeqSamples/CPDV150986-35/CPDV150986-35.R1.fastq.gz \
/project/cpdlab/Solid2/HiSeqSamples/CPDV150986-35/CPDV150986-35.R2.fastq.gz CCTAATCC \
/project/cpdlab/Solid2/HiSeqSamples/CPDV150986-35/CPDV150986-35.I2.fastq.gz /project/cpdlab/Solid2/HiSeqSamples/CPDV150986-35/


if you have immediate issues calls 610-331-8883 and email ashbigdeli@gmail.com