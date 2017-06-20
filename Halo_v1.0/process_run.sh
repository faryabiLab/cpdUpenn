#!/bin/sh

run_dir=$1
out_dir=$2
demux_sh='/data1/HiSeqRun/demultiplex.sh'

bsubs=$(ssh cpdlab@170.212.141.107 sh ${demux_sh} ${run_dir} ${out_dir})
echo ${bsubs}
./${bsubs}
