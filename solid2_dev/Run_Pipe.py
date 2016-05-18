import csv,sys,subprocess
from util import Paths

def parse_sheet (sample_sheet, out_dir):    
    run_args = []
    with open(sample_sheet, 'r') as f:
        lines = f.readlines()[1:] #skips headers
        for line in lines:
            sample_parse = line.split['\t']
            run_args.append(sample_parse[2] +'\t' + sample_parse[4] + '\t' + sample_parse[9])
    return run_args
    
def demultiplex(run_dir, sample_sheet):
    subprocess.Call(Paths.bcl2 + '--input-dir '+ run_dir + '/Data/Intensities/BaseCalls --output-dir ' + run_dir + 
                    ' Unaligned --sample-sheet ' + sample_sheet + ' --no-eamss --use-bases-mask Y150n,I8,Y10,Y150n --mismatches 1', shell = True)
    subprocess.Call ("'for j in `awk -F "," '{print $3}' " + sample_sheet + " | grep -v 'Sample_ID' | sort |uniq`" + 
                     " do zcat $j/$j*R1*gz > $j/$j.R1.fastq zcat $j/$j*R2*gz > $j/$j.R2.fastq zcat $j/$j*R3*gz > $j/$j.R3.fastq gzip $j/*fastq", shell = True)

def main():
    run_dir = sys.argv[1]
    sample_sheet = run_dir + '/Data/Intensities/BaseCalls/SampleSheet.csv'
    print sample_sheet
    sample_info = parse_sheet(sample_sheet)
    print 'demultiplexing...'
    demultiplex(run_dir, sample_sheet)
    for item in sample_info:
        print item
    