import glob,sys,subprocess
from util import Paths

def parse_sheet (sample_sheet):    
    run_args = []
    with open(sample_sheet, 'r') as f:
        lines = f.readlines()[1:] #skips headers
        for line in lines:
            sample_parse = line.split(',')
            run_args.append(sample_parse[2] +',' + sample_parse[4] + ',' + sample_parse[9])
    return run_args
    
def demultiplex(run_dir, sample_sheet):
    subprocess.Call(Paths.bcl2 + '--input-dir '+ run_dir + '/Data/Intensities/BaseCalls --output-dir ' + run_dir + 
                    ' Unaligned --sample-sheet ' + sample_sheet + ' --no-eamss --use-bases-mask Y150n,I8,Y10,Y150n --mismatches 1', shell = True)

# concat fasta
#
# @param1 = filenames: all files to be concatenated.
# @param2 = download_to: folder to concatenate files.
#
# Opens ALL files in folder and writes them to one .fasta file.
def concat_file( db_name, download_to):
    # create file list of all files to concatenate
    files_in_path = download_to 
    filenames = glob.glob(files_in_path)
    db_path = download_to + db_name + '.vcf'
    try:
        #open each and write to single file
        with open( db_path, 'w') as outfile:
            for filename in filenames:
                with open(filename) as infile:
                    for line in infile:
                        outfile.write(line)
        outfile.close()
    except IOError:
        print "There was an IO error concatinating your files, check your file path and authority!"
        
        
def main():
    run_dir = sys.argv[1]
    sample_sheet = run_dir + '/Data/Intensities/BaseCalls/SampleSheet.csv'
    sample_info = parse_sheet(sample_sheet)
    demultiplex(run_dir, sample_sheet)
    
main()