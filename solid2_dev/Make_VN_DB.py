# MakeVNDB.py
# Ashkan Bigdeli
# 
# Downloads 1000G data and concats it
# Usage: python MakeRedDB.py /path/to/place/sequences
#


import os, gzip, glob, os.path, sys
from ftplib import FTP


# download
#
# @param1 = filenames: list of all files in ftp available for download.
# @param2 = file_id: specify which files in folder to be downloaded.
# @param3 = download_to: folder to place downloaded files.
# @param4 = ftp: ftp site from which to download
#
# Will download all specified files to a designated folder.
def download(filenames, file_id, download_to,ftp):    
    for filename in filenames:
        if filename.endswith(file_id):
            try:
                local_filename = os.path.join(download_to, filename)
                file = open(local_filename, 'wb')
                ftp.retrbinary('RETR ' + filename, file.write)
                file.close()
            except IOError:
                print "There has been an IO error! Please check all file paths."

# unzip
#
# @param1 = download_to: folder to place downloaded files.
#
# Uses gzip import to unzip files within the same folder.
def unzip(download_to):
    
    try:
        # for each file in the directory
        for gzip_path in glob.glob(download_to + "/*"):
            if os.path.isdir(gzip_path) == False:
                in_file = gzip.open(gzip_path, 'rb')
            
                # uncompress the file into temp
                temp = in_file.read()
                in_file.close()

                # get the filename
                gzip_filename = os.path.basename(gzip_path)
            
                # get original filename and remove the extension
                filename = gzip_filename[:-3]
                uncompressed_path = os.path.join(download_to, filename)

                # write uncompressed file
                open(uncompressed_path, 'w').write(temp)
    except Exception:
        print "There was a problem unzipping your files :( Check GZip."

# remove
#
# @param1 = suffix: file type to be removed.
# @param2 = download_to: folder to remove downloaded files.

# Removes all given file types from folder.
def remove(suffix, download_to):
    try:
        files_in_path = download_to + suffix
        files = glob.glob(files_in_path)
        for f in files:
            os.remove(f)
    except IOError:
        print "There was an IO error deleting your files, check your file path and authority!"
                         
# concat fasta
#
# @param1 = filenames: all files to be concatenated.
# @param2 = download_to: folder to concatenate files.
#
# Opens ALL files in folder and writes them to one .fasta file.
def concat_fasta(db_name, download_to):
    # create file list of all files to concatenate
    files_in_path = download_to + '*'
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
    
    # get system specific information from user
    if len(sys.argv) !=2:
        print "Example: python  Make_VN_DB.py /project/cpdlab/ashkan/vn_test/"
    
    download_to =sys.argv[1]
    

    # file type to download
    file_id = '.genotypes.vcf.gz'
    
    # set DB name
    db_name = "1000G_All"

    # type of blast database, i.e 'nucl' = nucleotide
    
    # clear folder of previous entries or databases
    remove('*', download_to)
    
    try:
        # set ftp location
        ftp = FTP('ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/')
        
        # login to ftp
        # Omitting ftp.login('USERNAME', 'PASSWORD') will login as anonymous
        ftp.login()
    
        # Change directory in ftp to navigate to desired genome
        #ftp.cwd('/vol1/ftp/release/20130502/')
    
        # Create a list of the filenames in this location
        filenames = ftp.nlst()
        print filenames
        
        print "Downloading the latest Reference Genome ...."
        download(filenames, file_id, download_to,ftp)
        ftp.quit()
    except Exception:
        print "There has been a download issue, please check status of the ftp."
    
    print "Unzipping downloads.."
    unzip(download_to)
    
    # remove all compressed files
    #remove('*.gz', download_to)
    
    # concate files and remove once concatenated
    print "Concatinating files..."
    concat_fasta(db_name, download_to)
    #remove('*.fna', download_to)
    
    # profit!
    print "Complete! Please check " + download_to + " to begin use :)"
main()                                 