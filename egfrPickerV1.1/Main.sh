##V1.1
##Changes in V1.1: 1) Added the functionality to calculate the FDP, PercentageUsable and important ratios using CalculateRatio.py script which takes in the RunSummaryFinal.txt file as an input and gives out RunSummaryFinalCalculations.txt as output.
##		   2) Added the functionality to search for the EGFRvIII sequence in each samples's EGFR_exon1_8 btrimmed fastq sequences to confirm the presence of glysine residue at the junction of exon1 and exon8. This script also searches for EGFRwt sequence in EGFR_exon1_exon2 btrimmed fastq sequences.
##                 3) Removed 'Primers/Primer.PrimParts' and 'Primers/Primer_RC.PrimParts' arguments from the command line of 'egfrPick.sh' script.
##                 4) Changed 'RunFinalSummary' to 'RunSummaryFinal' on the 'echo' statement written prior to the command-to-generate the RunSummaryFinal.txt file.
##	           5) Changed 'EGFRv3' to 'EGFRvIII' in the grep statement used while processing SampleSheet.csv to generate Samplesheet_Tobeprocessed.txt. Also added an extra processing step,i.e to replace "_" with "-" in sample-names while processing SampleSheet.csv.
##		   6) 'cp *.py $i' statement is changed to 'cp MakingOfSummary.py $i' to not copy the other python scripts, i.e. CalculateRatio.py, seqToFasta.py and FastqToTableConverter.py into the individual sample folders. 
#!/bin/sh

# Select necessary columns from SampleSheet.csv
echo "Processing SampleSheet.csv to extract the SampleName,Barcodes and Description (FFPE_V1 or HEME_V1) ...\n"
sed -n '/Data/,//p' SampleSheet.csv | sed 1,2d | awk -F "," '{print $2,$6,$8,$10}'|grep EGFRvIII|sed 's/_/-/g' > Samplesheet_Tobeprocessed.txt


#Create a Sub-Directory for each Sample, copy the fastq files into respective sub-directories.

echo "Creating Sample Subdirectories ...\n"

for i in `awk '{print $1}' Samplesheet_Tobeprocessed.txt`
do
mkdir $i
for j in `ls|grep $i |grep fastq`
do
mv $j $i
done
cp egfrPick.sh $i
cp Count.sh $i
cp MakingOfSummary.py $i
cp -r Primers $i
done


#Run Pipe

echo "Running Pipe for each sample.."

for i in `awk '{print $1}' Samplesheet_Tobeprocessed.txt`
do
cd $i
nohup ./egfrPick.sh $i Primers/Manifest.txt Primers/Primer.txt Primers/Primer_RC.txt > $i.egfrPick.out
if [ -s $i.summary ]
then
cat $i.summary >> ../RunSummary.txt
else
echo "Summary File not found for Sample $i..."
fi
cd ..
done


echo "Preparing RunSummaryFinal File.."
sort -k 1,1 -r RunSummary.txt|uniq > RunSummaryFinal.txt

echo "Preparing RunSummaryFinalCalculations File.."
python CalculateRatio.py RunSummaryFinal.txt RunSummaryFinalCalculations.txt

echo "Searching for EGFRvIII sequence in exon1_exon8 btrimmed fastq sequences to confirm the presence of glycine residue. Searching is also performed for EGFRwt sequence in exon1_exon2 btrimmed fastq sequences..."
nohup ./seqSearch.sh > seqSearch.out
























