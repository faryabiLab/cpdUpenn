#Tolga Ayazseven
#1/11/2018

#This program is made to read MastervarFinal, create a unique key to look for the same key through the file.
#This program will help to identify Cross Contaminations between Samples.
#The output of this program will be a csv file which has Sample_ID, chromosome_Occurance count_Position, Snp_ID, Alternative, Reference, and the Mutation type

import sys
import os
import csv
import pandas as pd

folder_name=sys.argv[1]
path_1='/project/cpdlab/HiSeqRun/'
path_in=str(path_1+folder_name)
for file in os.listdir(path_in):  # Read the files in the current directory
    if file.endswith("Run_masterVarFinal.txt"):  # find the final variance file
        final_var = str(file)  # get the file name
        sample = os.path.splitext(final_var)[0]  # get the Project name with "Run_masterVarFinal" at the end

f_name=str(path_in+'/'+final_var)
sample_id = sample.replace("Run_masterVarFinal",
                           "Contamination.csv")  # replacing "Run_masterVarFinal" with "Contamination" file name
f_out=str(path_in+'/'+sample_id)

mutation_type = []
sample_name = []
chrom = []
pos = []
ref = []
alt = []
key_l = []
allel_freq = []
first_col = []
snps_id = []
keys = set()

with open(f_name, "r") as f:
    csvreader = csv.reader(f, delimiter='\t')
    next(csvreader)  # Skip the first line
    for line in csvreader:
        key0 = line[1] + line[2] + line[3] + line[4]  # Unique Identifier for the contamination var_type != SNP
        # key1 = line[0] + line[1] +line[2] + line[3] + line[4]
        key_id = line[16]  # Mutation type
        s_name = line[0]  # Sample Name
        chromo = line[1]  # chromosome #
        posit = line[2]  # variant position
        reference = line[3]  # Ref.
        alter = line[4]  # Alternative
        allel = line[54]  # Allel Frequency
        snp_id = line[24]
        if key0 in keys:
            mutation_type.append(key_id)
            sample_name.append(s_name)
            chrom.append(chromo)
            pos.append(posit)
            ref.append(reference)
            alt.append(alter)
            key_l.append(str(key0))
            allel_freq.append(str(allel))
            snps_id.append(snp_id)
            first_col.append(str(line[1] + ' / ' + line[2] + ' / '))
        else:
            keys.add(key0)

data = {'Sample_ID': sample_name,
        'Chr/Pos/Count': first_col,
        'Key': key_l,
        'Allel_Freq': allel_freq,
        'Reference': ref,
        'SNP_ID': snps_id,
        'Alternative': alt,
        'Mutation_Type': mutation_type,
        'Value': 1}

# First data frame with all the values for the count step
df = pd.DataFrame(data, columns=['Chr/Pos/Count','Key','Sample_ID', 'Allel_Freq', 'Reference', 'SNP_ID', 'Alternative',
                                 'Mutation_Type', 'Value'])
# df=df.sort_values(by=['Chr/Pos/Count'])
df['Counts'] = df.groupby('Chr/Pos/Count')['Value'].transform(
    sum)  # sum the value'1' to get the count of each contaminated sample
df['Chr_Pos_Count'] = df[['Chr/Pos/Count', 'Counts']].apply(lambda x: ' '.join(x.dropna().astype(str)), axis=1)
df = df[df.Counts != 1]
# Create a final data table
df1 = df.drop(['Chr/Pos/Count', 'Value'], axis=1).sort_values(by=['Chr_Pos_Count'])
#df4.to_csv(f_out)  # Writing a csv file

master_key=[]
categ=[]
#Path to filemaker data file which has category information
path_master='/project/cpdlab/Scripts/SolidV2/solidv2_v1.2/files/Clean_master.csv'
with open(path_master, 'r') as p:
    csv_master = csv.reader(p, delimiter=',')
    next(csv_master)
    for line in csv_master:
        master_key.append(line[1] + line[2] + line[3] + line[4]) #Creating another key which has the same format as the first key
        categ.append(line[5])
        
data2={'Key':master_key,
      'Category':categ}
df2= pd.DataFrame(data2, columns=['Key', 'Category']) #Second data frame with the filemaker key and the categorical information
df2=df2.drop_duplicates(subset=['Key'], keep='first')#Drop the duplicates from the Clean Master file 
df3=pd.merge(df1,df2, on='Key')#Merge with the Matching key
df4=df3.drop(['Key'], axis=1).sort_values(by=['Chr_Pos_Count']) #Drop the key
df4.to_csv(f_out)#Write to csv
