##V1.2.4
##Run this script as : python flagAdjacentSNPs.py <Sample.canon_annotated> <Sample.canon_annotated_flagged>
##This script goes through a sample's cannon_annotated file and flags the adjacent snp positions (occuring in same Codon).



import sys
import re
import string
from collections import defaultdict


Input=open(sys.argv[1],'r')
Output=open(sys.argv[2],'w')

dict1={}


while 1:
	InputLine=Input.readline()
	if(InputLine):
		InputLine=InputLine.rstrip()
		InputCompo=InputLine.split("\t")
		if(re.match("exonic",InputCompo[25]) and re.match("SNP",InputCompo[5])):
			dict1.update({InputCompo[1]+"\t"+InputCompo[2]:InputCompo[19]+"\t"+InputCompo[13]})
	else:
		break


##print dict1.items()
dict2=defaultdict(list)
for key, value in dict1.iteritems():
	dict2[value].append(key)

for pair in dict2.items():
	if(len(pair[1])==1):
		del dict2[pair[0]]
	else:
		continue

Input.seek(0)    # resets the file pointer
##print dict2.items()
while 1:
        InputLine=Input.readline()
        if(InputLine):
                InputLine=InputLine.rstrip()
                InputCompo=InputLine.split("\t")
		CodonNumberWithTranscript=InputCompo[19]+"\t"+InputCompo[13]
		for key_codon in dict2.keys():
			if(key_codon==CodonNumberWithTranscript):
				InputCompo[len(InputCompo)-1]='1'
				break
			else:
				continue
		Output.write("\t".join(InputCompo[0:]))
		Output.write("\n")
		
	else:
		break

Output.close()
Input.close()

				
	
				
			





##for pair in dict1.itmes():
##	if pair[1] not in dict2.keys():
##		dict2[pair[1]]=[]
##	dict2[pair[1]].append(pair[0])




	
	



	
	
	
			
			
			
			
			
	
	
	
	
	
