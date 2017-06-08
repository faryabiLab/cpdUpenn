import sys
InFile=open(sys.argv[1],'r')
OutFile=open(sys.argv[2],'w')
InList=InFile.readlines()

##OutFile.write(header)

for line in InList:
	line=line.rstrip()
	Components=line.split("\t")
	AnnoPart1="\t".join(Components[0:5])   # this means bring columns 0 through 4 over
	AnnoPart2=Components[11]+"\t"+Components[5]
	AnnoPart3="\t".join(Components[26:28])
	AnnoPart4=Components[30]+"\t"+"\t".join(Components[49:53])+"\t"+Components[73]+"\t"+Components[89]+"\t"+Components[24]+"\t"+"\t".join(Components[31:38])+"\t"+"\t".join(Components[17:19])+"\t"+Components[104]+"\t"+Components[135]+"\t"+Components[138]+"\t"+"\t".join(Components[6:11])+"\t"+"\t".join(Components[12:17])+"\t"+"\t".join(Components[19:24])+"\t"+Components[25]+"\t"+"\t".join(Components[28:30])+"\t"+"\t".join(Components[31:41])+"\t"+Components[48]+"\t"+"\t".join(Components[53:73])+"\t"+"\t".join(Components[74:89])+"\t"+"\t".join(Components[90:104])+"\t"+"\t".join(Components[105:135])+"\t"+"\t".join(Components[136:138])+"\t"+"\t".join(Components[139:141])+"\n"
	OutFile.write(AnnoPart1+"\t"+AnnoPart2+"\t"+AnnoPart3+"\t"+AnnoPart4)
OutFile.close()
InFile.close()
