import subprocess
import sys

InputFileName=sys.argv[1]
Input=open(InputFileName,'r')
Output=open(sys.argv[2],'w')

n_lines=subprocess.Popen(['wc','-l',InputFileName],stdout=subprocess.PIPE).communicate()[0]
n_lines=int(n_lines.split()[0])

i=0

while(i<n_lines):
	i=i+1
	line=Input.readline()
	if(i%4==0):
		Output.write(line)
	else:
		line=line.rstrip()
		Output.write(line+"\t")



Output.close()
Input.close()

