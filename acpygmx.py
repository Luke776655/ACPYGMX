#!/usr/bin/python3

import sys

def helpprint():
	print("ACPYGMX v. 0.0\n")
	print("Program usage:\n")
	print("acpygmx.py -f _inputname.pdb_ [-o _topologyname.top_]")
	print("or")
	print("acpygmx.py -h\n")
	print("Options:\n\t-f _inputname.pdb_\tname of input PDB file\n\t-o _topologyname.top_\tname of output topology file\n\t-h\t\t\thelp")
	sys.exit()

if(len(sys.argv)>5):
	print("\nERROR: To much input values\n")
	helpprint()

pdb_name = ''
topol_name = 'topol.top'
residuals = {}


for i in range(len(sys.argv)):
	if(sys.argv[i] == '-h'):
		helpprint()
	elif(sys.argv[i] == '-f' and sys.argv[i+1][0] != '-'):
		pdb_name = sys.argv[i+1]
	elif(sys.argv[i] == '-o' and sys.argv[i+1][0] != '-'):
		topol_name = sys.argv[i+1]
	elif(sys.argv[i][0] == '-'):
		print("\nERROR:Unknown option: " + sys.argv[i] + "\n")
		helpprint()

print("PDB file name: " + pdb_name)
print("Topology file name: " + topol_name)

try:
	pdb_file = open(pdb_name, "r").read().split("\n")
except:
	print("\nERROR: No such file in directory\n")
	helpprint()

i = 0
for line in pdb_file:
	line = line.split()
	if(len(line)>=8  and (line[0] == 'ATOM' or line[0] == 'HETATM')):
		if((line[3] in residuals) == False):
			residuals[line[3]] = [i]
		else:
			residuals[line[3]].append(i)
	i+=1

for k in residuals:
	print(k, residuals[k], '\n')
