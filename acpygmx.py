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

def createGMXdatabase(path):
	residatabase = []
	try:
		aa_resi = open(path + "//aminoacids.hdb", "r").read().split("\n")
	except:
		print("\nERROR: No such file in directory\n")
		helpprint()
	for line in aa_resi:
		line = line.split("\t")
		if(len(line[0])>2 and len(line[0])<6 and line[0][0].isupper() == True):
			residatabase.append(line[0])
	residatabase.append('HIS')
	try:
		dna_resi = open(path + "//dna.hdb", "r").read().split("\n")
	except:
		print("\nERROR: No such file in directory\n")
		helpprint()
	for line in dna_resi:
		line = line.split("\t")
		if(len(line[0])>2 and len(line[0])<6 and line[0][0].isupper() == True):
			residatabase.append(line[0])
	try:
		rna_resi = open(path + "//rna.hdb", "r").read().split("\n")
	except:
		print("\nERROR: No such file in directory\n")
		helpprint()
	for line in rna_resi:
		line = line.split("\t")
		if(len(line[0])>2 and len(line[0])<6 and line[0][0].isupper() == True):
			residatabase.append(line[0])
	return residatabase

if(len(sys.argv)>5):
	print("\nERROR: To much input values\n")
	helpprint()

pdb_name = ''
topol_name = 'topol.top'
residuals = {}
nonstdresi = []


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

GMXPATH = "./amber99sb-ildn.ff"
GMXdatabase = createGMXdatabase(GMXPATH)

i = 0
for line in pdb_file:
	line = line.split()
	if(len(line)>=8  and (line[0] == 'ATOM' or line[0] == 'HETATM')):
		if(((line[3] in GMXdatabase) == False) and ((line[3] in residuals) == False)):
			if(line[4].isupper() == True):
				residuals[line[3]] = [(line[4])]
				residuals[line[3]].append(int(line[5]))
			else:
				residuals[line[3]] = [0]
				residuals[line[3]] = [int(line[4])]
			residuals[line[3]].append(i)
		elif((line[3] in residuals) == True) :
				if(line[4].isupper() == True):
					if(residuals[line[3]][0] == line[4] and residuals[line[3]][1] == int(line[5])):
						residuals[line[3]].append(i)
				else:
					if(residuals[line[3]][0] == 0 and residuals[line[3]][1] == int(line[4])):
						residuals[line[3]].append(i)
	i+=1

for k in residuals:
	print(k, residuals[k], '\n')
