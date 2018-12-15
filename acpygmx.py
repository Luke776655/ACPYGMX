#!/usr/bin/python3

import sys
import avogadro

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

def writenonstdpdb(resi, pdb_org):
	for i in resi:
		try:
			resi_pdb = open(i + '.pdb', "w+")
		except:
			print("\nERROR: Cannot create a file\n")
			helpprint()
		resi_pdb.write("HEADER " + i + '\n')
		for j in range(2, len(resi[i]), 1):
			resi_pdb.write(pdb_org[resi[i][j]] + '\n')
		resi_pdb.write("END")
		resi_pdb.close()

if(len(sys.argv)>5):
	print("\nERROR: To much input values\n")
	helpprint()

pdb_name = ''
topol_name = 'topol.top'
nonstdresi = {}


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
		if(((line[3] in GMXdatabase) == False) and ((line[3] in nonstdresi) == False)):
			if(line[4].isupper() == True):
				nonstdresi[line[3]] = [(line[4])]
				nonstdresi[line[3]].append(int(line[5]))
			else:
				nonstdresi[line[3]] = [0]
				nonstdresi[line[3]].append(int(line[4]))
			nonstdresi[line[3]].append(i)
		elif((line[3] in nonstdresi) == True) :
				if(line[4].isupper() == True):
					if(nonstdresi[line[3]][0] == line[4] and nonstdresi[line[3]][1] == int(line[5])):
						nonstdresi[line[3]].append(i)
				else:
					if(nonstdresi[line[3]][0] == 0 and nonstdresi[line[3]][1] == int(line[4])):
						nonstdresi[line[3]].append(i)
	i+=1

for k in nonstdresi:
	print(k, nonstdresi[k], '\n')

writenonstdpdb(nonstdresi, pdb_file)
