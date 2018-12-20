#!/usr/bin/python3.5

#importowanie modułów
import sys
import subprocess

#funkcja wyświetlająca pomoc
def print_help():
	print("ACPYGMX v. 0.0\n")
	print("Program usage:\n")
	print("acpygmx.py -f _inputname.pdb_ [-o _topologyname.top_]")
	print("or")
	print("acpygmx.py -h\n")
	print("Options:\n\t-f _inputname.pdb_\tname of input PDB file\n\t-o _topologyname.top_\tname of output topology file\n\t-h\t\t\thelp")
	sys.exit()

#funkcja pobierająca bazę danych reszt chemicznych z GROMACS
def get_gmx_resi_database(path_to_ff):
	resi_database = []
	try:
		aa_resi = open(path_to_ff + "//aminoacids.hdb", "r").read().split("\n")
	except:
		print("\nERROR: No such file in directory\n")
		print_help()
	for line in aa_resi:
		line = line.split("\t")
		if(len(line[0])>2 and len(line[0])<6 and line[0][0].isupper() == True):
			resi_database.append(line[0])
	resi_database.append('HIS')

	try:
		dna_resi = open(path_to_ff + "//dna.hdb", "r").read().split("\n")
	except:
		print("\nERROR: No such file in directory\n")
		print_help()
	for line in dna_resi:
		line = line.split("\t")
		if(len(line[0])>2 and len(line[0])<6 and line[0][0].isupper() == True):
			resi_database.append(line[0])

	try:
		rna_resi = open(path_to_ff + "//rna.hdb", "r").read().split("\n")
	except:
		print("\nERROR: No such file in directory\n")
		print_help()
	for line in rna_resi:
		line = line.split("\t")
		if(len(line[0])>2 and len(line[0])<6 and line[0][0].isupper() == True):
			resi_database.append(line[0])
	return resi_database

#funkcja zapisująca wybrane reszty chemiczne do plików PDB
def write_resi_to_pdb(resis, source_pdb):
	for i in resis:
		subprocess.run(["mkdir", i])
		try:
			resi_pdb = open('./'+i+'/'+i+'.pdb', "w+")
		except:
			print("\nERROR: Cannot create a file\n")
			print_help()
		resi_pdb.write("HEADER " + i + '\n')
		for j in range(len(resis[i].record_line)):
			resi_pdb.write(source_pdb[resis[i].record_line[j]] + '\n')
		resi_pdb.write("END")
		resi_pdb.close()

#funkcja tworząca topologię dla wybranych reszt chemicznych
def make_resi_topology(resis):
	for i in resis:
		#otwieranie skryptu dodającego wodory
		subprocess.run(["babel", "-ipdb", i+".pdb", "-opdb", i+"_h.pdb", "-h"], cwd='./'+i)
		#otwieranie skryptu generującego topologię
		subprocess.run(["acpype", "-i", i+"_h.pdb"], cwd='./'+i)

#klasa ze strukturą wiersza PDB
class PdbLine:
	def __init__(self):
		self.PdbLine = []
	is_record_atomic = False
	record_name = ''
	atom_id = 0
	atom_name = ''
	resi_name = ''
	chain = ''
	resi_id = 0
	#funkcja przypisująca wartości zmiennym
	def get_it(self, line):
		if(len(line)>=50):
			self.record_name = line[0:6].strip()
		if(self.record_name == 'ATOM' or self.record_name == 'HETATM'):
			self.is_record_atomic = True
			self.atom_id = int(line[6:11])
			self.atom_name = line[12:16].strip()
			self.resi_name = line[17:20].strip()
			self.chain = line[21]
			self.resi_id = int(line[22:26])

#klasa ze strukturą lokalizacji wybranej reszty chemicznej
class ResiLocation:
	def __init__(self):
		self.ResiLocation = []
	chain = ''
	resi_id = 0
	record_line = []


#START PROGRAMU

#Sprawdzanie poprawności argumentów wejściowych programu
if(len(sys.argv)>5):
	print("\nERROR: To much input values\n")
	print_help()

source_pdb_name = ''
topology_name = 'topol.top'
nonstandard_resi = {}

#obsługa opcji wejściowych
for i in range(len(sys.argv)):
	if(sys.argv[i] == '-h'):
		print_help()
	elif(sys.argv[i] == '-f' and sys.argv[i+1][0] != '-'):
		source_pdb_name = sys.argv[i+1]
	elif(sys.argv[i] == '-o' and sys.argv[i+1][0] != '-'):
		topology_name = sys.argv[i+1]
	elif(sys.argv[i][0] == '-'):
		print("\nERROR:Unknown option: " + sys.argv[i] + "\n")
		print_help()

print("PDB file name: " + source_pdb_name)
print("Topology file name: " + topology_name)

#otwieranie źródłowego pliku PDB
try:
	source_pdb = open(source_pdb_name, "r").read().split("\n")
except:
	print("\nERROR: No such file in directory\n")
	print_help()

ff_path = "./amber99sb-ildn.ff"
GMXdatabase = get_gmx_resi_database(ff_path)

#wyszukiwanie niestandardowych reszt chemicznych
for i in range(len(source_pdb)):
	pdb_line = PdbLine()
	pdb_line.get_it(source_pdb[i])
	if(pdb_line.is_record_atomic == True):
		if((pdb_line.resi_name in GMXdatabase) == False and (pdb_line.resi_name in nonstandard_resi) == False):
			nonstandard_resi[pdb_line.resi_name] = ResiLocation()
			nonstandard_resi[pdb_line.resi_name].chain = pdb_line.chain
			nonstandard_resi[pdb_line.resi_name].resi_id = pdb_line.resi_id
			nonstandard_resi[pdb_line.resi_name].record_line.append(i)
		elif((pdb_line.resi_name in nonstandard_resi) == True):
			if(nonstandard_resi[pdb_line.resi_name].chain == pdb_line.chain and nonstandard_resi[pdb_line.resi_name].resi_id == pdb_line.resi_id):
				nonstandard_resi[pdb_line.resi_name].record_line.append(i)


for k in nonstandard_resi:
	print(k, nonstandard_resi[k], '\n')

write_resi_to_pdb(nonstandard_resi, source_pdb)
make_resi_topology(nonstandard_resi)
