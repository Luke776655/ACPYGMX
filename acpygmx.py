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
		aa_resi = open(path_to_ff + "//aminoacids.rtp", "r").read().split("\n")
	except:
		print("\nERROR: No such file in directory\n")
		print_help()
	for line in aa_resi:
		if(len(line)>4 and line[0] == '['):
			line = line.split()
			if(len(line[1])>=1 and len(line[1])<=3):
				resi_database.append(line[1])
	resi_database.append('HIS')

	try:
		dna_resi = open(path_to_ff + "//dna.rtp", "r").read().split("\n")
	except:
		print("\nERROR: No such file in directory\n")
		print_help()
	for line in dna_resi:
		if(len(line)>4 and line[0] == '['):
			line = line.split()
			if(len(line[1])>=1 and len(line[1])<=3):
				resi_database.append(line[1])

	try:
		rna_resi = open(path_to_ff + "//rna.rtp", "r").read().split("\n")
	except:
		print("\nERROR: No such file in directory\n")
		print_help()
	for line in dna_resi:
		if(len(line)>4 and line[0] == '['):
			line = line.split()
			if(len(line[1])>=1 and len(line[1])<=3):
				resi_database.append(line[1])
	return resi_database

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
	charge = 0
	#funkcja przypisująca wartości zmiennym
	def get_it(self, line):
		if(len(line)>=70):
			self.record_name = line[0:6].strip()
		if(self.record_name == 'ATOM' or self.record_name == 'HETATM'):
			self.is_record_atomic = True
			self.atom_id = int(line[6:11])
			self.atom_name = line[12:16].strip()
			self.resi_name = line[17:20].strip()
			self.chain = line[21]
			self.resi_id = int(line[22:26])
			if(line[78].isdigit() == True):
				if(line[79] == '+'):
					self.charge = int(line[78])
				elif(line[79] == '-'):
					self.charge = -1*int(line[78])
			

#klasa ze strukturą lokalizacji wybranej reszty chemicznej
class ResiLocation:
	def __init__(self):
		self.ResiLocation = []
	chain = ''
	resi_id = 0
	record_lines = []

#funkcja wyszukująca niestandardowe reszty chemiczne w pliku PDB, zwraca je i ich pozycje jako słownik
def get_nonstandard_resis(pdb, database):
	resis = {}
	for i in range(len(pdb)):
		line = PdbLine()
		line.get_it(pdb[i])
		if(line.is_record_atomic == True):
			if((line.resi_name in database) == False and (line.resi_name in resis) == False):
				resis[line.resi_name] = ResiLocation()
				resis[line.resi_name].chain = line.chain
				resis[line.resi_name].resi_id = line.resi_id
				resis[line.resi_name].record_lines = []
				resis[line.resi_name].record_lines.append(i)
			elif((line.resi_name in resis) == True):
				if(resis[line.resi_name].chain == line.chain and resis[line.resi_name].resi_id == line.resi_id):
					resis[line.resi_name].record_lines.append(i)
	return resis

#funkcja wczytująca wybrane reszty chemiczne z pliku PDB i zapisująca je do osobnych plików PDB
def split_pdb_by_resi(pdb, resis):
	for i in resis:
		subprocess.run(["mkdir", i])
		try:
			resi_pdb = open('./'+i+'/'+i+'_no_H.pdb', "w+")
		except:
			print("\nERROR: Cannot create a file\n")
			print_help()
		resi_pdb.write("HEADER " + i + '\n')
		for j in range(len(resis[i].record_lines)):
			resi_pdb.write(source_pdb[resis[i].record_lines[j]] + '\n')
		resi_pdb.write("END")
		resi_pdb.close()

#funkcja odczytująca ładunek związku z pliku PDB
def get_charge_from_pdb(pdb_path):
	charge = 0
	try:
		molecule = open(pdb_path, "r").read().split("\n")
	except:
		print("\nERROR: No such file in directory\n")
		print_help()
	for i in molecule:
		line = PdbLine()
		line.charge = 0
		line.get_it(i)
		charge += line.charge
	return charge

#funkcja tworząca topologię dla wybranych reszt chemicznych
def make_resi_topology(resis):
	for i in resis:
		#otwieranie skryptu dodającego wodory
		subprocess.run(["babel", "-ipdb", i+"_no_H.pdb", "-opdb", i+".pdb", "-p"], cwd='./'+i)
		#odczytanie ładunku uwodornionej reszty chemicznej
		charge = get_charge_from_pdb('./'+i+'/'+i+'.pdb')
		#otwieranie skryptu generującego topologię
		subprocess.run(["acpype", "-i", i+".pdb", "-n", str(charge)], cwd='./'+i)


#START PROGRAMU

#Sprawdzanie poprawności argumentów wejściowych programu
if(len(sys.argv)>5):
	print("\nERROR: To much input values\n")
	print_help()

source_pdb_name = ''
topology_name = 'topol.top'
nonstandard_resis = {}

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

#pobieranie bazy danych reszt chemicznych z GROMACS
amberff_path = "./amber99sb-ildn.ff"
gmx_database = get_gmx_resi_database(amberff_path)

#wyszukiwanie niestandardowych reszt chemicznych
nonstandard_resis = get_nonstandard_resis(source_pdb, gmx_database)

#wyświetlanie niestandardowych reszt cheemicznych
for i in nonstandard_resis:
	print(i, nonstandard_resis[i].chain, nonstandard_resis[i].resi_id, nonstandard_resis[i].record_lines, '\n')

#wczytywanie niestandardowych reszt chemicznych z pliku PDB i zapisywanie ich do osobnych plików PDB
split_pdb_by_resi(source_pdb, nonstandard_resis)

#tworzenie topologii dla niestandardowych reszt chemicznych
make_resi_topology(nonstandard_resis)
