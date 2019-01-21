#!/usr/bin/python3.5

import sys
import subprocess

'''
	ACPYGMX

	Requirements:
		Python 3.5 or higher
		Gromacs 5.0 or higher
		OpenBabel
		AmberTools
		Acpype
                 
	It was inspired by Acpype

	This code is released under GNU General Public License V3.

	<<<  NO WARRANTY AT ALL!!!  >>>

	Łukasz Radziński

    	s174023@student.pg.edu.pl
'''

def print_help():
	'''Wświetlanie pomocy'''
	print("ACPYGMX v. 0.0\n")
	print("Program usage:\n")
	print("acpygmx.py -f _inputname.pdb_ [-o _topologyname.top_]")
	print("or")
	print("acpygmx.py -h\n")
	print("Options:\n\t-f _inputname.pdb_\tname of input PDB file\n\t-o _topologyname.top_\tname of output topology file\n\t-h\t\t\thelp")
	sys.exit()

def get_gmx_resi_database(path_to_ff):
	'''Wczytywanie bazy danych reszt chemicznych.
	Argumenty: ścieżka do katalogu z polem siłowym.
	Zwraca listę z resztami chemicznymi z bazy danych.'''
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

class PdbLine:
	'''Struktura wiersza PDB'''
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
	connected_atoms = []
	def get_it(self, line):
		'''Wczytywanie wartości z wiersza PDB i przypisywanie ich zmiennym ze struktury.
		Argumenty: wiersz PDB'''
		if(len(line)>=60):
			self.record_name = line[0:6].strip()
		if(self.record_name == 'ATOM' or self.record_name == 'HETATM'):
			self.is_record_atomic = True
			self.atom_id = int(line[6:11])
			self.atom_name = line[12:16].strip()
			self.resi_name = line[17:20].strip()
			self.chain = line[21]
			self.resi_id = int(line[22:26])
			if(len(line)>78 and line[78].isdigit() == True):
				if(line[79] == '+'):
					self.charge = int(line[78])
				elif(line[79] == '-'):
					self.charge = -1*int(line[78])
		if(self.record_name == 'CONECT'):
			line = line[6:80].split()
			self.connected_atoms = []
			for i in line:
				self.connected_atoms.append(int(i))
		
class ResiLocation:
	'''Struktura lokalizacji danej reszty chemicznej w PDB'''
	def __init__(self):
		self.ResiLocation = {}
	location = {}
	bounded = False

def get_nonstandard_resis(pdb, database):
	'''Wyszukiwanie niestandardowych reszt chemicznych w pliku PDB.
	Argumenty: tekst pliku PDB podzielony na wiersze, lista reszt chemicznych z bazy danych.
	Zwraca słownik z niestandardowymi resztami chemicznymi i ich pozycjami'''
	resis = {}
	termination_line = 0
	for i in range(len(pdb)):
		line = PdbLine()
		line.get_it(pdb[i])
		if(line.record_name == 'TER'):
			termination_line = i
	print("Termination", termination_line)
	for i in range(len(pdb)):
		line = PdbLine()
		line.get_it(pdb[i])
		if(line.is_record_atomic == True):
			if((line.resi_name in database) == False):
				if((line.resi_name in resis) == False):
					resis[line.resi_name] = ResiLocation()
					resis[line.resi_name].location = {}
				if((line.chain in resis[line.resi_name].location) == False):
					resis[line.resi_name].location[line.chain] = {}
				if((line.resi_id in resis[line.resi_name].location[line.chain]) == False):
					resis[line.resi_name].location[line.chain][line.resi_id] = []
				resis[line.resi_name].location[line.chain][line.resi_id].append(i)
				if(i < termination_line):
					resis[line.resi_name].bounded = True
	return resis

index = 0
def split_pdb_by_resi(pdb, resis):
	'''Wczytywanie wybranych reszt chemicznych z pliku PDB, zapisywanie ich do osobnych plików PDB, dodanie atomów wodoru do każdej.
	Argumenty: tekst pliku PDB podzielony na wiersze, słownik z resztami chemicznymi i ich pozycjami.
	Zwraca zmodyfikowany wejściowy plik PDB, do niestandardowych reszt chemicznych zostały dodane atomy wodoru.'''
	for i in resis:
		subprocess.run(["mkdir", i])

		for j in resis[i].location:
			subprocess.run(["mkdir", j], cwd='./'+i)

			for k in resis[i].location[j]:
				try:
					resi_pdb = open('./'+i+'/'+j+'/'+str(k)+'_no_H.pdb', "w+")
				except:
					print("\nERROR: Cannot create a file\n")
					print_help()
				resi_pdb.write("HEADER "+i+'/'+j+'/'+str(k)+'\n')
				for l in range(len(resis[i].location[j][k])):
					m = 0
					
					while(m<len(pdb)):
						line = PdbLine()
						line.get_it(pdb[m])
						if(i == line.resi_name and j == line.chain and k == line.resi_id):
							#print("INDEX:", index)
							global index
							index = m
							#print(m)
							resi_pdb.write(pdb[m] + '\n')
							del(pdb[m])
							m-=1
						m+=1
				resi_pdb.write("END")
				resi_pdb.close()
				add_hydrogens(str(k), "./"+i+"/"+j+"/", resis[i].bounded)
				hydrogened = open("./"+i+"/"+j+"/"+str(k)+".pdb", "r").read().split("\n")
				h_counter = 0
				for n in hydrogened:
					line = PdbLine()
					line.get_it(n)
					if(line.is_record_atomic == True):
						if(line.atom_name == 'H'):
							if(h_counter > 0):
								h_name = 'H'+str(h_counter)
								h_name = h_name + '   '
								n = n[0:13]+h_name[0:4]+n[17:]
							h_counter +=1
						pdb.insert(index, n)
						index+=1
	return pdb
						
def add_hydrogens(name, path, is_bounded):
	'''Dodawanie wodoru do reszty chemicznej. Korzystanie z programu zewnętrznego "babel"
	Argumenty: nazwa reszty chemicznej, ścieżka do pliku PDB z resztą chemiczną, informacja o tym czy reszta jest w stanie związanym.'''
	subprocess.run(["babel", "-ipdb", name+"_no_H.pdb", "-opdb", name+".pdb", "-p"], cwd=path)
	if(is_bounded == True):
		c_peptide_id = 0
		hc_peptide_id = 0
		n_peptide_id = 0
		hn_peptide_ids = []
		h_to_remove_ids = []
		file_h = open(path+"/"+name+".pdb", "r+")
		pdb = file_h.read().split("\n")
		for j in range(len(pdb)):
			line = PdbLine()
			line.get_it(pdb[j])
			if(line.atom_name == 'C'):
				c_peptide_id = line.atom_id
				pdb[j] = pdb[j][0:78] + "  "
			if(line.atom_name == 'N'):
				n_peptide_id = line.atom_id
				pdb[j] = pdb[j][0:78] + "  "
		for j in range(len(pdb)):
			line = PdbLine()
			line.get_it(pdb[j])
			if(line.record_name == 'CONECT' and len(line.connected_atoms) == 2):
				if((c_peptide_id in line.connected_atoms) == True):
					for k in line.connected_atoms:
						if(k != c_peptide_id):
							hc_peptide_id = k
				if((n_peptide_id in line.connected_atoms) == True):
					for k in line.connected_atoms:
						if(k != n_peptide_id):
							hn_peptide_ids.append(k)
		hn_peptide_ids.pop()
		h_to_remove_ids = hn_peptide_ids
		h_to_remove_ids.append(hc_peptide_id)
		file_h.truncate(0)
		file_h.seek(0)
		for j in range(len(pdb)):
			line = PdbLine()
			line.get_it(pdb[j])
			if(((line.atom_id in h_to_remove_ids) == False or line.atom_name[0] != 'H') and line.record_name != 'CONECT'):
				file_h.write(pdb[j] + '\n')
		file_h.close()

def get_charge_from_pdb(pdb_path):
	'''Odczytywanie ładunku związku w pliku PDB.
	Argumenty: ścieżka do pliku PDB.
	Zwraca wartość ładunku.'''
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

def make_resi_topology(resis):
	'''Tworzenie topologii cząstkowych reszt chemicznych w postaci plików *.itp. Korzystanie z programu zewnętrznego "acpype".
	Argumenty: słownik z resztami chemicznymi i ich pozycjami.'''
	for i in resis:
		path = i+'/'+list(resis[i].location.keys())[0]+'/'
		source_file = str(list(resis[i].location[list(resis[i].location.keys())[0]].keys())[0])+'.pdb'
		subprocess.run(["cp", source_file, i+'.pdb'], cwd=path)
		#odczytanie ładunku uwodornionej reszty chemicznej
		charge = get_charge_from_pdb(path+i+'.pdb')
		#otwieranie skryptu generującego topologię
		subprocess.run(["acpype", "-i", i+'.pdb', "-n", str(charge)], cwd=path)
		subprocess.run(["cp", "./"+i+".acpype/"+i+"_GMX.itp", "../"], cwd=path)

class ItpLine:
	'''Struktura wiersza pliku ITP'''
	def __init__(self):
		self.ItpLine = []
	atom_nr = 0
	atom_type = ''
	resi = 0
	res = ''
	atom_name = ''
	cgnr = 0
	charge = 0
	mass = 0
	def get_atom_line(self, line):
		'''Wczytywanie wartości z wiersza pliku *.itp i przypisywanie ich zmiennym ze struktury.
		Argumenty: wiersz pliku *.itp'''
		line = line.split()
		self.atom_nr = int(line[0])
		self.atom_type = line[1]
		self.resi = int(line[2])
		self.res = line[3]
		self.atom_name = line[4]
		self.cgnr = int(line[5])
		self.charge = float(line[6])
		self.mass = float(line[7])
	'''
	def get_bound_line(self, line):
		line = line.split()
		self.ai = int(line[0])
		self.aj = int(line[1])
		self.funct = int(line[2])
		self.r = float(line[6])
		self.mass = float(line[7])'''

def make_rtp(resis):
	'''Konwertowanie pliku *.itp do *.rtp. Plik *.rtp jest zapisywany w katalogu pola siłowego.
	Argumenty: słownik z resztami chemicznymi i ich pozycjami.'''
	for i in resis:
		print(i)
		try:
			itp = open('./'+i+'/'+i+'_GMX.itp', "r").read().split("\n")
		except:
			print("\nERROR: Cannot create a file\n")
			print_help()
		
		atoms_record_begin = 0
		bonds_record_begin = 0
		bonds_record_end = 0

		for k in range(len(itp)):
			line = itp[k].split()
			if(len(line) >= 3 and line[0] == "[" and line[1] == "atoms"):
				atoms_record_begin = k
			if(len(line) >= 3 and line[0] == "[" and line[1] == "bonds"):
				bonds_record_begin = k
			if(len(line) >= 3 and line[0] == "[" and line[1] == "pairs"):
				bonds_record_end = k
		try:
			rtp = open(amberff_path+'/'+i+'.rtp', "w+")
		except:
			print("\nERROR: Cannot create a file\n")
			print_help()
		rtp.write("[ bondedtypes ]" + '\n')
		rtp.write("; bonds  angles  dihedrals  impropers all_dihedrals nrexcl HH14 RemoveDih" + '\n')
		rtp.write("     1       1          9          4        1         3      1     0" + '\n\n')
		rtp.write("[ "+i+" ]" + '\n')
		rtp.write(" [ atoms ]" + '\n')
		j = 1
		for k in range(atoms_record_begin, bonds_record_begin):
			line = itp[k].split()
			if(len(line) > 10 and line[0] != ';'):
				line = ItpLine()
				line.get_atom_line(itp[k])
				rtp.write("\t"+line.atom_name+"\t"+line.atom_type+"\t"+'{0:.6f}'.format(line.charge)+"\t"+str(j)+"\n")
				j +=1
		rtp.close()
			


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
	print(i, nonstandard_resis[i].location, '\n')

#wczytywanie niestandardowych reszt chemicznych z pliku PDB, zapisywanie ich do osobnych plików PDB, uzyskanie układu z uwodornionymi resztami niestandardowymi
complex_pdb = split_pdb_by_resi(source_pdb, nonstandard_resis)

#zapisywanie complex_pdb do pliku complex.pdb
file_h = open("./complex.pdb", "w+")
for i in complex_pdb:
	file_h.write(i + '\n')
file_h.close()


#tworzenie topologii dla niestandardowych reszt chemicznych
make_resi_topology(nonstandard_resis)

#dodawanie plików *.rtp z topologiami związków niestandardowych do katalogu pola siłowego
make_rtp(nonstandard_resis)

#tworzenie topologii ogólnej, używanie programu zewnętrznego "gmx pdb2gmx"
subprocess.run(["gmx", "pdb2gmx", "-f", "complex.pdb", "-o", "complex.gro", "-p", topology_name])

print("\nTopology is written to " + topology_name)
