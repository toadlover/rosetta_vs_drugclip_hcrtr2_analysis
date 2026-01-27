#a quick download sanity check to confirm that all placement files are represented in the smiles file for each pocket, and will note if there is a discrepancy from either side

import os,sys

#iterate over the directory and go pocket directory by pocket directory
for r,d,f in os.walk(os.getcwd):
	for dire in d:
		if dire == "pocket1" or dire == "pocket2" or dire == "pocket3":
			#print what pocket directory is being looked at
			print("On pocket " + dire)

			#create lists for the ligands represented in the smiles file and placement files
			smiles_ligands = []
			placements_ligands = []

			#read over the smiles file in this folder
			smiles_file = open(dire + "/" + "HCRTR2_" + dire + "_selected_ligands.txt")
			for line in smiles_file.readlines():
				stripped_line = line.strip()
				if len(stripped_line.split()) == 2:
					#get the ligand in the second position
					ligand = stripped_line.split()[1]
					smiles_ligands.append(ligand)

			#now look through the directory and get ligands
			for r2,d2,f2 in os.walk(dire):
				for file in f2:
					if file.endswith(".cif"):
						#the ligand is after the last underscore, lop off the .cif at the end and the "#-" at the start (can't simply parse by - because of PV drugs)
						ligand = file.split("_")[len(file.split("_")) - 1].split(".cif")[0][2:]
						placements_ligands.append(ligand)

			#now iterate over each list and determine what ligands are absent from each list to determine what mismatch exists
			#smiles first
			for lig in smiles_ligands:
				if lig not in placements_ligands:
					print(lig + " not in placements")

			for lig in placements_ligands:
				if lig not in smiles_ligands:
					print(lig + " not in smiles file")