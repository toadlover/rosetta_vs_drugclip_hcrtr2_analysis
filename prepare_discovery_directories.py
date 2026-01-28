#the purpose of this script is to prepare test_params discovery directories of up to 15 conformers per ligand from the pocket ligand smiles file
#this file will use conformator to "rehydrate" the smiles into up to 15 conformers, which we will then use Rosetta's molfile to params script to convert the conformers to params files and organize them for discovery
#this script has hard-coded paths to programs on the umass chan hpc I am using, but paths can be adjusted if this needs to be used elsewhere

#imports
import os,sys

#root location
root_location = os.getcwd()

#iterate over each smiles text file, and for each text file, make a folder for each ligand with conformers, prepared test params, and an arg file for discovery
for r,d,f in os.walk(os.getcwd()):
	for file in f:
		if file.startswith("HCRTR2_pocket") and file.endswith("selected_ligands.txt"):
			#read the file for smiles strings and ligand names per line
			read_file = open(r + "/" + file, "r")
			for line in read_file.readlines():
				line_stripped = line.strip()
				smiles = line_stripped.split()[0]
				ligname = line_stripped.split()[1]

				#make a folder for the ligand in the respective location
				os.chdir(r)
				os.system("mkdir -p " + ligname)
				os.chdir(ligname)

				#write the smiles string to a smi file for conformator to read
				os.system("echo '" + smiles + "' > " + ligname + ".smi")

				#now, use conformator to make up to 15 ligand of the smiles string
				os.system("/pi/summer.thyme-umw/2024_intern_lab_space/conformator_1.2.1/conformator -i " + ligname + ".smi -o confs.mol2 --keep3d --hydrogens -n 15 -v 0")

				#split the conformers file into multiple files
				#example files made: individual_conf_1.mol2, individual_conf_2.mol2, individual_conf_3.mol2, ... individual_conf_250.mol2
				os.system("obabel -i mol2 confs.mol2 -o mol2 -O " + ligname + "_.mol2 -m")

				#convert each conformer to params in a for loop
				for r2,d2,f2 in os.walk(os.getcwd()):
					for file in f2:
						if file.endswith(".mol2") and ligname in file:
							os.system("singularity exec /pi/summer.thyme-umw/2024_intern_lab_space/ari_work/rosetta_discovery_benchmark_test_space_september_2025/rosetta_and_conformator.sif python /rosetta/source/scripts/python/public/molfile_to_params.py " + file + " -n " + file.split(".")[0] + " --keep-names --long-names --clobber --no-pdb")

				#make a test_params folder and move all params files into the test_params
				os.system("mkdir test_params")
				os.system("mv *.params test_params")
				os.chdir("test_params")
				#make necessary files
				os.system("touch exclude_pdb_component_list.txt patches.txt")

				#write the residue_types file
				restype_file = open("residue_types.txt", "w")
				restype_file.write("## the atom_type_set and mm-atom_type_set to be used for the subsequent parameter\n")
				restype_file.write("TYPE_SET_MODE full_atom\n")
				restype_file.write("ATOM_TYPE_SET fa_standard\n")
				restype_file.write("ELEMENT_SET default\n")
				restype_file.write("MM_ATOM_TYPE_SET fa_standard\n")
				restype_file.write("ORBITAL_TYPE_SET fa_standard\n")
				restype_file.write("## Params files\n")

				for r2,d2,f2 in os.walk(os.getcwd()):
					for file in f2:
						if file.endswith(".params"):
							restype_file.write(file + "\n")

				restype_file.close()
