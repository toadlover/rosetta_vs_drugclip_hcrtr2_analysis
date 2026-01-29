#the purpose of this script is to prepare rosetta arg files and submit rosetta jobs in LSF, only for the pocket1 (main HCRTR2 pocket, we are ignoring the other 2 pockets)

#imports
import os,sys

#root location
root_location = os.getcwd()

#move into the input_data/pocket1 location
os.chdir("input_data/pocket1")

#pocket root
pocket_root = os.getcwd()

#for the input files/directories, they need to be mapped to a container location for the Rosetta container (by default, will map to the /input location for the singularity call)
#here we will make strings for the mapping, which will go in the rosetta args file (the actual paths will be used in executing in the singularity image)
input_target_pdb = "/input/4s0v_receptor_only_dc_aligned.pdb"
input_motifs_file = "/input/FINAL_motifs_list_filtered_2_3_2023.motifs"
input_test_params_dir = "/input/test_params/"


target_pdb = "/home/ari.ginsparg-umw/rosetta_vs_drugclip_hcrtr2_analysis/input_data/4s0v_receptor_only_dc_aligned.pdb"
motifs_file = "/pi/summer.thyme-umw/enamine-REAL-2.6billion/FINAL_motifs_list_filtered_2_3_2023.motifs"

#iterate over each directory, and prepare to run discovery on that set of conformers for the ligand
for r,d,f in os.walk(os.getcwd()):
	for dire in d:
		#if the dire is at the r level, work with it
		if r == os.getcwd():
			
			#temporary test limit for only running on dire Z1840327359
			if dire != "Z1840327359":
				continue

			os.chdir(dire)

			test_params_dir = r + "/" + dire + "/test_params/"

			#write an arg file
			args_file = open("args","w")

			#these first few args are hard coded for housekeeping, but can be removed later if we really want these to be mutable/not included
			#if someone really wants to change these, they can include the arg again in the extra args file, and the arg will get overwritten with the desired value
			args_file.write("#keep seed constant\n")
			args_file.write("-constant_seed 1\n")
			args_file.write("#ignore unrecognized residues to help mitigate crashes\n")
			args_file.write("-ignore_unrecognized_res\n")
			args_file.write("#handle ligand repeats if using multiple anchor residues, will otherwise crash without this flag\n")
			args_file.write("-in::file::override_database_params true\n")
			args_file.write("#constrain coordinates\n")
			args_file.write("-constrain_relax_to_start_coords\n")
			args_file.write("#keep all placements; 0 means keep all, any other integer means keep up to that integer\n")
			args_file.write("-best_pdbs_to_keep 0\n")

			#user input dependent
			args_file.write("#mapped protein system\n")
			#/home/ari.ginsparg-umw/rosetta_vs_drugclip_hcrtr2_analysis/input_data/4s0v_receptor_only_dc_aligned.pdb
			args_file.write("-s " + input_target_pdb + "\n")
			args_file.write("#mapped motifs file\n")
			args_file.write("-motif_filename " + input_motifs_file + "\n")
			args_file.write("#mapped test_params directory\n")
			args_file.write("-params_directory_path " + input_test_params_dir + "\n")
			args_file.write("#rosetta-indexed anchor residue index/indices\n")
			args_file.write("-protein_discovery_locus 500 \n")
			args_file.write("#fa_atr cutoff\n")
			args_file.write("-fa_atr_cutoff = -2 \n")
			args_file.write("#fa_rep cutoff\n")
			args_file.write("-fa_rep_cutoff = 150 \n")
			args_file.write("#ddg cutoff\n")
			args_file.write("-ddg_cutoff = -9 \n")

			#close the file
			args_file.close()

			#run Rosetta
			#we now have the args file written, now call Rosetta discovery
			os.system("singularity exec --bind " + test_params_dir + ":" + input_test_params_dir + " --bind " + os.getcwd() + "/args:/input/args --bind " + target_pdb + ":" + input_target_pdb + " --bind " + motifs_file + ":" + input_motifs_file + " /pi/summer.thyme-umw/enamine-REAL-2.6billion/rosetta_condensed_6_25_2024.sif /rosetta/source/bin/ligand_discovery_search_protocol.linuxgccrelease @/input/args")

			#move all pdb files to a placements directory
			os.system("mkdir placements")

			os.system("mv *pdb placements")

			os.chdir("placements")

			#rename each pdb file by prepending the anchor residue string used by this script
			for r,d,f in os.walk(os.getcwd()):
				for file in f:
					if file.endswith(".pdb") and r == os.getcwd():
						os.system("mv " + file + " res" + anchor_residue_string + "_" + file)

			#now, call the placement analysis script
			os.system("python /pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/rosetta/score_placed_ligands_with_filtering.py")

			#copy the csv files up a level for easy accession outside of the to-be compressed placements directory
			os.system("cp *csv ..")

			#then, compress the placement files (this will move all pdb files to a directory called placements, so do not keep any important pdbs in here)
			os.chdir("..")



			#at end, return to the pocket1 root
			os.chdir(pocket_root)