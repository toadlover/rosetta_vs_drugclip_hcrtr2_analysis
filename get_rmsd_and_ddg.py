#the purpose of this script is to get the distance of the predicted placements (by ligand center of mass, not atom-atom distance due to atoms not preserving atom identifiers)
#for each system, the script will return a list of each system ranked by distance from the drugclip placement, along with predicted ddg
#there will also be a summary sheet that lists the closest placement for each system for the closest of all returned rosetta placements, closest of the top 10 ddg, and the single best ddg

#imports
import os,sys

#preliminary summary file opening
summary_file = open("pocket1_best_rosetta_placements.csv", "w")
#header
summary_file.write("system,best_all,best_10,best_1\n")

#proximity counter for the number of systems within 2 angstroms, withing 2-5 angstroms and beyond 5 angstroms for the closest placements for ddg cutoffs
#a 3x3 array to note counts by system for the closest placement to drugclip in the 3 energy fliters
proximity_counter = [[0,0,0],[0,0,0],[0,0,0]]

#root location
root_location = os.getcwd()

#now, iterate over each ligand folder to get the placement and comparison data to report
for r,d,f in os.walk(root_location):
	for dire in d:
		if r == root_location + "/input_data/pocket1":
			#we are in a directory to work with
			print(dire)

			dc_file = ""

			#get the ligand atom coordinates from the drugclip placement
			#find the corresponding cif file
			for r2,d2,f2 in os.walk(root_location + "/input_data/pocket1"):
				for file in f2:
					if file.endswith(dire + ".cif"):
						dc_file = r + "/" + file
						break

			print(dc_file)

			#read the dc file and get the coordinates
			orig_file = open(dc_file,"r")

			#open the file and make a list to track the atoms to get the center of mass
			#atom coords will be a list of 3 entry lists for the xyz coordinates, which will be averaged to get the approximate center of mass
			atom_coords = []

			for line in orig_file.readlines():
				if line.startswith("HETATM"):
					#get the coordinates
					atom_coords.append([float(line.split()[10]),float(line.split()[11]),float(line.split()[12])])

			#derive the drugclip center of mass
			dc_x = 0
			dc_y = 0
			dc_z = 0

			#get sum
			for atom in atom_coords:
				dc_x = dc_x + atom[0]
				dc_y = dc_y + atom[1]
				dc_z = dc_z + atom[2]

			#divide by number of atoms (average)
			dc_x = dc_x / len(atom_coords)
			dc_y = dc_y / len(atom_coords)
			dc_z = dc_z / len(atom_coords)

			#note, this also does not actually account for atom mass (which would only really be slightly shifted with sulphur or heavy halogens)
			dc_com = [dc_x,dc_y,dc_z]
			print(dc_com)

			#start a list to track the placements; will spit out to a file after we get all the palcement data and sort by ddg
			#list the placements by the file, ddg, com distance from dc
			placements_data = []

			#now, iterate over all placements from rosetta for this ligand
			for r2,d2,f2 in os.walk(root_location + "/input_data/pocket1/" + dire):
				for file in f2:
					if file.startswith("4s0v_receptor_only_dc_aligned_") and dire in file and file.endswith(".pdb"):
						#set up values to hold the ddg and atom coordinates
						ddg = 0

						r_atom_coords = []

						#read the file
						placement_file = open(r2 + "/" + file, "r")
						for line in placement_file.readlines():
							#append coordinates like with the drugclip placement from the cif, skip hydrogens though
							line_stripped = line.strip()

							#for ligand atoms
							if line_stripped.startswith("HETATM"):
								#skip hydrogens
								if line_stripped.endswith("H"):
									continue

								r_atom_coords.append([float(line_stripped.split()[6]),float(line_stripped.split()[7]),float(line_stripped.split()[8])])

							#grab the ddg
							if line_stripped.startswith("Scoring: Post-HighResDock system ddG:"):
								ddg = float(line_stripped.split()[len(line_stripped.split()) - 1])

						#derive the com
						r_x = 0
						r_y = 0
						r_z = 0

						#get sum
						for atom in r_atom_coords:
							r_x = r_x + atom[0]
							r_y = r_y + atom[1]
							r_z = r_z + atom[2]

						#divide by number of atoms (average)
						r_x = r_x / len(r_atom_coords)
						r_y = r_y / len(r_atom_coords)
						r_z = r_z / len(r_atom_coords)

						#note, this also does not actually account for atom mass (which would only really be slightly shifted with sulphur or heavy halogens)
						r_com = [r_x,r_y,r_z]

						#get the distance between the coms
						distance = ((r_com[0] - dc_com[0])**2 + (r_com[1] - dc_com[1])**2 + (r_com[2] - dc_com[2])**2) ** 0.5

						#append placement data
						placements_data.append([file,ddg,distance])

			#write the sorted placements to a file
			placements_data_file = open(root_location + "/input_data/pocket1/" + dire + "/" + dire +"_placements_data.csv","w")
			#write a header
			placements_data_file.write("file,ddg,distance\n")


			#behavior for if there were no placements
			if len(placements_data) == 0:
				summary_file.write(dire + ",n.a.,n.a.,n.a.\n")
				continue


			#sort the placements data to have the lowest ddg first
			sorted_placements = sorted(placements_data, key=lambda x: x[1])


			for placement in sorted_placements:
				placements_data_file.write(str(sorted_placements[0]) + "," + str(sorted_placements[1]) + "," + str(sorted_placements[2]) + "\n")

			#determine the best placement by distance out of all, top 10 ddg, and top ddg
			best_dist_all = 100
			best_dist_10 = 100
			best_dist_1 = 100

			#find the best distances out of all palcements
			for i in range(len(sorted_placements)):
				#best ddg
				if i == 0:
					best_dist_1 = float(sorted_placements[i][2])
				#best 10 ddg
				if i <= 9:
					if float(sorted_placements[i][2]) < best_dist_10:
						best_dist_10 = float(sorted_placements[i][2])

				#best all
				if float(sorted_placements[i][2]) < best_dist_all:
					best_dist_all = float(sorted_placements[i][2])		
					
			#write to the file
			summary_file.write(dire + "," + str(best_dist_all) + "," + str(best_dist_10) + "," + str(best_dist_1) + "\n")

			#adjust the proximity counter

			if best_dist_all <= 2:
				proximity_counter[0][0] = proximity_counter[0][0] + 1
			if best_dist_all > 2 and best_dist_all <= 5:
				proximity_counter[0][1] = proximity_counter[0][1] + 1
			if best_dist_all > 5:
				proximity_counter[0][2] = proximity_counter[0][2] + 1

			if best_dist_10 <= 2:
				proximity_counter[1][0] = proximity_counter[1][0] + 1
			if best_dist_10 > 2 and best_dist_10 <= 5:
				proximity_counter[1][1] = proximity_counter[1][1] + 1
			if best_dist_10 > 5:
				proximity_counter[1][2] = proximity_counter[1][2] + 1

			if best_dist_1 <= 2:
				proximity_counter[2][0] = proximity_counter[2][0] + 1
			if best_dist_1 > 2 and best_dist_1 <= 5:
				proximity_counter[2][1] = proximity_counter[2][1] + 1
			if best_dist_1 > 5:
				proximity_counter[2][2] = proximity_counter[2][2] + 1

#print the proximities
proximity_file = open("closest_placements_summary.csv", "w")
proximity_file.write(",0-2A,2-5A,>5A\n")
proximity_file.write("best_all," + str(proximity_counter[0][0]) + "," + str(proximity_counter[0][1]) + "," + str(proximity_counter[0][2]) + "\n")
proximity_file.write("best_10," + str(proximity_counter[1][0]) + "," + str(proximity_counter[1][1]) + "," + str(proximity_counter[1][2]) + "\n")
proximity_file.write("best_1," + str(proximity_counter[2][0]) + "," + str(proximity_counter[2][1]) + "," + str(proximity_counter[2][2]) + "\n")