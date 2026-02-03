#the purpose of this script is to go through each placement file that was made through rosetta discovery, and pull up to the top 3 placements per system
#the script will pull the top 3 systems with the best ddg that also have a real motifs ratio of at least 0.25 (25%)
#all files that pass will be moved to a new directory and listed in a spreadsheet in alphabetical order by file placement name, which can be utilized for viewing

#imports
import os,sys

#make a directory to copy the top placements into
os.system("mkdir top_pocket1_placements")

#start a list to hold all of the top entries per system, which we will sort at the end and then write to a file
top_placements_all = []

#look at all sorted placement summary csvs to collect relevant data for each placement
for r,d,f in os.walk(os.getcwd()):
	for file in f:

		#a counter to determine if the number of placements kept for a ligand has exceeded 3, will reset to 0 here
		placements_kept = 0

		if file.endswith("_placements_data.csv") and file.split("_")[0] in root and "input_data/pocket1" in root:
			#read the file and get up to the top 3
			read_file = open(r + "/" + file,"r")
			for line in read_file.readlines():
				#skip the header
				if line.startswith("file,ddg,distance,real_motif_ratio"):
					continue

				#otherwise, break up the line and interpret if the placement should be kept
				line_stripped = line.strip()

				#file name
				placement_file = line_stripped.split(",")[0]
				#system free energy
				ddg = line_stripped.split(",")[1]
				#placement distance from drugclip prediction
				distance = line_stripped.split(",")[2]
				#real motif ratio
				real_motif_ratio = float(line_stripped.split(",")[3])

				#if the real motif ratio is >= 0.25, keep the placement
				if real_motif_ratio >= 0.25:
					top_placements_all.append([placement_file,ddg,distance,real_motif_ratio])

					#increment the placements kept counter
					placements_kept = placements_kept + 1

					#copy the ligand placement to the top placements folder
					os.system("cp " + r + "/" + placement_file + " top_pocket1_placements")

					#if the counter is 3 (we got 3 placements), break and move to the next
					if placements_kept == 3:
						break

#we should now have all placements, make a list file of them, sorted in alphabetical order by file name
sorted_placements = sorted(top_placements_all, key=lambda x: x[0])

#write to file
write_file = open("top_pocket1_placements/top_placements_list.csv", "w")

write_file.write("top_placements_all\n")

for placement in sorted_placements:
	write_file.write(str(placement[0]) + "," + str(placement[1]) + "," + str(placement[2]) + "," + str(placement[3]) + "\n")