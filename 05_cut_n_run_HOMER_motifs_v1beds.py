# import packages required
import os
from os import listdir
import sys

# read in input arguments that are required
peaks_dir				= sys.argv[1] # where is seacr.py located? <path>/03_peak_calling_v1beds
homer_beds_extension	= sys.argv[2] # where are the normalized bed_files homerMotif.bed
output_base_dir			= sys.argv[3] # base directory where to dump the identified motifs <path>/04_peak_motifs_v1beds

# Move to the folder with the saved bam and stats files
os.chdir( peaks_dir )

# figure out which are the nuclear files, which are the regular files
peak_beds = []
for file_name in listdir():
	if file_name.endswith(homer_beds_extension):
		peak_beds.append(peaks_dir+ "/" +file_name)
# sort alphabetically
peak_beds.sort()

# folders/programs to be used for initializing the script
output_directory		=	output_base_dir+ "/homer_motifs/"

# create output directories if don't exists
if not os.path.exists(output_directory):
	os.makedirs(output_directory)

# Move to the output folder
os.chdir( output_directory )

# save commands to the output script
script_name = "04_cut_n_run_homer_motifs_v1beds.sh"
output_script = open( script_name, 'w' )

# for each file re-arrange to make in the format for Homer to call peaks
for sample_counter in range(len(peak_beds)):
	# update the user as to whats going on
	output_command = "echo \"Running Motif Discovery on v1 " +peak_beds[sample_counter].split("/")[8][:-35]+ "\""
	output_script.write(output_command)
	output_script.write("\n")
	# call Homer motif discovery
	output_command = "findMotifsGenome.pl " +peak_beds[sample_counter]+ " hg38 ./" +peak_beds[sample_counter].split("/")[8][:-35]+ "_given -size given -mask -p 20 -S 50"
	output_script.write(output_command)
	output_script.write("\n")
	output_command = "findMotifsGenome.pl " +peak_beds[sample_counter]+ " hg38 ./" +peak_beds[sample_counter].split("/")[8][:-35]+ " -size 50,100,200 -mask -p 20 -S 50"
	output_script.write(output_command)
	output_script.write("\n\n\n")

# close the script file
output_script.close()

# make script executable
os.system("chmod +x " +script_name)