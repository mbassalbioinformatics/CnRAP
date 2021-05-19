# import packages required
import os
from os import listdir
import sys
import numpy as np
import pandas as pd
from tabulate import tabulate # python3 -m pip install tabulate

# define function to round scaling factors - coz there isnt a built in way?
def truncate(number, decimals=0):
    multiplier = 10 ** decimals
    return int(number * multiplier) / multiplier

# read in input arguments that are required
aligned_bams_folder		= sys.argv[1] # where the aligned bams are saved <path>/01_read_trimming_alignment/
normalized_beds_folder	= sys.argv[2] # where to dump the pocessed bed files which ar ready for SEACR <path>/02_bam_processing_v1
chrom_sizes_txt			= sys.argv[3] # path and name of the chormosome sizes text file needed for the bed conversion <path>/hg38.chrom.sizes

# folders/programs to be used for initializing the script
bamToBed_folder			=	normalized_beds_folder+ "/bed_converted_bams/"
logs_folder				=	normalized_beds_folder+ "/logs/"

# create output directories if don't exists
if not os.path.exists(normalized_beds_folder):
	os.makedirs(normalized_beds_folder)

if not os.path.exists(bamToBed_folder):
	os.makedirs(bamToBed_folder)

if not os.path.exists(logs_folder):
	os.makedirs(logs_folder)

# Move to the folder with the saved bam and stats files
os.chdir( aligned_bams_folder )

# figure out which are the mapped stats files to pull out the number of mapped reads
mapped_stats_files = []
for file_name in listdir():
	if file_name.endswith("mappedStats.txt"):
		mapped_stats_files.append(file_name)
# sort them alphabetically
mapped_stats_files.sort()
# reshape so each SC bam name is matched row to hg bam - so can iterate through and hvae everythign in order
#	that matches the order of the mapped bam files below
mapped_stats_files = np.reshape( mapped_stats_files, (-1,2) )

# figure out which are the bam files to convert and normalize - need full paths
mapped_sorted_bam_files = []
for file_name in listdir():
	if file_name.endswith("hgAligned_stampyFiltered.mappedSorted.bam"):
		mapped_sorted_bam_files.append(aligned_bams_folder+ "/" +file_name)
# sort alphabetically
mapped_sorted_bam_files.sort()

# pull out all the mapped reads for each file and normlaize hg to sc mapped number of reads
mapped_reads_per_file = [[0 for x in range(3)] for y in range(len(mapped_stats_files))]
# for every row of the stats dataframe
for sample_counter in range(len(mapped_stats_files)):
	# read in the hg file and skip all the header rows
	hg_file_reads = pd.read_csv(mapped_stats_files[sample_counter,0], sep='\t', header=None, skiprows=12)
	# split the 7th row by whitespace to pull out just the number of reads and save
	mapped_reads_per_file[sample_counter][0] = hg_file_reads[0][0].split()[1]
	# read in the hg file and skip all the header rows
	sc_file_reads = pd.read_csv(mapped_stats_files[sample_counter,1], sep='\t', header=None, skiprows=12)
	# split the 7th row by whitespace to pull out just the number of reads and save
	mapped_reads_per_file[sample_counter][1] = sc_file_reads[0][0].split()[1]
	# calculate scaling factor which will be used for scaling all the counts at each bp position when normalizing
	# based on https://github.com/Henikoff/Cut-and-Run/blob/master/spike_in_calibration.sh
	# scale * (primary_genome_mapped_count_at_bp)/(spike-in_genome_total_of_mapped_fragments)
	# since paired end, # fragments = # mapped reads / 2
	mapped_reads_per_file[sample_counter][2] = 10000000 / ( int(mapped_reads_per_file[sample_counter][0]) / 2 )


# Move to the folder with the saved bam and stats files
os.chdir( bamToBed_folder )
# write the file of the mapped reads and normalization factors for future reference
reads_norm_factor_file = open('999_mapped_num_reads_per_file_norm_factor.txt', 'w')
reads_norm_factor_file.write(tabulate(mapped_stats_files))
reads_norm_factor_file.write(tabulate(mapped_reads_per_file))
reads_norm_factor_file.close()
# now have a listing of all the human bam files and also a listing of the normalization factors for each - now to make output script and do the work!

# basic function of the script is to do the following and prepr the bams for input into SEACR
# bedtools bamtobed -i <...> -bedpe > $output_dir/"$output_name"_properpair.bed
# bedtools sort -i $output_dir/"$output_name"_properpair.bed > $output_dir/"$output_name"_properpair_sortchr.bed
# awk '{ print $1"\t"$2"\t"$6}' $output_dir/"$output_name"_properpair_sortchr.bed > $output_dir/"$output_name"_properpair_sortchr_startend.bed
# bedtools genomecov -bg -i $output_dir/"$output_name"_properpair_sortchr_startend.bed -g </mm10.chrom.sizes> > $output_dir/"$output_name"_properpair_sortchr_startend.bg

# start saving output script
os.chdir( normalized_beds_folder )
script_name = "02_cut_n_run_bamToBed_normalize.sh"
output_script = open( script_name, 'w' )

# convert each bam to a bed file as recommended by https://github.com/FredHutch/SEACR
for sample_counter in range(len(mapped_sorted_bam_files)):
	output_command = "echo \"converting and normlaizing : " +mapped_sorted_bam_files[sample_counter].split("/")[-1][:-4]+ "\""
	output_script.write(output_command)
	output_script.write("\n")
	# convert bam to bed using bedpe flag - pull out only proper paired reads and then create bedgraph of everything so as to allow scaling of reads!
	output_command = "samtools view -b -f 2 -F 524 " +mapped_sorted_bam_files[sample_counter]+ " | bedtools genomecov -bg -scale " +str( truncate(mapped_reads_per_file[sample_counter][2], 2) )+ " -ibam stdin > " +normalized_beds_folder+ "/" +mapped_sorted_bam_files[sample_counter].split("/")[-1][:-3]+ "normalized.bed 2> " +logs_folder+ "" +mapped_sorted_bam_files[sample_counter].split("/")[-1][:-3]+ "normalize.err"
	output_script.write(output_command)
	
	output_script.write("\n\n\n")

# close the script file
output_script.close()

# make script executable
os.system("chmod +x " +script_name)
