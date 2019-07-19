# load required packages
library(tidyverse)

# required paramaters
input_beds_directory <- "<path>/03_peak_calling_v1beds/annotated_relaxed_peaks"
input_beds_suffix <- "rearrangedCols.txt"
output_directory <- "/<path>/03_peak_calling_v1beds/annotated_relaxed_peaks"
output_script_name <- "04_cut_n_run_meme_motifs_v1beds.sh"
genome_sequence <- "<path>/hg38_masked.fa"

# make working directory for intermediate files
dir.create(file.path(output_directory, "working_files"), showWarnings = FALSE)
setwd(file.path(output_directory, "working_files"))

# make list of the input files to load
input_beds <- list.files( path = input_beds_directory, pattern = input_beds_suffix,
                          full.names = TRUE )

# for each file in the input list
# current_file_num <- 1
for(current_file_num in 1:length(input_beds)) {

    # pull out the current files name for saving the output file
    current_file_name <- str_split(input_beds[current_file_num], "/")[[1]]
    current_file_name <- current_file_name[length(current_file_name)]
    current_file_name <- str_split(current_file_name, "\\.")[[1]][1]
    
    # read in current file
    current_file <- read_delim( input_beds[current_file_num], 
                               "\t", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE)
    
    # update the user whats going on
    message(paste0("processing - ", current_file_name))
    
    # label the columns
    colnames(current_file) <- c("chr", "start", "end", "peak_id", "max_signal", "total_signal")
    
    # determine midpoint
    current_file$midpoint <- round((current_file$end - current_file$start)/2, digits = 0)
    # reset start position
    current_file$start <- current_file$start + current_file$midpoint
    # make the end position of the summit the position after
    current_file$end <- current_file$start + 1
    
    # drop columns no longer needed
    current_file <- current_file[,1:5]
    
    # save the peak summit files
    write.table( current_file[sample(nrow(current_file), round(nrow(current_file)/4)), ], file = paste0(current_file_name, ".peakSummits.bed"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)    
}

# for all the summit files, run bedops on the coordintes and pad with 150bp each side!
summit_files <- list.files( path = getwd(), pattern = ".peakSummits.bed",
                            full.names = TRUE )

# start output bash script to run
sink( file = "04-1_bedops_rangePad.sh", append = FALSE, type = "output", split = FALSE )
cat(paste("#!/bin/bash"))
cat(paste("\n"))

# for each file in the input list
# current_file_num <- 1
for(current_file_num in 1:length(summit_files)) {
    
    # pull out the current files name for saving the output file
    current_file_name <- str_split(summit_files[current_file_num], "/")[[1]]
    current_file_name <- current_file_name[length(current_file_name)]
    current_file_name <- str_split(current_file_name, "\\.")[[1]][1]
    
    # run bedops
    bedops_command <- paste0("bedops --range 150 -u ", summit_files[current_file_num], " > ./", current_file_name, ".paddedSummit.bed")
    cat(bedops_command)
    cat(paste("\n"))
}
sink()

# run the bedops scrips
system("bash ./04-1_bedops_rangePad.sh")


# for all the padded summit files, run getFasta to generate file for meme
padded_files <- list.files( path = getwd(), pattern = ".paddedSummit.bed",
                            full.names = TRUE )

# start output bash script to run
sink( file = "04-2_bedtools_getFasta.sh", append = FALSE, type = "output", split = FALSE )
cat(paste("#!/bin/bash"))
cat(paste("\n"))

# for each file in the input list
# current_file_num <- 1
for(current_file_num in 1:length(padded_files)) {
    
    # pull out the current files name for saving the output file
    current_file_name <- str_split(padded_files[current_file_num], "/")[[1]]
    current_file_name <- current_file_name[length(current_file_name)]
    current_file_name <- str_split(current_file_name, "\\.")[[1]][1]
    
    # run bedops
    bedtools_command <- paste0("bedtools getfasta -fi ", genome_sequence, " -bed ", padded_files[current_file_num], " -fo ../", current_file_name, ".peakRegions.forMeme.fa")
    cat(bedtools_command)
    cat(paste("\n"))
}
sink()

# run the bedtools script
system("bash ./04-2_bedtools_getFasta.sh")
