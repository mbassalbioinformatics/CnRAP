# import packages required
import os
import sys

# read in input arguments that are required
sample_id		= sys.argv[1] # what to name the samples - easier to define than to subset
read1_fq_gz		= sys.argv[2] # read1.fq.gz
read2_fq_gz		= sys.argv[3] # read2.fq.gz
num_cores		= sys.argv[4] # how many cores to use for the analysis
aligned_folder	= sys.argv[5] # where to dump the aligned reads

# folders/programs to be used for initializing the script
trimmomatic_path_jar	=	"<path>/trimmomatic-0.36.jar"
adapter_path			=	"<path>/zzz_adapters/"
kseq_path_exec			=	"<path>/kseq_test"
bwa_path_exec			=	"<path>/bwa"
bwa_hg_masked_index		=	"/<path>/human_bwa_index_hg38_masked"
bwa_sc_masked_index		=	"<path>/sacCer_bwa_index"
stampy_hg_stidx_hash	=	"<path>/human_bwa_index_hg38_masked"
stampy_sc_stidx_hash	=	"<path>/sacCer_bwa_index"
samtools_path_exec		=	"<path>/samtools"
picard_path_jar			=	"<path>/picard.jar"
read_length_file		=	"<path>/zzz_read_length/length"
multimap_limit			=	str( 20 )

# pull out how long the reads are from the file- should be 42 though
with open(read_length_file, 'r') as file:
	read_length = str( int( file.read().strip() ) )

# create output directories if don't exist
if not os.path.exists(aligned_folder):
	os.makedirs(aligned_folder)

trim_folder	=	"" +aligned_folder+ "/trimmed_reads/"
if not os.path.exists(trim_folder):
	os.makedirs(trim_folder)

intermediates_folder	=	"" +aligned_folder+ "/intermediate_files/"
if not os.path.exists(intermediates_folder):
	os.makedirs(intermediates_folder)

unapired_trim_folder	=	"" +trim_folder+ "unpaired/"
if not os.path.exists(unapired_trim_folder):
	os.makedirs(unapired_trim_folder)

logs_folder = "" +trim_folder+ "logs/"
if not os.path.exists(logs_folder):
	os.makedirs(logs_folder)


# start actually doing work now
os.chdir( aligned_folder )
# start saving output script
script_name = "01_cut_n_run_" +sample_id+ "_trim_align.sh"
output_script = open( script_name, 'w' )

# make output command write to file
output_command = "echo \"Trimmomatic : " +sample_id+ "\""
output_script.write(output_command)
output_script.write("\n")


# trim the files with trimmomatic
# java -jar $trimmomaticbin/$trimmomaticjarfile PE -threads 1 -phred33 $dirname/"$base"_R1_001.fastq.gz $dirname/"$base"_R2_001.fastq.gz $trimdir/"$base"_1.paired.fastq.gz $trimdir/"$base"_1.unpaired.fastq.gz $trimdir/"$base"_2.paired.fastq.gz $trimdir/"$base"_2.unpaired.fastq.gz ILLUMINACLIP:$adapterpath/Truseq3.PE.fa:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25
output_command = "java -jar " +trimmomatic_path_jar+ " PE -threads " +num_cores+ " -phred33 " +read1_fq_gz+ " " +read2_fq_gz+ " " +trim_folder+ "" +sample_id+ "_1.paired_trimmomatic.fastq.gz " +unapired_trim_folder+ "" +sample_id+ "_1.unpaired.fastq.gz " +trim_folder+ "" +sample_id+ "_2.paired_trimmomatic.fastq.gz " +unapired_trim_folder+ "" +sample_id+ "_2.unpaired.fastq.gz ILLUMINACLIP:" +adapter_path+ "Truseq3.PE.fa:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25 2> " +logs_folder+ "" +sample_id+ "_trimmomatic.err"
output_script.write(output_command)
output_script.write("\n\n")



# make output command write to file
output_command = "echo \"kseq Trimming : " +sample_id+ "\""
output_script.write(output_command)
output_script.write("\n")

# trim the files with kseq for both reads
# $kseqbin/kseq_test $trimdir/"$base"_1.paired.fastq.gz $len $trimdir2/"$base"_1.paired.fastq.gz
output_command = "" +kseq_path_exec+ " " +trim_folder+ "" +sample_id+ "_1.paired_trimmomatic.fastq.gz " +read_length+ " " +trim_folder+ "" +sample_id+ "_1.paired_kseqTrimmed.fastq.gz 2> " +logs_folder+ "" +sample_id+ "_kseqTrim_r1.err"
output_script.write(output_command)
output_script.write("\n")

# $kseqbin/kseq_test $trimdir/"$base"_2.paired.fastq.gz $len $trimdir2/"$base"_2.paired.fastq.gz
output_command = "" +kseq_path_exec+ " " +trim_folder+ "" +sample_id+ "_2.paired_trimmomatic.fastq.gz " +read_length+ " " +trim_folder+ "" +sample_id+ "_2.paired_kseqTrimmed.fastq.gz 2> " +logs_folder+ "" +sample_id+ "_kseqTrim_r2.err"
output_script.write(output_command)
output_script.write("\n\n")



# make output command write to file
output_command = "echo \"Alignning using BWA and Stampy hg : " +sample_id+ "\""
output_script.write(output_command)
output_script.write("\n")

# align the reads with bwa
# bwa aln -q10 -t8 hg18 reads_1.fastq > 1.sai
output_command = "" +bwa_path_exec+ " aln -t " +num_cores+ " " +bwa_hg_masked_index+ " " +trim_folder+ "" +sample_id+ "_1.paired_kseqTrimmed.fastq.gz > " +aligned_folder+ "/" +sample_id+ "_hgAligned_r1.sai 2> " +logs_folder+ "" +sample_id+ "_bwa_hg_r1sai.err"
output_script.write(output_command)
output_script.write("\n")

# bwa aln -q10 -t8 hg18 reads_2.fastq > 2.sai
output_command = "" +bwa_path_exec+ " aln -t " +num_cores+ " " +bwa_hg_masked_index+ " " +trim_folder+ "" +sample_id+ "_2.paired_kseqTrimmed.fastq.gz > " +aligned_folder+ "/" +sample_id+ "_hgAligned_r2.sai 2> " +logs_folder+ "" +sample_id+ "_bwa_hg_r2sai.err"
output_script.write(output_command)
output_script.write("\n")

# bwa sampe hg18 1.sai 2.sai reads_1.fastq reads_2.fastq | samtools view -Sb - > bwa.bam
output_command = "" +bwa_path_exec+ " sampe -n " +multimap_limit+ " " +bwa_hg_masked_index+ " " +aligned_folder+ "/" +sample_id+ "_hgAligned_r1.sai " +aligned_folder+ "/" +sample_id+ "_hgAligned_r2.sai " +trim_folder+ "" +sample_id+ "_1.paired_kseqTrimmed.fastq.gz " +trim_folder+ "" +sample_id+ "_2.paired_kseqTrimmed.fastq.gz | samtools view -Sb - > " +aligned_folder+ "/" +sample_id+ "_hgAligned.bam 2> " +logs_folder+ "" +sample_id+ "_bwa_hg_sampe.err"
output_script.write(output_command)
output_script.write("\n")

# ./stampy.py -g hg18 -h hg18 -t 8 --bamkeepgoodreads -M bwa.bam 
output_command = "stampy.py -g " +stampy_hg_stidx_hash+ " -h " +stampy_hg_stidx_hash+ " -t " +num_cores+ " --sensitive -M " +aligned_folder+ "/" +sample_id+ "_hgAligned.bam | samtools view -Sb - > " +aligned_folder+ "/" +sample_id+ "_hgAligned_stampyFiltered.bam 2> " +logs_folder+ "" +sample_id+ "_stampy_hg.err"
output_script.write(output_command)
output_script.write("\n\n")
# get the mapping stats since bwa doesnt provide
output_command = "bamtools stats -in " +aligned_folder+ "/" +sample_id+ "_hgAligned_stampyFiltered.bam > " +aligned_folder+ "/" +sample_id+ "_hgAligned_stampyFiltered.allStats.txt"
output_script.write(output_command)
output_script.write("\n\n")

# make output command write to file
output_command = "echo \"Removing Unmapped reads hg : " +sample_id+ "\""
output_script.write(output_command)
output_script.write("\n")
# remove unmapped reads
output_command = "bamtools filter -isMapped true -in " +aligned_folder+ "/" +sample_id+ "_hgAligned_stampyFiltered.bam > " +aligned_folder+ "/" +sample_id+ "_hgAligned_stampyFiltered.mapped.bam 2> " +logs_folder+ "" +sample_id+ "_bamtools_removeUnmapped_hg.err"
output_script.write(output_command)
output_script.write("\n")
# make output command write to file
output_command = "echo \"Sort, Index, Alignment stats hg : " +sample_id+ "\""
output_script.write(output_command)
output_script.write("\n")
# samtools sort -l 9 -O bam -@ " +num_cores+ " -o " +output_dir+ "" +sample_name+ ".sorted.bam " +output_dir+ "" +sample_name+ "Aligned.sortedByCoord.out.bam"
output_command = "samtools sort -l 9 -O bam -n -@ " +num_cores+ " -o " +aligned_folder+ "/" +sample_id+ "_hgAligned_stampyFiltered.mappedSorted.bam " +aligned_folder+ "/" +sample_id+ "_hgAligned_stampyFiltered.mapped.bam"
output_script.write(output_command)
output_script.write("\n")
# fixmates
output_command = "samtools fixmate " +aligned_folder+ "/" +sample_id+ "_hgAligned_stampyFiltered.mappedSorted.bam " +aligned_folder+ "/" +sample_id+ "_hgAligned_stampyFiltered.mappedSortedFixed.bam"
output_script.write(output_command)
output_script.write("\n")
# sort fixed mates
output_command = "samtools sort -l 9 -O bam -@ " +num_cores+ " -o " +aligned_folder+ "/" +sample_id+ "_hgAligned_stampyFiltered.mappedSorted.bam " +aligned_folder+ "/" +sample_id+ "_hgAligned_stampyFiltered.mappedSortedFixed.bam"
output_script.write(output_command)
output_script.write("\n")
# index
output_command = "samtools index -b " +aligned_folder+ "/" +sample_id+ "_hgAligned_stampyFiltered.mappedSorted.bam"
output_script.write(output_command)
output_script.write("\n")
# genome coverage map
output_command = "bamCoverage -b " +aligned_folder+ "/" +sample_id+ "_hgAligned_stampyFiltered.mappedSorted.bam -o " +aligned_folder+ "/" +sample_id+ "_hgAligned_stampyFiltered.mappedSorted.bw 2> " +logs_folder+ "" +sample_id+ "_bamCoverage.err"
output_script.write(output_command)
output_script.write("\n")
# sort fixed mates in way to convert to beds! - sort by names since thats whats required
output_command = "samtools sort -l 9 -O bam -n -@ " +num_cores+ " -n -o " +aligned_folder+ "/" +sample_id+ "_hgAligned_stampyFiltered.mappedSortedBed.bam " +aligned_folder+ "/" +sample_id+ "_hgAligned_stampyFiltered.mappedSortedFixed.bam"
output_script.write(output_command)
output_script.write("\n")
# move intermediate files
output_command = "mv " +sample_id+ "_hgAligned.bam " +intermediates_folder
output_script.write(output_command)
output_script.write("\n")
output_command = "mv " +sample_id+ "_hgAligned_stampyFiltered.bam " +intermediates_folder
output_script.write(output_command)
output_script.write("\n")
output_command = "mv " +sample_id+ "_hgAligned_stampyFiltered.mappedSortedFixed.bam " +intermediates_folder
output_script.write(output_command)
output_script.write("\n")
output_command = "mv " +sample_id+ "_hgAligned_stampyFiltered.mapped.bam " +intermediates_folder
output_script.write(output_command)
output_script.write("\n")
output_command = "mv " +sample_id+ "_hgAligned_r1.sai " +intermediates_folder
output_script.write(output_command)
output_script.write("\n")
output_command = "mv " +sample_id+ "_hgAligned_r2.sai " +intermediates_folder
output_script.write(output_command)
output_script.write("\n")
# get the mapping stats since bwa doesnt provide
output_command = "bamtools stats -in " +aligned_folder+ "/" +sample_id+ "_hgAligned_stampyFiltered.mappedSorted.bam > " +aligned_folder+ "/" +sample_id+ "_hgAligned_stampyFiltered.mappedStats.txt"
output_script.write(output_command)
output_script.write("\n\n")


# make output command write to file
output_command = "echo \"Alignning using BWA and Stampy sc : " +sample_id+ "\""
output_script.write(output_command)
output_script.write("\n")

# align the reads with bwa
# bwa aln -q10 -t8 hg18 reads_1.fastq > 1.sai
output_command = "" +bwa_path_exec+ " aln -t " +num_cores+ " " +bwa_sc_masked_index+ " " +trim_folder+ "" +sample_id+ "_1.paired_kseqTrimmed.fastq.gz > " +aligned_folder+ "/" +sample_id+ "_scAligned_r1.sai 2> " +logs_folder+ "" +sample_id+ "_bwa_sc_r1sai.err"
output_script.write(output_command)
output_script.write("\n")

# bwa aln -q10 -t8 hg18 reads_2.fastq > 2.sai
output_command = "" +bwa_path_exec+ " aln -t " +num_cores+ " " +bwa_sc_masked_index+ " " +trim_folder+ "" +sample_id+ "_2.paired_kseqTrimmed.fastq.gz > " +aligned_folder+ "/" +sample_id+ "_scAligned_r2.sai 2> " +logs_folder+ "" +sample_id+ "_bwa_sc_r2sai.err"
output_script.write(output_command)
output_script.write("\n")

# bwa sampe hg18 1.sai 2.sai reads_1.fastq reads_2.fastq | samtools view -Sb - > bwa.bam
output_command = "" +bwa_path_exec+ " sampe -n " +multimap_limit+ " " +bwa_sc_masked_index+ " " +aligned_folder+ "/" +sample_id+ "_scAligned_r1.sai " +aligned_folder+ "/" +sample_id+ "_scAligned_r2.sai " +trim_folder+ "" +sample_id+ "_1.paired_kseqTrimmed.fastq.gz " +trim_folder+ "" +sample_id+ "_2.paired_kseqTrimmed.fastq.gz | samtools view -Sb - > " +aligned_folder+ "/" +sample_id+ "_scAligned.bam 2> " +logs_folder+ "" +sample_id+ "_bwa_sc_sampe.err"
output_script.write(output_command)
output_script.write("\n")

# ./stampy.py -g hg18 -h hg18 -t 8 --bamkeepgoodreads -M bwa.bam 
output_command = "stampy.py -g " +stampy_sc_stidx_hash+ " -h " +stampy_sc_stidx_hash+ " -t " +num_cores+ " --sensitive -M " +aligned_folder+ "/" +sample_id+ "_scAligned.bam | samtools view -Sb - > " +aligned_folder+ "/" +sample_id+ "_scAligned_stampyFiltered.bam 2> " +logs_folder+ "" +sample_id+ "_stampy_sc.err"
output_script.write(output_command)
output_script.write("\n\n")

# get the mapping stats since bwa doesnt provide
output_command = "bamtools stats -in " +aligned_folder+ "/" +sample_id+ "_scAligned_stampyFiltered.bam > " +aligned_folder+ "/" +sample_id+ "_scAligned_stampyFiltered.allStats.txt"
output_script.write(output_command)
output_script.write("\n\n")

# make output command write to file
output_command = "echo \"Removing Unmapped reads sc : " +sample_id+ "\""
output_script.write(output_command)
output_script.write("\n")
# remove unmapped reads
output_command = "bamtools filter -isMapped true -in " +aligned_folder+ "/" +sample_id+ "_scAligned_stampyFiltered.bam > " +aligned_folder+ "/" +sample_id+ "_scAligned_stampyFiltered.mapped.bam 2> " +logs_folder+ "" +sample_id+ "_bamtools_removeUnmapped_sc.err"
output_script.write(output_command)
output_script.write("\n")
# make output command write to file
output_command = "echo \"Sort, Index, Alignment stats hg : " +sample_id+ "\""
output_script.write(output_command)
output_script.write("\n")
# samtools sort -l 9 -O bam -@ " +num_cores+ " -o " +output_dir+ "" +sample_name+ ".sorted.bam " +output_dir+ "" +sample_name+ "Aligned.sortedByCoord.out.bam"
output_command = "samtools sort -l 9 -O bam -n -@ " +num_cores+ " -o " +aligned_folder+ "/" +sample_id+ "_scAligned_stampyFiltered.mappedSorted.bam " +aligned_folder+ "/" +sample_id+ "_scAligned_stampyFiltered.mapped.bam"
output_script.write(output_command)
output_script.write("\n")
# fixmates
output_command = "samtools fixmate " +aligned_folder+ "/" +sample_id+ "_scAligned_stampyFiltered.mappedSorted.bam " +aligned_folder+ "/" +sample_id+ "_scAligned_stampyFiltered.mappedSortedFixed.bam"
output_script.write(output_command)
output_script.write("\n")
# sort fixed mates
output_command = "samtools sort -l 9 -O bam -@ " +num_cores+ " -o " +aligned_folder+ "/" +sample_id+ "_scAligned_stampyFiltered.mappedSorted.bam " +aligned_folder+ "/" +sample_id+ "_scAligned_stampyFiltered.mappedSortedFixed.bam"
output_script.write(output_command)
output_script.write("\n")
# index
output_command = "samtools index -b " +aligned_folder+ "/" +sample_id+ "_scAligned_stampyFiltered.mappedSorted.bam"
output_script.write(output_command)
output_script.write("\n")
# genome coverage map
output_command = "bamCoverage -b " +aligned_folder+ "/" +sample_id+ "_scAligned_stampyFiltered.mappedSorted.bam -o " +aligned_folder+ "/" +sample_id+ "_scAligned_stampyFiltered.mappedSorted.bw 2> " +logs_folder+ "" +sample_id+ "_bamCoverage.err"
output_script.write(output_command)
output_script.write("\n")
# sort fixed mates in way to convert to beds! - sort by names since thats whats required
output_command = "samtools sort -l 9 -O bam -n -@ " +num_cores+ " -n -o " +aligned_folder+ "/" +sample_id+ "_scAligned_stampyFiltered.mappedSortedBed.bam " +aligned_folder+ "/" +sample_id+ "_scAligned_stampyFiltered.mappedSortedFixed.bam"
output_script.write(output_command)
output_script.write("\n")
# move intermediate files
output_command = "mv " +sample_id+ "_scAligned.bam " +intermediates_folder
output_script.write(output_command)
output_script.write("\n")
output_command = "mv " +sample_id+ "_scAligned_stampyFiltered.bam " +intermediates_folder
output_script.write(output_command)
output_script.write("\n")
output_command = "mv " +sample_id+ "_scAligned_stampyFiltered.mappedSortedFixed.bam " +intermediates_folder
output_script.write(output_command)
output_script.write("\n")
output_command = "mv " +sample_id+ "_scAligned_stampyFiltered.mapped.bam " +intermediates_folder
output_script.write(output_command)
output_script.write("\n")
output_command = "mv " +sample_id+ "_scAligned_r1.sai " +intermediates_folder
output_script.write(output_command)
output_script.write("\n")
output_command = "mv " +sample_id+ "_scAligned_r2.sai " +intermediates_folder
output_script.write(output_command)
output_script.write("\n")
# get the mapping stats since bwa doesnt provide
output_command = "bamtools stats -in " +aligned_folder+ "/" +sample_id+ "_scAligned_stampyFiltered.mappedSorted.bam > " +aligned_folder+ "/" +sample_id+ "_scAligned_stampyFiltered.mappedStats.txt"
output_script.write(output_command)
output_script.write("\n\n")

# close the script file
output_script.close()

# make script executable
os.system("chmod +x " +script_name)
