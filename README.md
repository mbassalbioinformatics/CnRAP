# CnRAP (Cut & Run Analysis Pipeline)
Analytical pipeline developed to anlayze Cut and Run data.  Inspired by both Henikoff (SEACR) and Orkin (Cut&amp;RunTools) lab pipelines. 

_**CnRAP requires both Python2 and Python3. Python2 has reached end of life and is no longer developed or maintained. A new version of CnRAP is in development to improve upon every aspect of its usage, including dependencies and compatibility. This repository is provided for compatibility with our [publication](https://www.cell.com/cell-reports/pdf/S2211-1247(20)31563-1.pdf) and its corresponding STAR Protocol.**_ 

The python scripts are written to generate bash scripts (tested on Ubuntu 16.04). This was done to allow tweaking of bash scripts later if needed without having to re-run most steps etc... This approach works better for my workflow so thats why I did it this way.

The tools required to be used in this pipeline include:
- trimmomatic - v0.36 tested
- kseq trimmer (kseq_test) - developed by the Orkin lab - tool can be found on their [bitbucket page](https://bitbucket.org/qzhudfci/cutruntools/src/master/)
- bwa - v0.7.17-r1188 tested
- [Stampy](https://www.well.ox.ac.uk/research/research-groups/lunter-group/lunter-group/stampy) - v1.0.32 tested
- bamtools - v2.5.1 tested
- samtools - v1.5 tested
- deeptools - v2.5.7 tested
- bedGraphToBigWig - v4 tested
- SEACR - v1.1 tested
- HOMER - v4.10 tested
- MEME - v5.0.5 tested
- Picard - v2.21.2 tested
- python3 - v3.6.1 tested
- python2 - v2.7+ tested
- R - v3.6.1 tested
- ChIPseeker - v1.20.0 tested

_**CnRAP is written in Python3 and so all scripts are to be called using Python3. Stampy only supports Python2. Prior to running CnRAP, make sure you have working versions of both Python3 and Python2 setup in your system path or use Anaconda to manage the different python versions. If setting up CnRAP as per the following instructions, simply typing "python" will run Python2, and you would need to type "python3" to trigger Python3. CnRAP was written around this assumption (that python2 is the default python installation in the path).**_

# Installation
1. Download and install [Anaconda 3](https://www.anaconda.com/products/individual) on your system. Install conda accepting defaults unless you are a more experienced user. Once complete, restart the shell to initialize the conda installation.
```
wget https://repo.anaconda.com/archive/Anaconda3-2021.05-Linux-x86_64.sh

bash Anaconda3-2021.05-Linux-x86_64.sh
```


2. Most of CnRAP's tools can be setup through conda directly. To avoid version incompatibilities, it’s best to install most of CnRAP’s tools in their own conda environment. This can be done by opening up the terminal and typing the command:
```
conda create -n cnrap_env python=2
```
This will create a new python2 environment named _“cnrap_env”_ for which to continue the setup. When prompted, accept the listing of the additional components needed to be installed by pressing _“Y”_.


3. Activate the created environment and install the tools available from conda.
```
conda activate cnrap_env

conda install -c bioconda trimmomatic bwa picard samtools deeptools seacr homer meme ucsc-bedgraphtobigwig bioconductor-biocinstaller bioconductor-chipseeker
```
This will go through and setup most of the tools required.


4.  Next, install the required genome for HOMER. This is done by running the following command
```
perl ~/anaconda3/pkgs/homer-<v#>/share/homer/configureHomer.pl -install <req_genome>
```
 where _"<v#>"_ is the installed version number of HOMER and _"<req_genome>"_ is the required genome required. 
  
 At the time of writing this documentation, using the latest available packages from conda, the above command would look something along the lines of...
```
perl ~/anaconda3/pkgs/homer-4.11-pl526hc9558a2_3/share/homer/configureHomer.pl -install hg38
perl ~/anaconda3/pkgs/homer-4.11-pl526hc9558a2_3/share/homer/configureHomer.pl -install sacCer3
```
User your native file-browser or tab-auto-complete on the terminal to figure out which version of HOMER is installed.
_For the [SALL4 publication](https://www.cell.com/cell-reports/pdf/S2211-1247(20)31563-1.pdf), both the hg38 and the sacCer3 genomes are needed._
  

  To get a feel for what genomes are available for HOMER to install, type the following [command](http://homer.ucsd.edu/homer/introduction/configure.html)
```
perl ~/anaconda3/pkgs/homer-<v#>/share/homer/configureHomer.pl -list
```


5. Next, you need to install Stampy manually by going to [the Stampy website](https://www.well.ox.ac.uk/research/research-groups/lunter-group/lunter-group/stampy) where you need to [register and download](https://www.well.ox.ac.uk/forms/software-download-registration) Stampy. The [Stampy Documentation](https://www.rdm.ox.ac.uk/files/research/lunter-group/stampyreadme.txt) described its installation and usage instructions. You will then need to make sure that the Stampy path is in your ~/.bashrc file by including a statement along the lines of the following wherein _"<path_to_stampy>"_ is the full path to the installation directory of Stampy.
```
export PATH=$PATH:<path_to_stampy>/stampy-1.0.32
```
_The python2 restriction on CnRAP is due to Stampy. As such, future releases of CnRAP (in-development) will not use Stampy as its no longer maintained and updated by its developers._


6. The final step in the installation of CnRAP is to setup _kseq_trimmer_ from [Cut&Run Tools](https://bitbucket.org/qzhudfci/cutruntools/src/master/). Navigate to the bitbucket repository, and download the following files to their own folder:
```
  kseq.h
  kseq_test.c
```
Next, you need to compile the files for your system by navigating to the folder were you downloaded the 2 files and at the terminal typing...
```
gcc -O2 kseq_test.c -lz -o kseq_test
```
_(Optional)_ You can now add the path to _kseq_test_ to your ~/.bashrc file or, more conveniently, just pass the full path to the CnRAP script.
  

## Genome Prep Prior to Running CnRAP

In order to use this version of CnRAP, you need to create genome indexes for your genomes of interest using both bwa and Stampy. In our case that was the masked genomes for hg38 and sacCer3 which can be downloaded from [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html). Once your genomes of interest are downloaded, build the genome indexes using the following command...
```
bwa index -p <index_name> <masked_genome_fasta>.fa
stampy.py --species=<species> --assembly=<assembly_name> -G <index_name> <masked_genome_fasta>.fa
stampy.py -g <index_name> -H <index_name>
```
... in our case our commands were...
```
bwa index -p human_bwa_index_hg38_masked hg38_masked.fa
stampy.py --species=human --assembly=hg38_masked -G human_bwa_index_hg38_masked ./hg38_masked.fa
stampy.py -g human_bwa_index_hg38_masked -H human_bwa_index_hg38_masked
```
In our case, these indexes needed to be built for both hg38 and sacCer3.

  
## Script 01 Configuration
Prior to running script 01 for each pair of fq files (R1 and R2 fq files per sample), there are a number of parameters that need to be configured internal to the script - these are found in lines 13-24. Below is a description of what each line requires as input. The script itself also contains similar information.

**Script 01 Variable** | **Required Input**
-------------------|---------------
**trimmomatic_path_jar** | full path to the trimmomatic jar file _OR_ if installed as outlined above, simply write "trimmomatic"
**adapter_path** | full path to the bbmap adapters.fa file (this can be found in ~/anaconda3/pkgs/bbmap-<v#>/opt/bbmap-<v#>/resources/)
**kseq_path_exec** | full path to the kseq_test application _OR_ if included in .bashrc path, simply write "kseq_test"
**bwa_path_exec** | full path to bwa _OR_ if installed as outlined above, simply write "bwa"
**bwa_hg_masked_index** | full path to the *folder* of the masked hg index
**bwa_sc_masked_index** | full path to the *folder* of the masked sacCer index
**stampy_hg_stidx_hash** | full path to the *folder* of the stampy hg index
**stampy_sc_stidx_hash** | full path to the *folder* of the stampy sacCer index
**samtools_path_exec** | full path to samtools _OR_ if installed as outlined above, simply write "samtools"
**picard_path_jar** | full path to the picard jar file _OR_ if installed as outlined above, simply write "picard"
**read_length_file** | full path to a text file which only contains a single number - the length of the reads eg 42
**multimap_limit** | as C&R reads are shorter in length, the likelihood of multi-mapping increases. This value simply specifies out many locations a read is allowed to map to. For unique alignments only, simple enter "1", otherwise, enter the number required.


# How to run?
Each script is very well annotated so hopefully they are easy enough to follow and you can follow my logic in what I did. 

- At the start of each python script, there are a number of arguments that need to be passed into the script when calling. I hope the arguments are self-explanatory. 
- _*Python script 01 needs to be run <b>seperately</b> for every pair of fastq files.*_  This is because each run for a pair of sequencing files can be quite time-consuming in case you want to run 1 or all of the pairs at the same time - depending on if your setup can take it - how many cores to throw at it etc...  A simple bash script can be written to run them sequentially or in parallel - again depending on your setup. 
Something that reared its ugly head when trying to run the individual scripts in parallel was sometimes 1 of the scripts would fail with a broken pipe error in python.  Not sure whether that was python3 calling it quits because of pushing the system or something with Ubuntu failed to pipe things properly - dunno tbh. If that happens just re-run the script and it should work fine.

- The remaining scripts (both python and R) then assume that all the output files are together in a single directory for each step (follow the directory structure from the comments.  Again, doing it this way gives the flexibility to run 1 or all parts depending on what you want.

Oh 1 last thing...and this is quite bad of me... these scripts have NO error handling really. No checking of input conditions or files or anything.  Big no no(!!) I know but still...  So just make sure you feed in what you want carefully into the scripts.


- **06/05/21 Update -**
Ok soooo in developing our ChIP-Seq Analysis Pipeline (ChIP-AP) (https://github.com/JSuryatenggara/ChIP-AP), a number of people have asked about Cut&Run and its use with that suite of tools. So keeping this interest in mind, we are developing CnRAPv2 which takes many design elements from ChIP-AP and tailors the analysis experience for C&R that builds on this pipeline but greatly improves upon it. The code in this repo will remain as is (with only bug fixes being made) owing to its use in our [publication](https://www.cell.com/cell-reports/pdf/S2211-1247(20)31563-1.pdf) and a companion STAR Protocol (work-in-progress). However, when CnRAPv2 is ready, will link it here for future publications and uses.
