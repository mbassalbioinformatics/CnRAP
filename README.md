# CnRAP (Cut & Run Analysis Pipeline)
Analytical pipeline developed to anlayze Cut and Run data.  Inspired by both Henikoff (SEACR) and Orkin (Cut&amp;RunTools) lab pipelines. 

The python scripts are written to generate bash scripts (tested on Ubuntu 16.04). This was done to allow tweaking of bash scripts later if needed without having to re-run most steps etc...  This approach works better for my workflow so thats why I did it this way.

The tools required to be used in this pipeline include (in order of use)
- trimmomatic - v0.36 tested
- kseq trimmer - developed by the Orkin lab - tool can be found on their bitbucket page (https://bitbucket.org/qzhudfci/cutruntools/src/master/)
- bwa - v0.7.17-r1188 tested
- stampy - v1.0.32 tested
- bamtools - v2.5.1 tested
- samtools - v1.5 tested
- deeptools - v2.5.7 tested
- bedGraphToBigWig - v4 tested
- SEACR - v1.1 tested
- HOMER - v4.10 tested
- MEME - v5.0.5
- python3 - v3.6.1 tested from Anaconda
- R - v3.6.1 tested
- ChIPseeker - v1.20.0 tested from Bioconductor.

If any version seems a little out of date its just because it got updated online and I didnt update it on my setup.  Do I intend to update which versions are tested? Not planning on it for now... We'll see though. I do like to keep my tools updated though.  If any changes I'll update here.

How to run?
Pretty basic really.  Each script is very well annotated so hopefully they are easy enough to follow and you can follow my logic in what I did.  At the start of each python script, there are a number of arguments that need to be passed into the script when calling. I hope the arguments are self-explanatory. Python script 01 needs to be run <b>seperately</b> for every pair of fastq files.  This is because each run for a pair can be quite time-consuming in case you want to run 1 or all of the pairs at the same time - depending on if your setup can take it - how many cores to throw at it etc...  A simple bash script can be written to run them sequentially or in parallel - again depending on your setup. Something that reared its ugly head when trying to run the individual scripts in parallel was sometimes 1 of the scripts would fail with a broken pipe error in python.  Not sure whether that was python3 calling it quits because of pushing the system or something with Ubuntu failed to pipe things properly - dunno tbh. If that happens just re-run the script and it should work fine.

The remaining scripts (both python and R) then assume that all the output files are together in a single directory for each step (follow the directory structure from the comments.  Again, doing it this way gives the flexibility to run 1 or all parts depending on what you want.

Overall though, if any issues I will try to help to the best of my ability but unless there's a need for the work in my lab to develop/add features I probably won't be doing that, sorry :(

Oh 1 last thing...and this is quite bad of me... these scripts have NO error handling really. No checking of input conditions or files or anything.  Big no no(!!) I know but still...  So just make sure you feed in what you want carefully into the scripts.
