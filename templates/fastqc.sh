#!/bin/sh

####################
##### PBS PARAMS
####################
# Add a custom jobname
#$ -N $jobname

# pass all environment variables to the job
#$ -V

# Make sure that the .e and .o file arrive in the working directory
#$ -cwd

#Merge the standard out and standard error to one file
#$ -j y

####################
##### TOOL PARAMS
####################
echo "Started fastqc on: $(date)" 

THREADS=4

# absolute path to the folder containing your fastq files
indir="/dsgmnt/llfs_methyl/work/pkuntala/other/mini-test"

# output folder
outdir="../${workorder}/0_fastqc"
mkdir -p $outdir

fastqc -o $outdir --noextract --nogroup -t $THREADS ${indir}/$INPUT1 ${indir}/$INPUT2

echo "Finished fastqc on : $(date)"