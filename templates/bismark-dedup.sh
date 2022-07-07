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

# output folder
outdir="../${workorder}/3_bismark"
tmp_dir="../${workorder}/bismark_temp"
mkdir -p $outdir
mkdir -p $tmp_dir

# 2. deduplicate files
bam_pe=${outdir}/${sample}_bismark_bt2_pe.bam
log_dedup_pe=${outdir}/${sample}.deduplicate_bismark.log


### COMMANDS
# Dedpulicate reads
echo "-- 2. deduplicate_bismark started : $(date)"
deduplicate_bismark -p --bam $bam_pe --output_dir $outdir &>$log_dedup_pe
echo "-- Finished on : $(date)"
