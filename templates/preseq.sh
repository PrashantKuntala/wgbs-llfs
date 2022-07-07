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

outdir="../${workorder}/5_preseq"
mkdir -p $outdir

bam="../${workorder}/3_bismark/${sample}_bismark_bt2_pe.bam"
bam_sorted=${outdir}/${sample}.sorted.bam
output=${outdir}/${sample}.preseq_lc_extrap.txt

### COMMANDS
THREADS=6

echo "Started on : $(date)"

echo "1. Started samtools sort: $(date)"
samtools sort -m 2G -o $bam_sorted -@ $THREADS $bam

echo "2. Started preseq: $(date)"
preseq lc_extrap -o $output -B -P -D $bam_sorted
rm $bam_sorted

echo "Finished on : $(date)"