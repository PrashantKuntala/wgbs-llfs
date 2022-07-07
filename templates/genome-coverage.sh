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

echo "Started genome coverage on : $(date)"

# INPUT
cx_me="../${workorder}/3_bismark/${sample}_bismark_bt2.CXme.txt"

cnt_c=$(( 598683433+600854940 ))	# Watson strand + Crick strand
cnt_c=$(( $cnt_c - 171823*2 )) 		# Discard chrEBV
cnt_cg=$(( 29303965 * 2 ))


# OUTPUT
outdir="../${workorder}/7_coverage"
cov_genome=${outdir}/${sample}.genome_cov.txt


# COMMANDS
mkdir -p $outdir
c_cov=$( cat $cx_me | awk -F"\t" -v c=$cnt_c 'BEGIN{s=0} {s+=$2+$3} END{print s/c}' )
cg_cov=$( cat $cx_me | awk -F"\t" -v c=$cnt_cg 'BEGIN{s=0} $1=="CG" {s+=$2+$3} END{print s/c}' )

echo -e "${sample},$c_cov,$cg_cov" > $cov_genome

echo "Finished genome coverage on : $(date)"