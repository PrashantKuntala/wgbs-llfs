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


# 5. bismark2bedgraph files
log_bismark2bg=${outdir}/${sample}.bismark2bedGraph.log
bedGraph=${sample}.bedGraph.gz
cg_pe=${outdir}/CpG_context_${sample}_bismark_bt2_pe.deduplicated.txt.gz
ch_pe=${outdir}/Non_CpG_context_${sample}_bismark_bt2_pe.deduplicated.txt.gz
merged=${outdir}/${sample}_bismark_bt2.extracted.txt.gz


### COMMANDS
# Generate bedGraph file
echo "-- 5. Generate bedGraph started : $(date)"
mv $ch_pe $merged
cat $cg_pe >>$merged
bismark2bedGraph --dir $outdir --cutoff 1 --CX_context --buffer_size=50G --scaffolds \
                 -o $bedGraph $merged &>$log_bismark2bg

#rm ${outdir}/$bedGraph # $merged

echo "-- Finished on : $(date)"

