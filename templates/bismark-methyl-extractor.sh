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

threads=2

# output folder
outdir="../${workorder}/3_bismark"
tmp_dir="../${workorder}/bismark_temp"
mkdir -p $outdir
mkdir -p $tmp_dir


# 3. extract files & # 4. bismark2report
log_methx_pe=${outdir}/${sample}.bismark_methylation_extractor.log
bam_dedup_pe=${outdir}/${sample}_bismark_bt2_pe.deduplicated.bam


### COMMANDS
# Run methylation extractor for the sample
echo "-- 3. bismark_methylation_extractor started : $(date)"
bismark_methylation_extractor --paired-end --no_overlap --comprehensive --merge_non_CpG --report \
                               --gzip --parallel $threads -o $outdir $bam_dedup_pe &>$log_methx_pe
echo "done : $(date)"

# Generate HTML Processing Report
echo "-- 4. bismark2report started : $(date)"
bismark2report -o ${sample}_bismark_bt2_PE_report.html --dir $outdir \
               --alignment_report ${outdir}/${sample}_bismark_bt2_PE_report.txt \
               --dedup_report ${outdir}/${sample}_bismark_bt2_pe.deduplication_report.txt \
               --splitting_report ${outdir}/${sample}_bismark_bt2_pe.deduplicated_splitting_report.txt \
               --mbias_report ${outdir}/${sample}_bismark_bt2_pe.deduplicated.M-bias.txt \
               --nucleotide_report ${outdir}/${sample}_bismark_bt2_pe.nucleotide_stats.txt &>${outdir}/${sample}.bismark2report.log
echo "-- Finished on : $(date)"

