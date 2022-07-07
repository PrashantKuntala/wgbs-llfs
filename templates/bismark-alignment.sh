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

genome_dir="/dsgmnt/llfs_methyl/work/pkuntala/genome/hg38/bismark"
min_insert=0
max_insert=2000
threads=2

# output folder
outdir="../${workorder}/3_bismark"
tmp_dir="../${workorder}/bismark_temp"
mkdir -p $outdir
mkdir -p $tmp_dir

# Input read files
read1_fq="../${workorder}/1_trimgalore/${INPUT1}"
read2_fq="../${workorder}/1_trimgalore/${INPUT2}"

# 1. bismark PE files
log_bismark_pe=${outdir}/${sample}.bismark.log

### COMMANDS
# Mapping with bismark/bowtie2
echo "-- 1. Mapping to reference with bismark/bowtie2 : $(date)"
bismark -q -I $min_insert -X $max_insert --parallel 2 \
        --bowtie2 -N 1 -L 28 --score_min L,0,-0.6 -p $threads \
        -o $outdir --temp_dir $tmp_dir --gzip --nucleotide_coverage \
        $genome_dir -1 $read1_fq -2 $read2_fq &>$log_bismark_pe

mv $outdir/${sample}_R1_bismark_bt2_pe.bam $outdir/${sample}_bismark_bt2_pe.bam
mv $outdir/${sample}_R1_bismark_bt2_pe.nucleotide_stats.txt $outdir/${sample}_bismark_bt2_pe.nucleotide_stats.txt
mv $outdir/${sample}_R1_bismark_bt2_PE_report.txt $outdir/${sample}_bismark_bt2_PE_report.txt

echo "-- Finished on : $(date)"

