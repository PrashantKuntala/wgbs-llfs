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

echo "Started trimgalore on: $(date)" 

trim_base1=10
trim_base2=15

cores=4

# input folder
indir="/dsgmnt/llfs_methyl/work/pkuntala/other/mini-test"

# output folder
outdir="../${workorder}/1_trimgalore"
mkdir -p $outdir

trim_galore -q 20 --phred33 --fastqc --fastqc_args "-o $outdir --noextract --nogroup" \
            --illumina --stringency 1 -e 0.1 --length 20 \
            --clip_R1 $trim_base1 --clip_R2 $trim_base2 \
            -o $outdir \
            -j $cores \
            --paired --retain_unpaired -r1 21 -r2 21 ${indir}/$INPUT1 ${indir}/$INPUT2

# renaming the output files (removing _val_1 and _val_2 from trimmed fastq, fastqc.html files)
mv ${outdir}/${sample}_R1_val_1.fq.gz ${outdir}/${sample}_R1.fq.gz
mv ${outdir}/${sample}_R2_val_2.fq.gz ${outdir}/${sample}_R2.fq.gz

mv ${outdir}/${sample}_R1_val_1_fastqc.html ${outdir}/${sample}_R1_fastqc.html
mv ${outdir}/${sample}_R2_val_2_fastqc.html ${outdir}/${sample}_R2_fastqc.html

echo "Finished trimgalore on : $(date)"
