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

# output folder
outdir="../${workorder}/3_bismark"
tmp_dir="../${workorder}/bismark_temp"
mkdir -p $outdir
mkdir -p $tmp_dir

# 6. coverage2cytosine files
log_cov2c=${outdir}/${sample}.coverage2cytosine.log
cov=${outdir}/${sample}.bismark.cov.gz
cx_report=${sample}.CX_report.txt.gz
cx_me=${outdir}/${sample}_bismark_bt2.CXme.txt

### COMMANDS

# Calculate average methylation levels per each CN context
echo "-- 6. Generate cytosine methylation file started : $(date)"
coverage2cytosine -o ${sample} --dir $outdir --genome_folder $genome_dir --CX_context --gzip $cov &>$log_cov2c

echo "done : $(date)"

#rm ${outdir}/$cov
zcat ${outdir}/$cx_report |
    awk 'BEGIN{ca=0;cc=0;cg=0;ct=0;mca=0;mcc=0;mcg=0;mct=0}
         $7~/^CA/ {ca+=$5; mca+=$4}
         $7~/^CC/ {cc+=$5; mcc+=$4}
         $7~/^CG/ {cg+=$5; mcg+=$4}
         $7~/^CT/{ct+=$5; mct+=$4}
         END{printf("CA\t%d\t%d\t%.3f\n", ca, mca, mca/(ca+mca));
             printf("CC\t%d\t%d\t%.3f\n", cc, mcc, mcc/(cc+mcc));
             printf("CG\t%d\t%d\t%.3f\n", cg, mcg, mcg/(cg+mcg));
             printf("CT\t%d\t%d\t%.3f\n", ct, mct, mct/(ct+mct));}' >$cx_me
echo ""

# Print the files generated in descending order of file sizes.
echo "-- The results..."
ls -lSh ${outdir}/*${sample}*
echo ""
echo "-- Finished on : $(date)"

