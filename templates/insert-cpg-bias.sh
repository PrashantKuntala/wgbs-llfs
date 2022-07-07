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

# additional QC Rscripts
dir_qc="/dsgmnt/llfs_methyl/work/pkuntala/genome/hg38/wgbs_qc"
bg_chr1_1kb=${dir_qc}/CpGs.hg38_chr1_1kb_win.bg.gz
rscript_insert=${dir_qc}/density_insert_length.R
rscript_cpgbias=${dir_qc}/CpGbias_1kb.R

# input
bam_dedup="../${workorder}/3_bismark/${sample}_bismark_bt2_pe.deduplicated.bam"

# output
outdir="../${workorder}/6_insert_cpg_bias"
mkdir -p $outdir

bam_tmp=${outdir}/${sample}.tmp.bam
insert_tmp=${outdir}/${sample}.tmp.insert.txt.gz
cov_chr1=${outdir}/${sample}.tmp.CpG.cov_chr1_1kb_win.txt.gz

out_insert=${outdir}/${sample}.insert_length.txt
rlog_insert=${outdir}/${sample}.density_insert_length.R.log
rlog_cpgbias=${outdir}/${sample}.CpGbias_1kb.R.log


# COMMANDS
THREADS=4

# move into the output directory
cd $outdir

echo "Started insert-cpg-bias on : $(date)"

echo "1. Started samtools view $(date)"
# select the first 100K alignments
samtools view -h -@ $THREADS $bam_dedup |
    head -100000197 | samtools view -b -o $bam_tmp -@ $THREADS

echo "2. Started bamtobed $(date)"
# make temporary insert length txt file
# Write BAM alignments in BEDPE format | awk extract columns values 1,2,6, & (6-2) | write compressed version to the tempfile.
bedtools bamtobed -bedpe -i $bam_tmp | awk -vOFS="\t" '{print $1,$2,$6,$6-$2}' | gzip -nc > $insert_tmp

echo "3. Started bamtobed coverage $(date)"
# select only chr1
bedtools bamtobed -bedpe -i $bam_tmp | awk '$1=="chr1"' | bedtools coverage -counts -a $bg_chr1_1kb -b stdin | awk '$4>0 && $5>0' | gzip -nc > $cov_chr1

echo "4. Started Rscript insert $(date)"
# makes a plot for the insert sizes and writes out density
/dsgmnt/dev/prashant/R-4.2.0/R_BUILD/bin/Rscript $rscript_insert $insert_tmp $out_insert &> $rlog_insert

echo "5. Started Rscript cpgbias $(date)"
# Does pearson and spearman correlation, plots a graph.
/dsgmnt/dev/prashant/R-4.2.0/R_BUILD/bin/Rscript $rscript_cpgbias $cov_chr1 ${sample} &> $rlog_cpgbias

# This is calculating average insert size for the sample.
insert_tmp2=${outdir}/${sample}.tmp.insert.txt.gz
out_insert2=${outdir}/${sample}.avg.insert.txt
python ${dir_qc}/average_insert.py $insert_tmp2 $out_insert2

echo "Finished on : $(date)"

# remove temporary files
#rm $bam_tmp $insert_tmp $cov_chr1