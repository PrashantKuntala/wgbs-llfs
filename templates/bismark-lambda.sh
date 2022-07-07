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

genome_dir="/dsgmnt/llfs_methyl/work/pkuntala/genome/lambda"
min_insert=0
max_insert=2000
threads=2

# output folder
outdir="../${workorder}/4_bismark_lambda"
tmp_dir="../${workorder}/bismark_temp"
mkdir -p $outdir
mkdir -p $tmp_dir

# Input read files
read1_fq="../${workorder}/1_trimgalore/${INPUT1}"
read2_fq="../${workorder}/1_trimgalore/${INPUT2}"

# 1. bismark PE files
log_bismark_pe=${outdir}/${sample}.bismark.log

# 2. deduplicate files
bam_pe=${outdir}/${sample}_bismark_bt2_pe.bam
log_dedup_pe=${outdir}/${sample}.deduplicate_bismark.log

# 3. extract files
log_methx_pe=${outdir}/${sample}.bismark_methylation_extractor.log
bam_dedup_pe=${outdir}/${sample}_bismark_bt2_pe.deduplicated.bam

# 4. bismark2report

# 5. bismark2bedgraph files
log_bismark2bg=${outdir}/${sample}.bismark2bedGraph.log
bedGraph=${sample}.bedGraph.gz
cg_pe=${outdir}/CpG_context_${sample}_bismark_bt2_pe.deduplicated.txt.gz
ch_pe=${outdir}/Non_CpG_context_${sample}_bismark_bt2_pe.deduplicated.txt.gz
merged=${outdir}/${sample}_bismark_bt2.extracted.txt.gz

# 6. coverage2cytosine files
log_cov2c=${outdir}/${sample}.coverage2cytosine.log
cov=${outdir}/${sample}.bismark.cov.gz
cx_report=${sample}.CX_report.txt.gz
cx_me=${outdir}/${sample}_bismark_bt2.CXme.txt

### COMMANDS
echo "-- Started on $(date)"
echo ""


# Mapping with bismark/bowtie2
echo "-- 1. Mapping to reference with bismark/bowtie2 : $(date)"
bismark -q -I $min_insert -X $max_insert --parallel 2 \
        --bowtie2 -N 1 -L 28 --score_min L,0,-0.6 -p $threads \
        -o $outdir --temp_dir $tmp_dir --gzip --nucleotide_coverage \
        $genome_dir -1 $read1_fq -2 $read2_fq &>$log_bismark_pe

mv $outdir/${sample}_R1_bismark_bt2_pe.bam $outdir/${sample}_bismark_bt2_pe.bam
mv $outdir/${sample}_R1_bismark_bt2_pe.nucleotide_stats.txt $outdir/${sample}_bismark_bt2_pe.nucleotide_stats.txt
mv $outdir/${sample}_R1_bismark_bt2_PE_report.txt $outdir/${sample}_bismark_bt2_PE_report.txt

echo "done : $(date)"

# Dedpulicate reads
echo "-- 2. deduplicate_bismark started : $(date)"
deduplicate_bismark -p --bam $bam_pe --output_dir $outdir &>$log_dedup_pe
echo "done : $(date)"

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
echo "done : $(date)"

# Generate bedGraph file
echo "-- 5. Generate bedGraph started : $(date)"
mv $ch_pe $merged
cat $cg_pe >>$merged
bismark2bedGraph --dir $outdir --cutoff 1 --CX_context --buffer_size=50G --scaffolds \
                 -o $bedGraph $merged &>$log_bismark2bg

#rm ${outdir}/$bedGraph # $merged
echo "done : $(date)"

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
