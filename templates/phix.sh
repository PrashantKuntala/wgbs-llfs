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
echo "Started phiX on : $(date)"

outdir="../${workorder}/2_phix"
mkdir -p $outdir

fq1="../${workorder}/1_trimgalore/${INPUT1}"
fq2="../${workorder}/1_trimgalore/${INPUT2}"

genome_path="/dsgmnt/llfs_methyl/work/pkuntala/genome/phiX174/bwa_index/phiX174.fa"

CPU=4

# total number of reads
total=$( zcat $fq1 | grep -c "^@" )

bam=${outdir}/${sample}.bwa_phiX.bam
output=${outdir}/${sample}.bwa_phiX.txt

# Mapping to the phiX genome to check the contamination, reporting only mapped reads
bwa mem -t 4 $genome_path $fq1 $fq2 | samtools view -b -o $bam -F 4 -@ $CPU

# Count reads mapped to the phiX genome
phi=$( samtools view -c -@ $CPU $bam )

# Calculate the rate of phiX mapped reads
phi_rate=$( awk -v c=$phi -v t=$total 'BEGIN{ rate=sprintf("%.7g", c/(t*2)); print rate}' )

# Print out the results
echo -e "${sample},$total,$phi,$phi_rate" > $output

echo "Finished phiX on : $(date)"
