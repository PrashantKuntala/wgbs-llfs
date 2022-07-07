# wgbs-llfs
Whole Genome Bisulfite Sequencing Analysis

### **Dependencies**
- samtools & htslib (>=1.15.1)
- bedtools (v2.26.0)
- pigz (2.7)
- trimgalore (0.6.6)
- bismark (0.23.1)
  - perl GD::Graph package. It used for generating M-bias plots.
- fastqc (v0.11.9)
- bowtie2 (2.2.9)
- bwa (0.7.12-r1039)
- python (>=3.7)

## Usage
> NOTE: Before using `fastqc.sh` or `trimgalore.sh` templates, please make sure the `indir` is set correctly within the template file. It should be set to the absolute path to the folder containing your input Read1 & Read2 fastq files. Scripts assumes these filenames end with _R1.fastq.gz & _R2.fastq.gz. To test-run the pipeline on the toy example provided (see `example/data`) the existing value of `indir` will work and no modifications are required.

### **`getPBSFromSamples.py`**

Lets you create pbs scripts autmatically from a template and sample list. Useful when you have hundreds of samples on which to run the WGBS pipeline.

```bash
python getPBSFromSamples.py templates/trimgalore.sh samples.txt demo
```

 Takes three inputs: 
 - one of the pbs template scripts from `templates/` folder 
 - sample list (one sample per line, see `example/samples.txt`) for which a pbs script needs to be created and 
 - workorder or output folder-name to store the results of running various tools of the pipeline.

### **`createBatches.py`**

Lets you create batches for processing. Useful to control how many jobs you submit at a time and also not to overwhelm the HPC queues which can result in job failures due to insufficient memory.

```bash
python createBatches.py bismark_pbs_scripts/ 2 normal
```

 Takes three inputs: 
 - Absolute path to the folder containing your pbs scripts. 
 - Number of chunks or scripts to submit in a batch.
 - PBS queue name to which you want to submit these jobs.
 
 ## Tool User Guide

|Tool | Reference URL |
|-----|-----|
|`Bismark` | [user guide](https://github.com/FelixKrueger/Bismark/tree/master/Docs)|
|`TrimGalore` | [user guide](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md) |
| `preseq` | [docs](http://smithlabresearch.org/software/preseq/)|
