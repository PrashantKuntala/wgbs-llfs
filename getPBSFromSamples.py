#!/bin/python

"""
Program to generate PBS shell scripts for each tool in the WGBS pipeline using templates and a samplelist.

example:
python getPBSFromSamples.py templates/trimgalore.sh samples.txt demo

"""

import argparse
import os
from string import Template

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('template',help='Specify a template file')
    parser.add_argument('samples', help='List of filenames in a txt or csv format')
    parser.add_argument('workorder', help='Workorder or folder name from MGI that is being processed')
    args = parser.parse_args()

    # read PBS template
    with open(args.template, 'r') as f:
        template_src = Template(f.read())

    job_prefix = args.template.strip().split("/")[1].strip(".sh")

    # create output folder for pbs_scripts
    outfolder = "./{}_pbs_scripts".format(job_prefix)
    if not os.path.exists(outfolder):
        os.mkdir(outfolder)
        print("created : {}".format(outfolder))

    sampleList = open(args.samples,'r').readlines()

    for sample in sampleList:
        temp = sample.strip().split(',')[0]
        print(temp)

        # substitute dictionary for fastqc & trimgalore vs others.
        # trimgalore uses (or automatically renames) the output extension to .fq.gz 
        if job_prefix.endswith('fastqc') or job_prefix.endswith('trimgalore'):
            template_data = {
            'jobname': '{}-{}'.format(job_prefix,temp),
            'sample': '{}'.format(temp),
            'INPUT1':'{}_R1.fastq.gz'.format(temp),
            'INPUT2':'{}_R2.fastq.gz'.format(temp),
            'workorder': args.workorder
            }
        else:
            template_data = {
            'jobname': '{}-{}'.format(job_prefix,temp),
            'sample': '{}'.format(temp),
            'INPUT1':'{}_R1.fq.gz'.format(temp),
            'INPUT2':'{}_R2.fq.gz'.format(temp),
            'workorder': args.workorder
            }        

        # code that subsitutes values that are required.
        filename = outfolder+"/"+template_data['jobname']+".sh"
        result = template_src.safe_substitute(template_data)
        with open(filename,'w') as fq:
            fq.write(result)
