#!/bin/python

"""
Program to generate a batch submission script. A script that submits (N) jobs. 

example:
python createBatches.py bismark_pbs_scripts/ 2 normal
"""

import argparse
import os
from string import Template


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

header = """#!/bin/sh

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
"""

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('scripts',help='Path to the folder containing PBS scripts')
    parser.add_argument('chunks', help='Number of jobs in each batch')
    parser.add_argument('queue', help='Queue name to which the jobs are submitted')
    args = parser.parse_args()

    files = os.listdir(args.scripts)
    if '.DS_Store' in files:
        files.remove('.DS_Store')
    batches = list(chunks(files, int(args.chunks)))

    for i in range(0, len(batches)):
        outfile = 'batch-{}.sh'.format(i)
        with open(outfile,'w') as f:
            headstring = header.replace('$jobname','batch-{}'.format(i))
            f.write(headstring+"\n")
            for j in batches[i]:
                f.write('qsub -q {} {}\n'.format(args.queue,os.path.abspath(args.scripts)+"/"+j))
