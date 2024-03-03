#! /bin/bash -login
#SBATCH -p med2                                 # partition requested
#SBATCH -J index                                # job name
#SBATCH -t 0-5:00:00                            # requested wall time (D-HH:MM)
#SBATCH -N 1                                    # number of nodes           
#SBATCH -n 1                                    # number of cores
#SBATCH -c 4                                    # number of cpus per task
#SBATCH --mem=15gb                              # memory pool to all cores
#SBATCH -e slurm.index.bwameth.j%j.err          # STANDARD ERROR FILE TO WRITE TO
#SBATCH -o slurm.index.bwameth.j%j.out          # STANDARD OUTPUT FILE TO WRITE TO
#SBATCH --mail-user=msleeper@ucdavis.edu        # YOUR EMAIL ADDRESS
#SBATCH --mail-type=ALL                         # NOTIFICATIONS OF SLURM JOB STATUS 

'''
------------------------
description and usage
------------------------
0_index-ref-genome-bwameth.sh
Created: 7/25/23 MS

this script will index a reference genome for alignment with bwa-meth

input reference genome in fasta.gz format:
  - example input: 
      - GRCh38_no_alt_analysis_set.fa.gz

output is a set of index files:
  - example outputs:
      - GRCh38_no_alt_analysis_set.fa.gz.bwameth.c2t
      - GRCh38_no_alt_analysis_set.fa.gz.bwameth.c2t.bwt
      - GRCh38_no_alt_analysis_set.fa.gz.bwameth.c2t.sa
      - GRCh38_no_alt_analysis_set.fa.gz.bwameth.c2t.pac
      - GRCh38_no_alt_analysis_set.fa.gz.bwameth.c2t.amb
      - GRCh38_no_alt_analysis_set.fa.gz.bwameth.c2t.ann

options:
  -r indicates the reference genome to use (use full path to reference genome file)
  - does not take positional arguments

usage:
sbatch 0_index-ref-genome-bwameth.sh -r [reference-genome]

resources used:
indexing human hg38 analysis set (GRCh38_no_alt_analysis_set.fa.gz)
  - time: 3 hours
  - memory: 9 gb
-------------------------------------------
'''

# activate conda in general
. "/home/msleeper/miniconda3/etc/profile.d/conda.sh"

# activate a specific conda environment, if you so choose
conda activate bwa

# make things fail on errors
set -o nounset
set -o errexit
set -x

# defining flags for arguments passed in
while getopts ':o:w:r:' flag
do
    case "${flag}" in
        r) reference=${OPTARG}
            ;;
        \?) echo "$0: Error: Invalid option: -${OPTARG}" >&2; exit 1
            ;;
        :) echo "$0: Error: option -${OPTARG} requires an argument" >&2; exit 1
            ;;
    esac
done

# Index specified file as reference genome for bwa-meth use
bwameth.py index $reference  #Indexes with BWA-MEM (default)

# Print out values of select the current jobs SLURM environment variables
env | grep SLURM

