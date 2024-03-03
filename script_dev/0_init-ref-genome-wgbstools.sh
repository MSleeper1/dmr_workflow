#! /bin/bash -login
#SBATCH -p med2                                 # partition requested
#SBATCH -J wgbstools-init                       # job name
#SBATCH -t 0-0:30:00                            # requested wall time (D-HH:MM)
#SBATCH -N 1                                    # number of nodes           
#SBATCH -n 1                                    # number of cores
#SBATCH -c 4                                    # number of cpus per task
#SBATCH --mem=3000                               # memory pool to all cores
#SBATCH -e slurm.init.wgbstools.j%j.err          # STANDARD ERROR FILE TO WRITE TO
#SBATCH -o slurm.init.wgbstools.j%j.out          # STANDARD OUTPUT FILE TO WRITE TO
#SBATCH --mail-user=msleeper@ucdavis.edu        # YOUR EMAIL ADDRESS
#SBATCH --mail-type=ALL                         # NOTIFICATIONS OF SLURM JOB STATUS 

'''
------------------------
description and usage
------------------------
0_init-ref-genome-wgbstools.sh
Created: 1/25/24 MS

this script will initialize a reference genome for use with wgbstools

input: 
  path to reference genome file (can be fasta.gz or fasta) and name to use for reference genome

output:
  directory of files created by wgbstools located in the wgbstools program references ".../wgbstools/references/[name-of-reference-genome]"

required arguments:
  -r indicates the reference genome to use (use full path to reference genome fasta file)
  -n indicates the name of the reference genome (ex: "hg38"). This name will be the name of the directory created by wgbstools.
  - does not take positional arguments

usage:
  sbatch 0_init-ref-genome-wgbstools.sh -r [reference-genome-file-path] -n [name-of-reference-genome]

resources used:
  indexing human hg38 analysis set (GRCh38_no_alt_analysis_set.fa.gz)
    - time: 5 minutes
    - memory: 2 gb

-------------------------------------------
'''

# activate conda in general
. "/home/msleeper/miniconda3/etc/profile.d/conda.sh"

# activate a specific conda environment, if you so choose
conda activate dmr # see wgbstools.yaml for specific environment details

# make things fail on errors
set -o nounset
set -o errexit
set -x

# defining flags for arguments passed in
while getopts ':o:w:r:n:' flag
do
    case "${flag}" in
        r) reference=${OPTARG}
            ;;
        n) genome_name=${OPTARG}
            ;;
        \?) echo "$0: Error: Invalid option: -${OPTARG}" >&2; exit 1
            ;;
        :) echo "$0: Error: option -${OPTARG} requires an argument" >&2; exit 1
            ;;
    esac
done


# check if reference genome and name of reference genome were provided
if [ -z "$reference" ]
then
    echo "Error: no reference genome provided. please provide a reference genome with the -r flag followed by the path to the reference genome file"
    exit 1
fi

if [ -z "$genome_name" ]
then
    echo "Error: no name for reference genome provided. please provide a name for the reference genome with the -n flag followed by the name of the reference genome"
    exit 1
fi

# if refernence ends in .gz, unzip it
if [[ $reference == *.gz ]]
then
    gunzip $reference
    reference=${reference%.gz}
fi

# -f will overwrite any existing reference genome with the same name
wgbstools init_genome -f $genome_name --fasta_path $reference 

# Print out values of select the current jobs SLURM environment variables
env | grep SLURM

