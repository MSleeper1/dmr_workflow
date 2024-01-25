#! /bin/bash -login
#SBATCH -p med2                                 # partition requested
#SBATCH -J index                                # job name
#SBATCH -t 0-5:00:00                            # requested wall time (D-HH:MM)
#SBATCH -N 1                                    # number of nodes           
#SBATCH -n 1                                    # number of cores
#SBATCH -c 4                                    # number of cpus per task
#SBATCH --mem=18000                             # memory pool to all cores
#SBATCH -e slurm.index.bismark.j%j.err          # STANDARD ERROR FILE TO WRITE TO
#SBATCH -o slurm.index.bismark.j%j.out          # STANDARD OUTPUT FILE TO WRITE TO
#SBATCH --mail-user=msleeper@ucdavis.edu        # YOUR EMAIL ADDRESS
#SBATCH --mail-type=ALL                         # NOTIFICATIONS OF SLURM JOB STATUS 


#-----------------------------------------------------------------------------------#
# description and usage
#-----------------------------------------------------------------------------------#
# 0_index-ref-genome.sh
# Created: 1/24/24 MS
#
# this script will index a reference genome for alignment with bismark
#
# input reference genome
#
# options:
#   -r indicates the directory containing reference genome to use (use full path to directory)
# does not take positional arguments
#
# sbatch 0_index-ref-genome-bismark.sh -r [reference-genome]
#-----------------------------------------------------------------------------------#

# activate conda in general
. "/home/msleeper/miniconda3/etc/profile.d/conda.sh"

# activate a specific conda environment, if you so choose
conda activate bismark

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

# Index specified file as reference genome for bismarkuse
bismark_genome_preparation --verbose --bowtie2 --parallel 6 $reference

# Print out values of select the current jobs SLURM environment variables
env | grep SLURM

