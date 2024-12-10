#! /bin/bash -login
#SBATCH -p med2                                 # partition requested
#SBATCH -J merge                                # job name
#SBATCH -t 0-10:00:00                            # requested wall time (D-HH:MM)
#SBATCH -N 1                                    # number of nodes           
#SBATCH -n 1                                    # number of cores
#SBATCH -c 4                                    # number of cpus per task
#SBATCH --mem=10gb                              # memory pool to all cores
#SBATCH -e slurm.merge.j%j.err                  # STANDARD ERROR FILE TO WRITE TO
#SBATCH -o slurm.merge.j%j.out                  # STANDARD OUTPUT FILE TO WRITE TO
#SBATCH --mail-user=msleeper@ucdavis.edu        # YOUR EMAIL ADDRESS
#SBATCH --mail-type=ALL                         # NOTIFICATIONS OF SLURM JOB STATUS 
 
# OBSERVED RESOURCE USAGE: 3 gb and 3 hours to merge 16 files ~6-10 gb each
#  - 12 files at 15-29 gb: 4.5 hrs, 5 gb

#-----------------------------------------------------------------------------------#
# description and usage
#-----------------------------------------------------------------------------------#
# 4_merge.sh
# Created: 3/15/23 MS
#
# input accession numbers for aligned bam files 
#   - this should be the file name without .bam extension output from bwameth-workflow.sh
# outputs: 
#   merged bam of all sorted and dup-marked bams - name will be designated by -o argument
#
# takes args:
#   -w indicates where to navigate to for access to file
#   -o indicated output file name without extension for file type
#   non-option args should be the accession numbers for files being processed 
#
# sbatch 4_merge.sh -w [directory/path] -o [output-filename] [accession-num-1] [accession-num-2] [accession-num-3]...
#
#-----------------------------------------------------------------------------------#
 
# activate conda 
. "/home/msleeper/miniconda3/etc/profile.d/conda.sh"
 
# activate a specific conda environment
conda activate bamba
 
# make things fail on errors
set -o nounset
set -o errexit
set -x

# defining flags for arguments passed in
while getopts ':o:w:' flag

do
    case "${flag}" in
        o) output=${OPTARG}
           ;;
        w) where=${OPTARG}
           ;;
        \?) echo "$0: Error: Invalid option: -${OPTARG}" >&2; exit 1
           ;;
        :) echo "$0: Error: option -${OPTARG} requires an argument" >&2; exit 1
           ;;
    esac
done

# navigate to directory specified by -w flag for where
echo "Going to directory: $where"
cd $where

# needed to be able to use positional vaiables without flags
shift "$((OPTIND - 1))" # now the  positional variables have the non-option arguments

# creating an array to keep track of marked bam files that are created
FILES=()

# check timing implementation for errors by creating control values
/usr/bin/time -o ~/benchmarking/00_control/control.txt --append -f "%C,%E,%P,%K,%x" sleep 22s
/usr/bin/time -o ~/benchmarking/00_control/control.txt --append -f "%C,%E,%P,%K,%x" sleep 55s

# iterate over files to convert, sort, and mark duplicates
for i; do

    # creating array for files to use
    FILES+=($i'-marked.bam')

done

# merging all number-marked.bams
echo "Merging ${FILES[@]} into $output.bam"
/usr/bin/time -o ~/benchmarking/4_merge/merge.txt --append -f "%C,%E,%P,%K,%x" sambamba merge $output.bam ${FILES[@]}

# Print out values of select the current jobs SLURM environment variables
env | grep SLURM


