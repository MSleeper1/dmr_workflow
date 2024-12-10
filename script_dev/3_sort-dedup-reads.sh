#! /bin/bash -login
#SBATCH -p med2                                 # partition requested
#SBATCH -J bamba                                # job name
#SBATCH -t 0-10:00:00                            # requested wall time (D-HH:MM)
#SBATCH -N 1                                    # number of nodes           
#SBATCH -n 1                                    # number of cores
#SBATCH -c 4                                    # number of cpus per task
#SBATCH --mem=10gb                              # memory pool to all cores
#SBATCH -e slurm.bamba.j%j.err                  # STANDARD ERROR FILE TO WRITE TO
#SBATCH -o slurm.bamba.j%j.out                  # STANDARD OUTPUT FILE TO WRITE TO
#SBATCH --mail-user=msleeper@ucdavis.edu        # YOUR EMAIL ADDRESS
#SBATCH --mail-type=ALL                         # NOTIFICATIONS OF SLURM JOB STATUS 
 
# OBSERVED RESOURCE USAGE: 2-3 gb and ~1 hour per 10 gb bam file processed
#                          3 gb and 3 hours (max and rounded up) for size 20-44 gb bam files

#-----------------------------------------------------------------------------------#
# description and usage
#-----------------------------------------------------------------------------------#
# 3_sort-dedup-reads.sh
# Created: 3/15/23 MS
# Modified: 3/25/23 MS
#
# input accession numbers for aligned bam files 
#   - this should be the file name without .bam extension output from bwameth-workflow.sh
# outputs: 
#   sorted bam - name will be accession-number-sorted.bam
#   sorted and dup-marked bam - name will be accession-number-marked.bam
# 
# options:
#   -w indicates where to navigate to for access to file
#   non-option args should be the accession numbers for files being processed 
#
# sbatch 3_sort-dedup-reads.sh -w [directory/path] [accession-num-1] [accession-num-2] [accession-num-3]...
# 
#-----------------------------------------------------------------------------------#
 
# activate conda 
. "/home/msleeper/miniconda3/etc/profile.d/conda.sh"
 
# activate a specific conda environment
conda activate bamba
 
# make things fail on errors
set -o nounset  # same as 'set -u' which treats unset variables as an error when substituting
set -o errexit  # same as 'set -e' which exits immediately if a command exits with a non-zero status
set -x          # same as set '-o xtrace' which prints commands and their arguments as they are executed

# defining flags for arguments passed in
while getopts ':w:' flag

do
    case "${flag}" in
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
    
    # for sort -o designates file output
    echo "Sorting $i.bam"
    /usr/bin/time -o ~/benchmarking/3_sort-dedup-reads/bamba_sort.txt --append -f "%C,%E,%P,%K,%x" sambamba sort -o $i-sorted.bam $i.bam 
    
    # marking duplicates in each bam file
    echo "Marking duplicates in $i.bam to output $i-marked.bam"
    /usr/bin/time -o ~/benchmarking/3_sort-dedup-reads/bamba_markdup.txt --append -f "%C,%E,%P,%K,%x" sambamba markdup $i-sorted.bam $i-marked.bam
    FILES+=($i'-marked.bam')

done

# letting the user know where to find the sorted and marked files
echo "all files have been sorted and marked"
echo "${FILES[@]} can be found in $where"

# Print out values of select the current jobs SLURM environment variables
env | grep SLURM


