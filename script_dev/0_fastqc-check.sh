#! /bin/bash -login
#SBATCH -p med2                                 # partition requested
#SBATCH -J qc-check                             # job name
#SBATCH -t 0-8:00:00                            # requested wall time (D-HH:MM)
#SBATCH -N 1                                    # number of nodes           
#SBATCH -n 1                                    # number of cores
#SBATCH -c 4                                    # number of cpus per task
#SBATCH --mem=8gb                               # memory pool to all cores
#SBATCH -e slurm.fastqc.j%j.err                 # STANDARD ERROR FILE TO WRITE TO
#SBATCH -o slurm.fastqc.j%j.out                 # STANDARD OUTPUT FILE TO WRITE TO
#SBATCH --mail-user=msleeper@ucdavis.edu        # YOUR EMAIL ADDRESS
#SBATCH --mail-type=ALL                         # NOTIFICATIONS OF SLURM JOB STATUS 


#-----------------------------------------------------------------------------------#
# description and usage
#-----------------------------------------------------------------------------------#
# 0_fastqc-check.sh
# Created: 7/26/23 MS
#
# this script will accept fastq files and run fastqc on them
#
# input: 
#   accession numbers for fastq files
# outputs: 
#   fastqc files
#
# options:
#   -o indicates outfolder that will be created at location -w
#      output fastqc files will be found here
#   -w indicates where to navigate to for access to files
#   
# non-option args should be the accession numbers for fastq files being processed 
#
# sbatch 0_fastqc-check.sh -w [directory/path] -o [folder-name] [accession-num-1] [accession-num-2] [accession-num-3]...
#
#-----------------------------------------------------------------------------------#

# activate conda in general
. "/home/msleeper/miniconda3/etc/profile.d/conda.sh"

# activate a specific conda environment, if you so choose
conda activate fastqc

# make things fail on errors
set -o nounset
set -o errexit
set -x

# defining flags for arguments passed in
while getopts ':o:w:' flag
do
    case "${flag}" in
        o) outfolder=${OPTARG}
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

if [ ! -d "qc" ]; then
    mkdir qc
fi

# create directory with name specified by -o flag for outfolder
if [ ! -d "$outfolder" ]; then
    echo "The directory $outfolder does not exits"
    echo "Creating driectory named: $outfolder"
    mkdir $outfolder
fi

# create array to hold processed accession numbers
FILES=()

# needed to be able to use positional vaiables without flags
shift "$((OPTIND - 1))" # now the  positional variables have the non-option arguments

# check timing implementation for errors by creating control values
/usr/bin/time -o ~/benchmarking/00_control/control.txt --append -f "%C,%E,%P,%K,%x" sleep 22s
/usr/bin/time -o ~/benchmarking/00_control/control.txt --append -f "%C,%E,%P,%K,%x" sleep 55s

# iterate over files to run fastqc on
for i; do

    /usr/bin/time -o ~/benchmarking/0_fastqc-check/fastqc.txt --append -f "%C,%E,%P,%K,%x" fastqc -o qc/$outfolder $i"_1.fastq" $i"_2.fastq" # run fastqc on fastq file
    
    # add completed accession number to FILES array
    FILES+=($i)
    echo "fastqc complete for $i"

done


# echo array of completed accession numbers 
echo "${FILES[@]} have been processed with fastqc"
echo "fastqc files can be found in $where/qc/$outfolder"

# Print out values of select the current jobs SLURM environment variables
env | grep SLURM





