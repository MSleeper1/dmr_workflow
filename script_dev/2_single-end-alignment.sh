#! /bin/bash -login
#SBATCH -p med2                                 # partition requested
#SBATCH -J align                                # job name
#SBATCH -t 1-0:00:00                            # requested wall time (D-HH:MM)
#SBATCH -N 1                                    # number of nodes           
#SBATCH -n 1                                    # number of cores
#SBATCH -c 4                                    # number of cpus per task
#SBATCH --mem=30gb                              # memory pool to all cores
#SBATCH -e slurm.alignment.j%j.err              # STANDARD ERROR FILE TO WRITE TO
#SBATCH -o slurm.alignment.j%j.out              # STANDARD OUTPUT FILE TO WRITE TO
#SBATCH --mail-user=msleeper@ucdavis.edu        # YOUR EMAIL ADDRESS
#SBATCH --mail-type=ALL                         # NOTIFICATIONS OF SLURM JOB STATUS 

# OBSERVED RESOURCE USAGE:

#-----------------------------------------------------------------------------------#
# 2_single-end-alignment.sh
# Created: 3/20/23 MS
# Last Updated: 8/20/23 MS - added time component for benchmarking and removed absolute path to bwameth.py
#
# this script will align single read fastq files to a reference and convert resulting sam file to bam file
#
# input accession numbers for aligned bam files 
#   - this should be the file name without extension output from sratools-workflow.sh
# outputs: 
#   alignment files in sam and bam format
#
# options:
#   -o indicates outfolder that will be created at location -w
#      output bam and sam files will be found here
#   -w indicates where to navigate to for access to file
#   -r indicates the reference genome to use
#   non-option args should be the accession numbers for fastq files being processed 
#
# sbatch 2_single-end-alignment.sh -w [directory/path] -o [folder-name] -r [reference-genome] [accession-num-1] [accession-num-2] [accession-num-3]...
#
#-----------------------------------------------------------------------------------#

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
        o) outfolder=${OPTARG}
            ;;
        w) where=${OPTARG}
            ;;
        r) reference=${OPTARG}
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

# create directory with name specified by -o flag for outfolder
if [ ! -d "$outfolder" ]; then
    echo "The directory $outfolder does not exits"
    echo "Creating driectory named: $outfolder"
    mkdir $outfolder
fi

# needed to be able to use positional vaiables without flags
shift "$((OPTIND - 1))" # now the  positional variables have the non-option arguments

# creating an array to keep track of marked bam files that are created
FILES=()

# check timing implementation for errors by creating control values
/usr/bin/time -o ~/benchmarking/00_control/control.txt --append -f "%C,%E,%P,%K,%x" sleep 22s
/usr/bin/time -o ~/benchmarking/00_control/control.txt --append -f "%C,%E,%P,%K,%x" sleep 55s

# iterate over files to align and convert sam to bam
for i; do

    # Align
    echo "aligning $i.fastq to $reference"

    # added time component and removed absolute path to bwameth.py
    /usr/bin/time -o ~/benchmarking/2_single-end-alignment/bwa_align.txt --append -f "%C,%E,%P,%K,%x" bwameth.py --reference $reference  $i.fastq > $outfolder/$i.sam


    # converting sam to bam
    echo "converting $i.sam to $i.bam"
    samtools view -S -b $outfolder/$i.sam > $outfolder/$i.bam
    # /usr/bin/time -o ~/benchmarking/2_single-end-alignment/samtools_convert.txt --append -f "%C,%E,%P,%K,%x" samtools view -S -b $outfolder/$i.sam > $outfolder/$i.bam


    # add completed accession number to FILES array
    FILES+=($i)
    echo "alignment and conversion complete for $i"

done

# echo array of completed accession numbers 
echo "${FILES[@]} have been properly aligned"
echo "sam and bam files can be found in $where/$outfolder"

# Print out values of select the current jobs SLURM environment variables
env | grep SLURM

