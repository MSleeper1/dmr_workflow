#! /bin/bash -login
#SBATCH -p med2                                 # partition requested
#SBATCH -J bam2pat                              # job name
#SBATCH -t 1-0:00:00                            # requested wall time (D-HH:MM)
#SBATCH -N 1                                    # number of nodes           
#SBATCH -n 1                                    # number of cores
#SBATCH -c 4                                    # number of cpus per task
#SBATCH --mem=20gb                              # memory pool to all cores
#SBATCH -e slurm.bam2pat.j%j.err                # STANDARD ERROR FILE TO WRITE TO
#SBATCH -o slurm.bam2pat.j%j.out                # STANDARD OUTPUT FILE TO WRITE TO
#SBATCH --mail-user=msleeper@ucdavis.edu        # YOUR EMAIL ADDRESS
#SBATCH --mail-type=ALL                         # NOTIFICATIONS OF SLURM JOB STATUS 

# OBSERVED RESOURCE USAGE: used 5-6 gb and 14-15 hours to process 3 merged bams that are around 100 gb each
#                          out of memory error when using 10 gb mem for 207 gb bam file

#-----------------------------------------------------------------------------------#
# description and usage
#-----------------------------------------------------------------------------------#
# 5_convert-bams.sh
# Created: 3/27/23 MS
# Last Updated: 8/28/23 MS - added flags for ref genome, temp dir, and mbias plot generation. added timing implementation
# Future Updates: add flags for qc
# this script will take bam files and process with wgbstools to:
#   - generate pat and beta files for bam files
#   - generate mbias plots for bam files
#
# input bam files
# outputs pat and beta files
#
# options:
#   -o indicates outfolder that will be created at location -w
#      output bam and sam files will be found here
#   -w indicates where to navigate to for access to file
#   -r indicates reference genome name to use (must be in name from wgbstools genome list) 
#   -t indicates temp directory to use
# non-option args are bam file names with .bam extension
#
# sbatch 5_convert-bams.sh -w [directory/path] -o [folder-name] -t [temp_dir] -r [ref_name] [bam file] [bam file]...
# resource usage: used 5-6 gb and 14-15 hours to process 3 merged bams that are around 100 gb each
#-----------------------------------------------------------------------------------#

# activate conda in general
. "/home/msleeper/miniconda3/etc/profile.d/conda.sh"

# activate a specific conda environment, if you so choose
conda activate cfDNA_env

# go to a particular directory
cd /home/msleeper/programs/wgbs_tools

# add wgbstools to path
export PATH=${PATH}:$PWD

# make things fail on errors
set -o nounset
set -o errexit
set -x

# defining flags for arguments passed in
while getopts ':o:w:r:t:' flag
do
    case "${flag}" in
	    o) outfolder=${OPTARG}
            ;;
        w) where=${OPTARG}
            ;;
        r) ref_name=${OPTARG}
            ;;
        t) temp_dir=${OPTARG}
            ;;
        \?) echo "$0: Error: Invalid option: -${OPTARG}" >&2; exit 1
            ;;
        :) echo "$0: Error: option -${OPTARG} requires an argument" >&2; exit 1
            ;;
    esac
done

# set region of genome to investigate and translate genomic loci to CpG-index range
# change region to an input option after trial run
#region=chr3:119527929-119531943
#wgbstools convert -r $region

# navigate to directory specified by -w flag for where
echo "Going to directory: $where"
cd $where

# create directory with name specified by -o flag for outfolder if it does not already exits
if [[ -d "$outfolder" ]] ;
then
    echo "Directory $outfolder exists"
else
    echo "The directory $outfolder does not exist"
    echo "Creating driectory named: $outfolder"
    mkdir $outfolder
fi

# create directory with name specified by -t flag for temp_dir if it does not already exits
if [[ -t "$temp_dir" ]] ;
then
    echo "Directory $temp_dir exists"
else
    echo "The directory $temp_dir does not exist"
    echo "Creating driectory named: $temp_dir"
    mkdir $temp_dir
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

    # generate pat and beta files for all bam files in the directory bams (add -f to force overwrite)
    /usr/bin/time -o ~/benchmarking/5_convert-bams/bam2pat.txt --append -f "%C,%E,%P,%K,%x" wgbstools bam2pat -f --out_dir $outfolder --mbias --genome $ref_name --temp_dir $temp_dir $i 

    # add completed accession number to FILES array
    FILES+=($i)
    echo "bam2pat complete for $i"

done

# echo array of completed accession numbers
echo "${FILES[@]} have been processed with bam2pat"
echo "output files can be found $where/$outfolder"

# Print out values of the current jobs SLURM environment variables
env | grep SLURM
