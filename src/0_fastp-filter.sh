#! /bin/bash -login
#SBATCH -p med2                                 # partition requested
#SBATCH -J fastp                                # job name
#SBATCH -t 3-0:00:00                            # requested wall time (D-HH:MM)
#SBATCH -N 1                                    # number of nodes           
#SBATCH -n 1                                    # number of cores
#SBATCH -c 4                                    # number of cpus per task
#SBATCH --mem=40gb                              # memory pool to all cores
#SBATCH -e slurm.fastp.j%j.err                  # STANDARD ERROR FILE TO WRITE TO
#SBATCH -o slurm.fastp.j%j.out                  # STANDARD OUTPUT FILE TO WRITE TO
#SBATCH --mail-user=msleeper@ucdavis.edu        # YOUR EMAIL ADDRESS
#SBATCH --mail-type=ALL                         # NOTIFICATIONS OF SLURM JOB STATUS 

#-----------------------------------------------------------------------------------#
# description and usage
#-----------------------------------------------------------------------------------#
# 0_fastp-filter.sh
# Created: 8/3/23 MS
#
# this script will accept paired end fastq files and run fastp quality control and filtering on them including:
#   - length filtering (minimum length 16 bp)
#   - adapter trimming  (detected by per-read overlap analysis)
#   - quality trimming  (sliding window trimming, window size 4, minimum quality 20)
#   - deduplication (remove duplicate reads)
#   - generate json and html output files (containing quality control information)
#   - generate filtered fastq files (containing reads that passed quality control)
#
# input: 
#   accession numbers for paired end fastq files
# outputs: 
#   fastp html and json files (named [accession-num]_fastp.html and [accession-num]_fastp.json)
#   filtered fastq files (named [accession-num]_1_fastp-filtered.fastq and [accession-num]_2_fastp-filtered.fastq)
#
# options:
#   -o indicates outfolder that will be created at location -w
#      output fastqc files will be found here
#   -w indicates where to navigate to for access to files
#   
# non-option args should be the accession numbers for paired end fastq files being processed 
#
# sbatch 0_fastp-filter.sh -w [directory/path] -o [folder-name] [accession-num-1] [accession-num-2] [accession-num-3]...
#
#-----------------------------------------------------------------------------------#

# activate conda in general
. "/home/msleeper/miniconda3/etc/profile.d/conda.sh"

# activate a specific conda environment, if you so choose
conda activate fastp

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

    # run fastp on fastq file
    /usr/bin/time -o ~/benchmarking/0_fastp-filter/fastqc.txt --append -f "%C,%E,%P,%K,%x" fastp -i $i"_1.fastq" -I $i"_2.fastq" -o $i"_1_fastp-filtered.fastq" -O $i"_2_fastp-filtered.fastq" --detect_adapter_for_pe --dedup --json $outfolder/$i"_fastp.json" --html $outfolder/$i"_fastp.html" 
    
    # add completed accession number to FILES array
    FILES+=($i)
    echo "fastp complete for $i"

done


# echo array of completed accession numbers 
echo "${FILES[@]} have been processed with fastp"
echo "fastp html and json files can be found in $where/qc/$outfolder"
echo "fastp filtered fastq files can be found $where"

# Print out values of select the current jobs SLURM environment variables
env | grep SLURM
