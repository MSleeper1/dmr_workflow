#! /bin/bash -login
#SBATCH -p med2                                 # partition requested
#SBATCH -J sra                                  # job name
#SBATCH -t 1-0:00:00                            # requested wall time (D-HH:MM)
#SBATCH -N 1                                    # number of nodes
#SBATCH -n 1                                    # number of cores
#SBATCH -c 4                                    # number of cpus per task
#SBATCH --mem=10gb                              # memory pool to all cores
#SBATCH -e slurm.sra.j%j.err                # STANDARD ERROR FILE TO WRITE TO
#SBATCH -o slurm.sra.j%j.out                # STANDARD OUTPUT FILE TO WRITE TO
#SBATCH --mail-user=msleeper@ucdavis.edu        # YOUR EMAIL ADDRESS
#SBATCH --mail-type=ALL                         # NOTIFICATIONS OF SLURM JOB STATUS

# OBSERVED RESOURCE USAGE:
#-----------------------------------------------------------------------------------#
# 1_sra_download.sh
# Created: 2/25/23 MS
#
# Downloads sra files based on accession numbers passed in as non-flagged arguments
# Flags -w and -o are requires and must preceed accesion numbers passed in
# -w indicates where to navigate to before downloading SRA files
# -o indicates the outfolder name and is used to create a directory in the location specified by -w
# all SRA files will be downloaded within the folder -o in the location -w
# 
# sbatch 1_sra_download.sh -w [directory/path] -o [folder-name] [accession-num-1] [accession-num-2] [accession-num-3]...
# 
#-----------------------------------------------------------------------------------#

# activate conda in general
. "/home/msleeper/miniconda3/etc/profile.d/conda.sh"

# activate a specific conda environment, if you so choose
conda activate sra-download

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

# create directory with name specified by -o flag for outfolder
if [ ! -d "$outfolder" ]; then
    echo "The directory $outfolder does not exits"
    echo "Creating driectory named: $outfolder"
    mkdir $outfolder
fi

cd $outfolder

# needed to be able to use positional vaiables without flags
shift "$((OPTIND - 1))" # now the  positional variables have the non-option arguments

# check timing implementation for errors by creating control values
/usr/bin/time -o ~/benchmarking/00_control/control.txt --append -f "%C,%E,%P,%K,%x" sleep 22s
/usr/bin/time -o ~/benchmarking/00_control/control.txt --append -f "%C,%E,%P,%K,%x" sleep 55s

# fetching data and requesting files from SRA with remaining arguments passed without a flag
for i; do

    echo "prefetching:"
    echo $i
    /usr/bin/time -o ~/benchmarking/1_sra-download/prefetch.txt --append -f "%C,%E,%P,%K,%x" prefetch $i

    echo "requesting file:"
    echo $i
    /usr/bin/time -o ~/benchmarking/1_sra-download/fasterq-dump.txt --append -f "%C,%E,%P,%K,%x" fasterq-dump $i

done

# Print out values of select the current jobs SLURM environment variables
env | grep SLURM
