#! /bin/bash -login
#SBATCH -p med2                                 # partition requested
#SBATCH -J dmr                                  # job name
#SBATCH -t 1-0:00:00                            # requested wall time (D-HH:MM)
#SBATCH -N 1                                    # number of nodes           
#SBATCH -n 1                                    # number of cores
#SBATCH -c 4                                    # number of cpus per task
#SBATCH --mem=15gb                              # memory pool to all cores
#SBATCH -e slurm.dmr.j%j.err                    # STANDARD ERROR FILE TO WRITE TO
#SBATCH -o slurm.dmr.j%j.out                    # STANDARD OUTPUT FILE TO WRITE TO
#SBATCH --mail-user=msleeper@ucdavis.edu        # YOUR EMAIL ADDRESS
#SBATCH --mail-type=ALL                         # NOTIFICATIONS OF SLURM JOB STATUS 

# OBSERVED RESOURCE USAGE: 

#-----------------------------------------------------------------------------------#
# description and usage (needs to be updated)
#-----------------------------------------------------------------------------------#
# 6_segment-find-dmrs.sh
# Created: 3/30/23 MS
# Modified: 7/25/23 MS # added timing element
# Modified: 8/3/23 MS # added flags for arguments
#
# this script will take pat and beta files produced by wgbstools-bam2pat.sh
# and process with wgbstools to:
#   - segment into similarly methylated chunks
#   - send segmented chunk info to tables in csv
#   - produce index file for visualizing methylation
#   - look for DMRs and send to csv
#
# usage: bash 6_segment-find-dmrs.sh -w /path/to/working/directory -g /path/to/groups/file -o /path/to/output/directory
#
# create a script that takes arguments to create the groups file
#-----------------------------------------------------------------------------------#

# activate conda in general
. "/home/msleeper/miniconda3/etc/profile.d/conda.sh"

# activate a specific conda environment, if you so choose
conda activate dmr

# make things fail on errors
set -o nounset
set -o errexit
set -x

# defining flags for arguments passed in
while getopts ':w:g:o:' flag
do
    case "${flag}" in
        w) where=${OPTARG}
            ;;
        g) groups=${OPTARG}
            ;;
        o) out=${OPTARG}
            ;;
        \?) echo "$0: Error: Invalid option: -${OPTARG}" >&2; exit 1
            ;;
        :) echo "$0: Error: option -${OPTARG} requires an argument" >&2; exit 1
            ;;
    esac
done

# # go to a particular directory
# cd /home/msleeper/programs/wgbs_tools

# # add wgbstools to path
# export PATH=${PATH}:$PWD

# navigate to directory specified by -w flag for where
echo "Going to directory: $where"
cd $where

# create project output directory
mkdir -p $out

# needed to be able to use positional vaiables without flags
shift "$((OPTIND - 1))" # now the  positional variables have the non-option arguments

# # check timing implementation for errors by creating control values
# /usr/bin/time -o ~/benchmarking/00_control/control.txt --append -f "%C,%E,%P,%K,%x" sleep 22s
# /usr/bin/time -o ~/benchmarking/00_control/control.txt --append -f "%C,%E,%P,%K,%x" sleep 55s

# segmenting into homogenously methylated chunks of CpG sites according to the tutorial settings ( removed --min_cpg 4)
# /usr/bin/time -o ~/benchmarking/6_segment-find-dmrs/segment.txt --append -f "%C,%E,%P,%K,%x" 
wgbstools segment --betas betas/*beta --max_bp 5000 -o $out/blocks.bed

# compress bed file and generate corresponding index file for visualization steps
# /usr/bin/time -o ~/benchmarking/6_segment-find-dmrs/index.txt --append -f "%C,%E,%P,%K,%x" 
wgbstools index blocks.bed

# collapse beta files to homogenously methylated chunks found in segmentation step and writes to csv
# /usr/bin/time -o ~/benchmarking/6_segment-find-dmrs/beta2table.txt --append -f "%C,%E,%P,%K,%x" 
wgbstools beta_to_table $out/blocks.bed.gz -c 10 --betas betas/*beta | column -t > $out/meth-segments.csv

# find DMRs for all 3 samples (lung, pancreas, and colon) according to tutorial settings
# note that pval is set to 1 and change this after testing phase
# /usr/bin/time -o ~/benchmarking/6_segment-find-dmrs/find-dmrs.txt --append -f "%C,%E,%P,%K,%x" 
wgbstools find_markers --blocks_path $out/blocks.bed.gz -c 25 --min_cpg 5 --min_bp 10 --max_bp 1500 --betas *beta --groups_file groups --pval 1 --targets cancer -o $out

# writes markers to csv file
cat $out/Markers.*.bed > $out/DMR_markers.csv

# Print out values of the current jobs SLURM environment variables
env | grep SLURM
