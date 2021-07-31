#!/bin/bash 

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --nodes=1
#SBATCH --job-name=read2_hist
#SBATCH --error=read2_hist.err
#SBATCH --output=read2_hist.out
#SBATCH --time=0-10:00:00
#SBATCH --cpus-per-task=8

conda activate bgmp_py39
/projects/bgmp/wrosales/bioinfo/Bi622/Demultiplex/Assignment-the-first/qual_score_dist_read2.py