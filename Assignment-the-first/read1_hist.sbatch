#!/bin/bash 

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --nodes=1
#SBATCH --job-name=read1_hist
#SBATCH --error=read1_hist.err
#SBATCH --output=read1_hist.out
#SBATCH --time=0-10:00:00
#SBATCH --cpus-per-task=8

conda activate bgmp_py39
/projects/bgmp/wrosales/bioinfo/Bi622/Demultiplex/Assignment-the-first/qual_score_dist_read1.py