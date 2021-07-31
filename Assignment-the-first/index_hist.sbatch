#!/bin/bash 

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --nodes=1
#SBATCH --job-name=index_hist
#SBATCH --error=index_hist.err
#SBATCH --output=index_hist.out
#SBATCH --time=0-10:00:00
#SBATCH --cpus-per-task=8

/projects/bgmp/wrosales/bioinfo/Bi622/Demultiplex/Assignment-the-first/qual_score_dist_index.py