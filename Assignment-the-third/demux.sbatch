#!/bin/bash 

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --nodes=1
#SBATCH --job-name=demux
#SBATCH --error=demux.err
#SBATCH --output=demux.out
#SBATCH --time=0-12:00:00
#SBATCH --cpus-per-task=8

conda activate bgmp_py39
/usr/bin/time -v /projects/bgmp/wrosales/bioinfo/Bi622/Demultiplex/Assignment-the-third/demux.py -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -bc /projects/bgmp/shared/2017_sequencing/indexes.txt -o /projects/bgmp/wrosales/bioinfo/Bi622/Demultiplex/Assignment-the-third/demux_out/