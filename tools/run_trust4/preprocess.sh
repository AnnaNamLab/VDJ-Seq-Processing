#!/bin/bash
#SBATCH --partition=scu-gpu   # cluster-specific
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=code_slurm
#SBATCH --time=2-00:00:00   # HH/MM/SS
#SBATCH --mem=64G   # memory requested, units available: K,M,G,T

conda activate /athena/namlab/scratch/sam4032/trust
cd /athena/namlab/scratch/sam4032/HL10_s1s2


run-trust4 -f /athena/namlab/scratch/tsa4002/trust4/trust4_basefiles/hg38_bcrtcr.fa \
          --ref /athena/namlab/scratch/tsa4002/trust4/trust4_basefiles/human_IMGT+C.fa \
          -u /athena/namlab/scratch/sam4032/HL8R/TCR/all_HL8R_R2.fastq.gz \
          --barcode /athena/namlab/scratch/sam4032/HL8R/TCR/all_HL8R_R1.fastq.gz \
          --readFormat bc:0:15 \
          --barcodeWhitelist /athena/namlab/scratch/tsa4002/trust4/trust4_basefiles/737K-august-2016.txt \
          --UMI /athena/namlab/scratch/sam4032/HL8R/TCR/all_HL8R_R1.fastq.gz \
          --readFormat um:16:25 \
          --od TCR_HL8R_T4_Output_$(date +“%Y_%m_%d”) \
          -t 8 \
          --repseq