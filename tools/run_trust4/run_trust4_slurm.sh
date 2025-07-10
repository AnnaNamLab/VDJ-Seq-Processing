#!/bin/bash
#SBATCH --partition=scu-gpu            # Replace with your cluster's partition
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=trust4_run          # You can override with --job-name=<your_job_name>
#SBATCH --time=48:00:00                # Max run time (HH:MM:SS)
#SBATCH --mem=64G                      # Requested memory

# Activate the environment
conda activate /path/to/your/conda/env

# ----------- Input arguments --------------
SAMPLE=$1
READ1=$2
READ2=$3
OUTPUT_PREFIX=$4

# Validate input
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <sample_name> <read1.fastq.gz> <read2.fastq.gz> <output_prefix>"
  exit 1
fi

# ----------- TRUST4 paths -----------------
TRUST4_EXEC="run-trust4"  # assumes it's in PATH or change to full path
REF_FASTA="/path/to/hg38_bcrtcr.fa"
VJ_REF="/path/to/human_IMGT+C.fa"
BC_WHITELIST="/path/to/737K-august-2016.txt"

# ----------- Create working directory -----
cd /path/to/your/working/directory/$SAMPLE

# ----------- Run TRUST4 -------------------
$TRUST4_EXEC \
  -f "$REF_FASTA" \
  --ref "$VJ_REF" \
  -u "$READ2" \
  --barcode "$READ1" \
  --readFormat bc:0:15 \
  --barcodeWhitelist "$BC_WHITELIST" \
  --UMI "$READ1" \
  --readFormat um:16:25 \
  --od "${OUTPUT_PREFIX}_${SAMPLE}_T4_Output_$(date +%Y_%m_%d)" \
  -t 8 \
  --repseq
