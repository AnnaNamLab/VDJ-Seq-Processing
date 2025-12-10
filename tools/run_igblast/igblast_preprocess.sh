#!/bin/bash

#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=igblast
#SBATCH --time=48:00:00
#SBATCH --mem=128G

# -----------------------------
# Get command line arguments
# -----------------------------
sample="$1"         # e.g., HL8
chain="$2"          # e.g., IGL
working_dir="$3"    # e.g., /athena/.../HL8_T4_Output_2024_01_04; the folder that contains the *annot.fa from trust4

if [[ -z "$sample" || -z "$chain" || -z "$working_dir" ]]; then
  echo "Usage: sbatch run_igblast.sh <sample> <chain> <working_dir>"
  exit 1
fi

# -----------------------------
# Move to working directory
# -----------------------------
cd "$working_dir" || {
  echo "Error: Directory not found: $working_dir"
  exit 1
}

# -----------------------------
# Locate *_annot.fa file
# -----------------------------
annot_fa=$(ls singleline*_annot.fa 2>/dev/null | head -n 1)  ### trust4 annot.fa file

if [[ -z "$annot_fa" ]]; then
  echo "Error: No *_annot.fa file found in $working_dir"
  exit 1
fi

echo "Found FASTA file: $annot_fa"

# -----------------------------
# Clean FASTA headers
# -----------------------------
sed '/^>/ s/ .*//' "$annot_fa" > fix_my_sample.fa

# -----------------------------
# Subset contigs
# -----------------------------
contig_file="contigs_${sample}_${chain}.txt"
if [[ ! -f "$contig_file" ]]; then
  echo "Error: Contig file not found: $contig_file"
  exit 1
fi

grep -w -A 1 -f "$contig_file" fix_my_sample.fa --no-group-separator > sub_fix_my_sample_${chain}.fa

# -----------------------------
# Run IgBLAST
# -----------------------------


if [[ "$chain" == "IGL" ]]; then
    echo "Running IgBlast for IGL chain..."
    igblastn \
        -germline_db_V /home/sam4032/share/igblast/database/imgt_human_IGLV_clean.fasta \
        -num_alignments_V 3 \
        -germline_db_D /home/sam4032/share/igblast/database/empty_IGLD.fasta \
        -num_alignments_D 3 \
        -germline_db_J /home/sam4032/share/igblast/database/imgt_human_IGLJ_clean.fasta \
        -num_alignments_J 3 \
        -query sub_fix_my_sample_${chain}.fa \
        -organism human \
        -outfmt 19 \
        -out ${sample}_igblast_output2_${chain}.txt

elif [[ "$chain" == "IGK" ]]; then
    echo "Running IgBlast for IGK chain..."
    igblastn \
        -germline_db_V /home/sam4032/share/igblast/database/imgt_human_IGKV_clean.fasta \
        -num_alignments_V 3 \
        -germline_db_D /home/sam4032/share/igblast/database/empty_IGKD.fasta \
        -num_alignments_D 3 \
        -germline_db_J /home/sam4032/share/igblast/database/imgt_human_IGKJ_clean.fasta \
        -num_alignments_J 3 \
        -query sub_fix_my_sample_${chain}.fa \
        -organism human \
        -outfmt 19 \
        -out ${sample}_igblast_output2_${chain}.txt

elif [[ "$chain" == "IGH" ]]; then
    echo "Running IgBlast for IGH chain..."
    igblastn \
        -germline_db_V /home/sam4032/share/igblast/database/imgt_human_IGHV_clean.fasta \
        -num_alignments_V 3 \
        -germline_db_D /home/sam4032/share/igblast/database/imgt_human_IGHD_clean.fasta \
        -num_alignments_D 3 \
        -germline_db_J /home/sam4032/share/igblast/database/imgt_human_IGHJ_clean.fasta \
        -num_alignments_J 3 \
        -query sub_fix_my_sample_${chain}.fa \
        -organism human \
        -outfmt 19 \
        -out ${sample}_igblast_output2_${chain}_test1.txt

else
    echo "ERROR: Unsupported chain '$chain'. Must be IGH, IGK, or IGL."
    exit 1
fi




echo "IgBLAST completed: ${sample}_igblast_output2_${chain}.txt"




