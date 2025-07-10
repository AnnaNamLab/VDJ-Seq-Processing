
# VDJ processing

VDJ (Variable, Diversity, and Joining) processing is a key mechanism in the generation of diverse B cell receptors (BCRs), which are essential for the adaptive immune response. The BCR is an antibody anchored on the surface of B cells, recognizing specific antigens through highly variable regions, generated during B cell development. The generation of a diverse B cell receptor (BCR) repertoire is fundamental to the immune system's ability to recognize and respond to a vast array of pathogens. This repository contains all codes related to generating BCR/TCR sequencing, and all downstream analysis.

## Requirements

1) TRUST4 (https://github.com/liulab-dfci/TRUST4)
2) igblast (https://ncbi.github.io/igblast/cook/How-to-set-up.html)
3) PhyloPart_v2.1 (https://sourceforge.net/projects/phylopart/)

## Method

Building the phylogenetic tree from BCR data has 5 main steps: generating BCR repertoire using TRUST4, removing the doublets, detecting the dominant clone, correcting the VDJ assignments errors in Hodgkin and Reed-Sternberg (HRS) cells, and building the phylogenetic tree

<div align="left">
  <img width="500" height="280" alt="BCR_steps" src="https://github.com/user-attachments/assets/2d5769f9-1764-46a4-9117-1ef7365244b8" />
</div>

### BCR processing steps:

1. Generating BCR Repertoire Using TRUST4<br>running trust4.sh code  
1. Removing the Doublets<br>running the script to detect singlet/ doublets/ not-assigned cell status using the doublet-finding code in BCR
1. Detecting the Dominant Clone
1. Correcting the V/D/J Assignments Errors
1. Building the Phylogenetic Tree

   ## ▶️ How to Run the doublet finding Script

Make sure you have R and the required packages installed.

Then run the script from the terminal using for BCR:

```bash
Rscript doublet_finding_after_refinement_July2025_github.R \
  --input /path/to/trust4_output_dir \
  --metadata /path/to/cell_metadata.csv \
  --sample HL1 \
  --file HL1_cdr3.out

```

This command runs the BCR doublet detection using TRUST4 output and associated metadata.
