
# Single-Cell BCR Sequencing Analysis

VDJ (Variable, Diversity, and Joining) processing is a key mechanism in the generation of diverse B cell receptors (BCRs), which are essential for the adaptive immune response. The BCR is an antibody anchored on the surface of B cells, recognizing specific antigens through highly variable regions, generated during B cell development. The generation of a diverse B cell receptor (BCR) repertoire is fundamental to the immune system's ability to recognize and respond to a vast array of pathogens. This repository contains all codes related to generating BCR/TCR sequencing, and all downstream analysis.

## Requirements

1) TRUST4 (https://github.com/liulab-dfci/TRUST4) [ref]
2) igblast (https://ncbi.github.io/igblast/cook/How-to-set-up.html) [ref]
3) PhyloPart_v2.1 (https://sourceforge.net/projects/phylopart/) [ref]

## Method

Building the phylogenetic tree from BCR data has 5 main steps: generating BCR repertoire using TRUST4, removing the doublets, detecting the dominant clone, correcting the VDJ assignments errors in Hodgkin and Reed-Sternberg (HRS) cells, and building the phylogenetic tree.

### BCR processing steps:

1. Generating BCR Repertoire Using TRUST4
1. Removing the Doublets
1. Detecting the Dominant Clone
1. Correcting the V/D/J Assignments Errors
1. Building the Phylogenetic Tree

<div align="left">
  <img width="500" height="280" alt="BCR_steps" src="https://github.com/user-attachments/assets/2d5769f9-1764-46a4-9117-1ef7365244b8" />
</div>

More details about the steps are provided in the workflow folder.



### ▶️ Run TRUST4 via SLURM

You can submit TRUST4 jobs to your HPC cluster using the generalized SLURM batch script.
Before submitting the job, make sure the script is executable:

```bash
chmod +x run_trust4_slurm.sh
```
Then run run_trust4_slurm.sh
```bash
sbatch run_trust4_slurm.sh <sample_name> <read1.fastq.gz> <read2.fastq.gz> <output_prefix>
```
The cdr3.out file from TRUST4 is then applied as the input of trust-cluster.py for clustering similar CDR3 sequences.
To do this, first move to the directory where you have you copied trust-cluster.py, and provide cdr3.out as input.
The output is an tsv file.

```bash
python trust-cluster.py <input_cdr3_file> > <output_file>

```
The outcome file (cluster_clone.tsv) is used for more downstream analysis and generating the BCR phylogenetic tree.

### ▶️ How to Run the doublet finding Script

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
For TCR: 
```bash
Rscript TCRdoublet_finding_after_refinement_May2025.R \
  --input /path/to/trust4_output_dir \
  --metadata /path/to/cell_metadata.csv \
  --sample HL1 \
  --file HL1_cdr3.out

```
### ▶️ Detecting the dominant clone 
To detect the cancer clone the heatmap of the CDR3 sequences helps to detcte the most dominant V,D,J.
With ethis code: 
```bash
Rscript heatmap_before_VDJcorrection_March2024.R
--sample HL1 \
--metadata /path/to/cell_metadata.csv \
--input /path/to/trust4_output_dir \
--cluster_file /path/to/cluster_clone.tsv \
--light_doublet /path/to/light chain doublets \
--heavy_doublet /path/to/heavy chain doublets \

```
An example if run in on local terminal:
```bash
Rscript heatmap_before_VDJcorrection_March2024.R \
  --sample HL1 \
  --metadata /Users/saramoein/Documents/new_run_HL_May2025/2024-11-26_CellMetadata_HL1-24incHL8R_RetainedCellsOnly_MainCellTypeAndSubtypeNames.csv \
  --input /Users/saramoein/Documents/new_run_HL_May2025/HL1 \
  --cluster_file out_clone_sara.tsv \
  --light_doublet /Users/saramoein/Documents/new_run_HL_May2025/FINAL_doublets_BCR_thre07_ent08/HL1_Rawdata_IGK_IGL.csv \
  --heavy_doublet /Users/saramoein/Documents/new_run_HL_May2025/FINAL_doublets_BCR_thre07_ent08/HL1_Rawdata_IGH.csv

```
### ▶️ Correcting the V/D/J Assignments Errors (optional)

This step uses the igblast tool and fasta file as extra resources to survive the HRS cells with different V,D,J from the dominant clone.
To run this part, a list of HRS contigs are required, based on this pattern "contigs_${sample}_${chain}.txt".

For example, we subset the HRS contigs in "contigs_HL1_IGL.txt". IGL is  the clonal chain.

```bash
AAACGGGCAAAGTGCG_21556
AAACGGGCAGTCAGAG_21980
AAACGGGGTATCAGTC_21656
AAAGATGGTTGGTTTG_22094
AAAGCAAAGTGGTAGC_21648
AAATGCCTCTCTGTCG_21695
AACGTTGGTAAGTTCC_21525
AACTGGTCAAACGTGG_22089
AACTGGTTCTGCTGTC_21948

```
This file should be saved to the directory that *annot.fa (output from trust4) is located.
Then we should run submit_igblast.sh:

```bash
chmod +x submit_igblast.sh
./submit_igblast.sh <sample> <clone_chain> <path to *annot.fa file as an output of TRUST4>
```

Here is an example: 
```bash
./submit_igblast.sh HL1 IGL /athena/namlab/scratch/sam4032/HL1_s1s2/HL1_T4_Output_2024_01_04
```
After running igblast, FIX_VDJ_BCR_step1_step2_igblast.R is run. 
After running this code, the new heatmap after potential correction of VDJs can be generated.

### ▶️ Heatmap after VDJ correction 
If you have run the step for correcting the VDJ assignment, then you will have a heatmap based on the corrected VDJs. For that, run:  
```bash
Rscript heatmap_after_VDJcorrection_March2024.R
--sample HL1 \
--metadata /path/to/cell_metadata.csv \
--input /path/to/VDJ correction output.csv \
--DominantChain dominant chian that contains the clone \
--light_doublet /path/to/light chain doublets.csv \
--heavy_doublet /path/to/heavy chain doublets.csv \
--output /path/to/output_folder
```
An example if run in on local terminal:
```bash
Rscript heatmap_after_VDJcorrection_March2024.R \
  --sample HL1 \
  --metadata /Users/saramoein/Documents/new_run_HL_May2025/2024-11-26_CellMetadata_HL1-24incHL8R_RetainedCellsOnly_MainCellTypeAndSubtypeNames.csv \
  --input /Users/saramoein/Documents/new_run_HL_May2025/HL1/HL1_FILTERED_out_clone_dominantChain_HL1.csv' \
  --DominantChain IGL \
  --light_doublet /Users/saramoein/Documents/new_run_HL_May2025/FINAL_doublets_BCR_thre07_ent08/HL1_Rawdata_IGK_IGL.csv \
  --heavy_doublet /Users/saramoein/Documents/new_run_HL_May2025/FINAL_doublets_BCR_thre07_ent08/HL1_Rawdata_IGH.csv \
  --output /Users/saramoein/Documents/new_run_HL_May2025/HL1
```
### ▶️ Generating phylogenetic tree
To generate the phylogenetic tree
