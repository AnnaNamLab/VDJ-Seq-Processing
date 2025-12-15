
# Single-Cell VDJ processing

The B cell receptor (BCR) is a membrane-bound antibody that recognizes specific antigens via highly variable regions formed during B cell development. This diversity, generated through V(D)J recombination, is essential for effective adaptive immunity. This repository contains all code for BCR/TCR sequence generation and downstream analysis.

## Requirements

1) TRUST4 (https://github.com/liulab-dfci/TRUST4) [1]
2) igblast (https://ncbi.github.io/igblast/cook/How-to-set-up.html) [2]
3) IgPhyML (https://github.com/immcantation/igphyml) [3, 4]

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
<div align="left">
  <a href="assets/myfile.pdf"> <!-- link to PDF -->
    <img width="500" height="280" alt="My PDF Preview" src="https://github.com/AnnaNamLab/VDJ-Seq-Processing/workflow/flowchaart_BCR.pdf?raw=true" />
  </a>
</div>

More details about the steps are provided in the workflow folder.



### Run TRUST4 via SLURM

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

### How to run the doublet finding step

Make sure you have R and the required packages installed.

Then run the script from the terminal using for BCR:

```bash
Rscript doublet_finding_BCR.R \
  --input /path/to/trust4_output_dir \
  --metadata /path/to/cell_metadata.csv \
  --sample HL1 \
  --file /path/to/HL1_cdr3.out

```

This command runs the BCR doublet detection using TRUST4 output and associated metadata.
For TCR: 
```bash
Rscript doublet_finding_TCR.R \
  --input /path/to/trust4_output_dir \
  --metadata /path/to/cell_metadata.csv \
  --sample HL1 \
  --file /path/to/HL1_cdr3.out

```

### Refining the V/D/J assignments 

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
chmod +x igblast_preprocess.sh
./igblast_preprocess.sh <sample> <clone_chain> <path to *annot.fa file as an output of TRUST4>
```

Here is an example: 
```bash
./igblast_preprocess.sh HL1 IGL /athena/namlab/scratch/sam4032/HL1_s1s2/HL1_T4_Output_2024_01_04
```
After running igblast, refining_VDJ_BCR.R is run. 



### Generating phylogenetic tree

To generate the phylogenetic tree, IgPhyML package is used. To install it, please refere to: https://igphyml.readthedocs.io/en/latest/install.html
The input of igphyml is generated from  running heatmap_after_VDJcorrection.R and  then IgPhyML_fasta_prepration.R. 
After installation, copy you CDR3_preprocessed_MAFFT_aligned.fa to your working directory and run the below script:

```bash
sbatch igphyml_run.sh HL1
```
make sure that igphyml_run.sh is modifed based on the directory where igphyml is installed.
you can interactively run the IgPhyML if the number of CDR3 sequences is low. Here is an example:

```
sample=HL1
/athena/namlab/scratch/sam4032/CDR3_tree/igphyml/src/igphyml \
 -i /athena/namlab/scratch/sam4032/CDR3_tree/singleCDR3/$sample/sampled_unique_sequences_${sample}.fasta \
 -m GY \
 --run_id singleCDR3_$sample
```
The igphyml tree downstream anlysis will be then performed by running igphyml_downStream.R. 

References
----------

1. Song L, Cohen D, Ouyang Z, Cao Y, Hu X, Liu XS. TRUST4: immune repertoire reconstruction from bulk and single-cell RNA-seq data. Nat Methods. 2021;18(6):627-630. doi:10.1038/s41592-021-01142-2.

2. Ye J, Ma N, Madden TL, Ostell JM. IgBLAST: an immunoglobulin variable domain sequence analysis tool. Nucleic Acids Res. 2013;41(W34–40). doi:10.1093/nar/gkt382.

3. Hoehn KB, Lunter G, Pybus OG. A Phylogenetic Codon Substitution Model for Antibody Lineages. Genetics. 2017;206(1):417–427. https://doi.org/10.1534/genetics.116.196303

4. Hoehn KB, Vander Heiden JA, Zhou JQ, Lunter G, Pybus OG, Kleinstein SH. Repertoire-wide phylogenetic models of B cell molecular evolution reveal evolutionary signatures of aging and vaccination. bioRxiv. 2019. https://doi.org/10.1101/558825
EOF
