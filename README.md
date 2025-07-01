
# VDJ processing

VDJ (Variable, Diversity, and Joining) processing is a key mechanism in the generation of diverse B cell receptors (BCRs), which are essential for the adaptive immune response. The BCR is an antibody anchored on the surface of B cells, recognizing specific antigens through highly variable regions, generated during B cell development. The generation of a diverse B cell receptor (BCR) repertoire is fundamental to the immune system's ability to recognize and respond to a vast array of pathogens. This repository contains all codes related to generating BCR/TCR sequencing, and all downstream analysis.

## Requirements

1) TRUST4 (https://github.com/liulab-dfci/TRUST4)
2) igblast (https://ncbi.github.io/igblast/cook/How-to-set-up.html)
3) PhyloPart_v2.1

## Method

Building the phylogenetic tree from BCR data has 5 main steps: generating BCR repertoire using TRUST4, removing the doublets, detecting the dominant clone, correcting the VDJ assignments errors in Hodgkin and Reed-Sternberg (HRS) cells, and building the phylogenetic tree

<div align="center">
  <img width="545" alt="BCR_steps" src="https://github.com/user-attachments/assets/2d5769f9-1764-46a4-9117-1ef7365244b8" />
</div>

### 1)  Generating BCR Repertoire Using TRUST4

To reconstruct the BCR repertoire, we have applied TRUST4. TRUST4 is an open-source algorithm for reconstruction of immune receptor repertoire in αβ/γδ B cells and T cells from bulk and single-cell RNA-seq data [1]. TRUST4 can assemble full sequences based on VDJ assignments. This algorithm extracts the good TCR/BCR candidate reads from BAM or FATSQ files with fast speed and has a de-novo approach for assembling the candidate reads. It implements a greedy strategy by aligning the candidate reads to an existing contig. For contigs annotation, TRUST4 uses ImmunoGeneTics (IMGT) database [2]. TRUST4 has multiple output files including the extracted contigs and corresponding nucleotide weight, a fasta file for annotation of the consensus assembly, and files with the focus on the CDR3 and constructed full sequences. Also, TRUST4 provides an extra step for clustering the contigs based on similarity of VDJs and CDR3, and during this step they remove any partial or incomplete CDR3. This clustering file is the output we used for downstream. Trust4 paper is in https://www.nature.com/articles/s41592-021-01142-2

### 2)	Removing the Doublets

To remove the doublets based on BCR data, we have used the CDR3 raw file from TRUST4. This file contains all the contigs with number of reads that are supporting that. Each cell can be identified as singlet, or not-assigned, or doublet. The histogram of read counts for contigs per each cell defines the status of each cell. For each cell: if a cell has only one contig, with any number of reads, that cell will be called singlet. If not, then using excess kurtosis, normalized entropy, and confidence level of VDJ assignment, a cell will be called doublet, not-assigned, or singlet. If there are  equal or less than 7 reads in cells with 2 or more contigs, then the cell can be called singlet if either of these conditions are met: either the normalized entropy is small (less than 0.8) or the CDR3 Levenshtein distance in contigs is less than 0.2. If there are equal or less than 7 reads in cells with 2 or more contigs, but none of the conditions for CDR3 similarity or small normalized-entorpy are met, then the cell is called not-assigned.
If the excess kurtosis is positive value (not +Inf), then the cell is called singlet. If the excess kurtosis is negative value (not -Inf), then the cell is called doublet. If the kurtosis in +Inf or -Inf, and the number of reads is more than 7, then if the normalized entropy is more than 0.8, the cell is called doublet.

There is an extra step to check the doublets based on the conidence levelof VDJ assignments. Do obtain the the confidence level of VDJ assignments, the length of CDR1 and CDR2 are interfered. If one of CDR1 or CDR2 were missing or both of them are incomplete, then that contig's VDJ assignment has medium-confident. If there is complete CDR1 and CDR2, then the contig's VDJ assignment is high_confident. If both CDR1 and CDR2 are missing then the contig's VDJ assignment is low_confident.

For any doublet, we have checked the confidence level of VDJ assignment to mae sure the doublets are coming from different origin. To do this step, for any doublet, the confidence level of VDJ-assignments for two contigs (with highest number of reads) are checked. If both contigs are high-confident, then the cell-status stays doublet. If the two contigs are both low-confident, or one of them is low-confident and another is medium-confident, then we calculate the similarity of the contigs CDR3 based on Levenshtein distance . If the CDR3s were similar (Levenshtein distance <= 0.2), then the cell status will change to singlet. But if the CDR3s were not similar, then the cell status stays doublet.

The singlets need to be intersected with HRS cells or B cells and we do this step per each chain (IGH and IGL/IGK). Then the clustering output from step 1 should be subsetted based on obtained singlets at this step. 


### 3)	Detecting the Dominant Clone

After removing the doublets from step 2, we identify the dominant clone per each chain. The dominant clone is the clone with the majority of HRS cells with same V/D/J assignment. We put no threshold for similarity of CDR3, since that causes to remove some of the HRS cells from the cancer clone, while we know that HRS cells are part of the clone. The dominant V/D/J on each chain will be recorded for the next step. 

### 4)	Correcting the V/D/J Assignments Errors

Form the output of the step 3 for detecting the dominant V/D/J, it happens that some of the HRS cells are carrying a different V/D/J comparing to dominant clone. One possible reason is about the errors during V/D/J gene identification after alignment of contigs to the sequences from IMGT database. To correct these errors, we extracted the contigs of HRS cells, and using the fasta file, which is one of the outcomes of TRUST4, started to extract new V/D/J assignments. From IgBlast, the top 3 V/D/J genes are extracted. Then for every HRS cell, if obtained IgBlast V/D/Js are the same of dominant V/D/J, then we correct the error from TRUST4. We, also explore the fasta file per each of the contigs, and if any of the top 3 V/D/Js are among the detected V/D/J for that contig, we correct that. That means per each contig any of the top 3 V/D/J obtained from IgBlast are neighbors of the dominant clone and that can still show that the cell belongs to the dominant clone. After correcting the V/D/Js per chain, we re-extracted all the contigs with similar V/D/J as the dominant clone and took it for generating the phylogenetic tree. 

### 5)	Building the Phylogenetic Tree

After extracting the HRS cells with similar V/D/J to the dominant clone and adding other HRS cells with not-similar V/D/J to the dominant clone, we use them to generate the phylogenetic tree. Phylogenetic tree is a diagram that shows the evolutionary relationships between CDR3/full-sequences of BCR. We used a maximum likelihood phylogenetic tree from Ape package in R [3]. The tree aims to put the sequences with highest similarity in the same clade based on genetic distances that are calculated using a substitution model. We used F81 substitution model for calculating the distances between nucleotides, since it is one of the wildly been used models for generating phylogenetic tree [4]. We generated the phylogenetic tree based on CDR3 and then annotated the cell states on each node on the tree. Also, we can identify the V/D/J assignments on the tree, or other BCR parameters to investigate their evolutionary process.



